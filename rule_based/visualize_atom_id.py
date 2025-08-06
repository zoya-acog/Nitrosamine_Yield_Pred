# import pandas as pd
# from rdkit import Chem
# from rdkit.Chem.Draw import rdMolDraw2D
# from pathlib import Path
# import logging

# # Set up logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)

# def visualize_denitrosated_atoms(input_csv, output_dir):
#     """Visualize denitrosated molecules with highlighted atom IDs."""
#     logger.info(f"Visualizing denitrosated atoms from {input_csv}")
    
#     # Load the CSV
#     try:
#         df = pd.read_csv(input_csv)
#     except Exception as e:
#         logger.error(f"Error reading input CSV {input_csv}: {e}")
#         raise
    
#     # Create visualization directory
#     vis_dir = output_dir / "visualizations"
#     vis_dir.mkdir(parents=True, exist_ok=True)
    
#     # Loop and draw for a few molecules
#     for idx, row in df.iterrows():
#         smiles = row.get("Recovered_SMILES")
#         atom_id = row.get("Recovered_Atom_ID")
        
#         if pd.isna(smiles):
#             continue
        
#         mol = Chem.MolFromSmiles(smiles)
#         if mol is None:
#             logger.warning(f"Invalid SMILES: {smiles}")
#             continue
        
#         Chem.rdDepictor.Compute2DCoords(mol)
        
#         drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
#         drawer.drawOptions().addAtomIndices = True
        
#         highlight_atoms = [int(atom_id)] if pd.notna(atom_id) else []
#         drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
#         drawer.FinishDrawing()
        
#         svg = drawer.GetDrawingText().replace('svg:', '')
        
#         # Save SVG
#         name = row.get('Name', f"Row_{idx}")
#         svg_file = vis_dir / f"denitrosated_{name}_{idx}.svg"
#         with open(svg_file, 'w') as f:
#             f.write(svg)
#         logger.info(f"Saved visualization for {name} (SMILES: {smiles}, Atom ID: {atom_id}) to {svg_file}")
        
#         # Limit display (optional)
#         if idx > 9:  # Show only first 10
#             break



import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path
import logging
import base64

# Set up logging
logging.basicConfig(level=logging.INFO, filename='rule_based_pipeline.log')
logger = logging.getLogger(__name__)

def visualize_denitrosated_atoms(input_csv, output_dir):
    """Visualize denitrosated molecules with highlighted atom IDs and return base64-encoded SVGs."""
    logger.info(f"Visualizing denitrosated atoms from {input_csv}")
    
    # Convert output_dir to Path if not already
    output_dir = Path(output_dir)
    
    # Load the CSV
    try:
        df = pd.read_csv(input_csv)
        logger.info(f"Loaded CSV with {len(df)} rows")
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise
    
    # Create visualization directory
    vis_dir = output_dir / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)
    
    # List to store base64-encoded SVGs
    base64_svgs = []
    
    # Loop and draw for each molecule
    for idx, row in df.iterrows():
        name = row.get('Name', f"Row_{idx}")
        smiles = row.get("Recovered_SMILES")
        atom_id = row.get("Recovered_Atom_ID")
        
        # Validate inputs
        if pd.isna(smiles) or not smiles:
            logger.warning(f"Row {idx} ({name}): Missing or invalid Recovered_SMILES")
            base64_svgs.append(None)
            continue
        
        if pd.isna(atom_id):
            logger.warning(f"Row {idx} ({name}): Missing Recovered_Atom_ID")
            base64_svgs.append(None)
            continue
        
        try:
            # Convert atom_id to int (handle float like 1.0)
            atom_id = int(float(atom_id))
        except (ValueError, TypeError) as e:
            logger.warning(f"Row {idx} ({name}): Invalid Recovered_Atom_ID '{atom_id}': {e}")
            base64_svgs.append(None)
            continue
        
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Row {idx} ({name}): Invalid SMILES '{smiles}'")
                base64_svgs.append(None)
                continue
            
            Chem.rdDepictor.Compute2DCoords(mol)
            
            drawer = rdMolDraw2D.MolDraw2DSVG(200, 150)
            drawer.drawOptions().addAtomIndices = True
            drawer.drawOptions().padding = 0.2
            
            highlight_atoms = [atom_id] if atom_id < mol.GetNumAtoms() else []
            drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
            drawer.FinishDrawing()
            
            svg = drawer.GetDrawingText().replace('svg:', '')
            
            # Save SVG to file for debugging
            svg_file = vis_dir / f"denitrosated_{name}_{idx}.svg"
            with open(svg_file, 'w') as f:
                f.write(svg)
            logger.info(f"Row {idx} ({name}): Saved SVG to {svg_file}")
            
            # Convert SVG to base64
            base64_svg = base64.b64encode(svg.encode('utf-8')).decode('utf-8')
            base64_svgs.append(base64_svg)
            
        except Exception as e:
            logger.error(f"Row {idx} ({name}): Failed to generate SVG for SMILES '{smiles}', Atom ID '{atom_id}': {e}")
            base64_svgs.append(None)
    
    logger.info(f"Generated {sum(1 for x in base64_svgs if x is not None)} SVGs out of {len(df)} rows")
    return base64_svgs