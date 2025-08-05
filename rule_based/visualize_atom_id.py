import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def visualize_denitrosated_atoms(input_csv, output_dir):
    """Visualize denitrosated molecules with highlighted atom IDs."""
    logger.info(f"Visualizing denitrosated atoms from {input_csv}")
    
    # Load the CSV
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise
    
    # Create visualization directory
    vis_dir = output_dir / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)
    
    # Loop and draw for a few molecules
    for idx, row in df.iterrows():
        smiles = row.get("Recovered_SMILES")
        atom_id = row.get("Recovered_Atom_ID")
        
        if pd.isna(smiles):
            continue
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES: {smiles}")
            continue
        
        Chem.rdDepictor.Compute2DCoords(mol)
        
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
        drawer.drawOptions().addAtomIndices = True
        
        highlight_atoms = [int(atom_id)] if pd.notna(atom_id) else []
        drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
        drawer.FinishDrawing()
        
        svg = drawer.GetDrawingText().replace('svg:', '')
        
        # Save SVG
        name = row.get('Name', f"Row_{idx}")
        svg_file = vis_dir / f"denitrosated_{name}_{idx}.svg"
        with open(svg_file, 'w') as f:
            f.write(svg)
        logger.info(f"Saved visualization for {name} (SMILES: {smiles}, Atom ID: {atom_id}) to {svg_file}")
        
        # Limit display (optional)
        if idx > 9:  # Show only first 10
            break