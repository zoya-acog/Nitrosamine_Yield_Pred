import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def visualize_nitrosamine_products(input_csv, output_dir):
    """Visualize nitrosamine products with highlighted N=O group."""
    logger.info(f"Visualizing nitrosamine products from {input_csv}")
    
    # Load the CSV
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise
    
    # SMARTS for N=O functional group
    pattern = Chem.MolFromSmarts("N(=O)")
    
    # Create visualization directory
    vis_dir = output_dir / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each row
    processed_count = 0
    for idx, row in df.iterrows():
        smiles = row['Nitrosated_SMILES']
        name = row.get('Name', f"Row_{idx}")
        n_type = row.get('Type_of_nitrogen', 'Unknown')
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logger.warning(f"Invalid SMILES at row {idx}: {smiles}")
            continue
        
        Chem.rdDepictor.Compute2DCoords(mol)
        
        matches = mol.GetSubstructMatch(pattern)
        if not matches:
            logger.warning(f"No N=O group found in row {idx}: {smiles}")
            continue
        
        # Draw molecule with N=O highlighted
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol, highlightAtoms=list(matches))
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        
        # Save SVG
        svg_file = vis_dir / f"nitrosamine_{name}_{idx}.svg"
        with open(svg_file, 'w') as f:
            f.write(svg)
        logger.info(f"Saved visualization for {name} (Type: {n_type}, N=O at atoms: {matches}) to {svg_file}")
        
        processed_count += 1
    
    logger.info(f"Processed {processed_count} SMILES for visualization")