import pandas as pd
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def visualize_nitrogen_centers(input_csv, output_dir):
    """Visualize nitrogen centers with highlighted atoms."""
    logger.info(f"Visualizing nitrogen centers from {input_csv}")
    
    # Load the file
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise
    
    # Create visualization directory
    vis_dir = output_dir / "visualizations"
    vis_dir.mkdir(parents=True, exist_ok=True)
    
    # Loop and draw each molecule with highlighted atom
    for idx, row in df.iterrows():
        smiles = row['canonical_SMILES']
        atom_id = int(row['Atom_ID'])
        name = row['Name']
        n_type = row['Type_of_nitrogen']
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            logger.warning(f"Invalid SMILES at row {idx}: {smiles}")
            continue
        
        Chem.rdDepictor.Compute2DCoords(mol)
        
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
        drawer.DrawMolecule(mol, highlightAtoms=[atom_id])
        drawer.FinishDrawing()
        
        svg = drawer.GetDrawingText()
        
        # Save SVG
        svg_file = vis_dir / f"nitrogen_center_{name}_{idx}.svg"
        with open(svg_file, 'w') as f:
            f.write(svg)
        logger.info(f"Saved visualization for {name} (Atom ID: {atom_id}, Type: {n_type}) to {svg_file}")