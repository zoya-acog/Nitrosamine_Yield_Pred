import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
import logging
from .utils import canonicalize_smiles

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def identify_nitrogen_centers(input_csv, output_csv):
    """Identify secondary and tertiary nitrogen centers and save to CSV."""
    logger.info(f"Identifying nitrogen centers from {input_csv}")
    
    # Step 1: Load processed SMILES
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise
    
    # Step 2: Expand SMILES based on qualifying nitrogen atoms
    expanded_data = []
    original_count = len(df)
    invalid_smiles_set = set()
    valid_with_output_set = set()
    
    for idx, row in df.iterrows():
        smiles = row['canonical_SMILES']
        mol = Chem.MolFromSmiles(smiles)
        
        if not mol:
            invalid_smiles_set.add(smiles)
            continue
        
        qualifying_found = False
        
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() != 7:  # Only nitrogen
                continue
            
            if atom.GetFormalCharge() != 0:
                continue  # skip charged nitrogen
            
            h_count = atom.GetTotalNumHs()
            neighbors = atom.GetNeighbors()
            atom_idx = atom.GetIdx()
            
            # Check if it's part of an aromatic system
            is_aromatic = atom.GetIsAromatic()
            
            # Secondary: bonded to 2 atoms + 1 H
            if len(neighbors) == 2 and h_count == 1:
                n_type = "Secondary"
            # Tertiary: bonded to 3 atoms, usually 0 H
            elif len(neighbors) == 3 and h_count == 0:
                n_type = "Tertiary"
            # Also consider tertiary if part of ring and has 3 heavy neighbors
            elif len(neighbors) == 3 and atom.IsInRing():
                n_type = "Tertiary"
            else:
                continue
            
            qualifying_found = True
            valid_with_output_set.add(smiles)
            
            new_row = row.to_dict()
            new_row['Atom_ID'] = atom_idx
            new_row['Type_of_nitrogen'] = n_type
            expanded_data.append(new_row)
        
        if not qualifying_found:
            invalid_smiles_set.add(smiles)
    
    # Step 3: Convert to final DataFrame
    expanded_df = pd.DataFrame(expanded_data)
    generated_count = len(expanded_df)
    invalid_count = len(invalid_smiles_set)
    
    # Step 4: Print Summary
    logger.info(f"{original_count} SMILES processed")
    logger.info(f"{generated_count} SMILES generated after assigning type of nitrogen")
    logger.info(f"{invalid_count} SMILES had NO secondary or tertiary nitrogen atoms and were skipped")
    
    # Step 5: Bar Plot for Nitrogen Types
    counts = expanded_df['Type_of_nitrogen'].value_counts().reindex(['Secondary', 'Tertiary'], fill_value=0)
    colors = ['#619CFF', '#F8766D']
    
    plt.figure(figsize=(6, 4))
    bars = plt.bar(counts.index, counts.values, color=colors)
    
    for bar, val in zip(bars, counts.values):
        plt.text(bar.get_x() + bar.get_width() / 2, val + 1, str(val),
                 ha='center', fontsize=12)
    
    plt.ylabel("Number of Nitrogen Atoms")
    plt.title("Distribution of Secondary and Tertiary Nitrogen Centers")
    plt.tight_layout()
    plt.savefig(output_csv.parent / "nitrogen_types_distribution.png")
    plt.close()
    logger.info(f"Nitrogen types plot saved to {output_csv.parent / 'nitrogen_types_distribution.png'}")
    
    # Step 6: Save final CSV
    expanded_df.to_csv(output_csv, index=False)
    logger.info(f"Saved nitrogen centers CSV to {output_csv}")