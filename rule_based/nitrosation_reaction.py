import pandas as pd
from rdkit import Chem
import matplotlib.pyplot as plt
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def nitrosate_secondary_n(mol, atom_id):
    """Nitrosate secondary nitrogen."""
    mol = Chem.RWMol(mol)
    new_n = Chem.Atom(7)  # N
    new_o = Chem.Atom(8)  # O
    new_n_idx = mol.AddAtom(new_n)
    new_o_idx = mol.AddAtom(new_o)
    mol.AddBond(atom_id, new_n_idx, Chem.rdchem.BondType.SINGLE)
    mol.AddBond(new_n_idx, new_o_idx, Chem.rdchem.BondType.DOUBLE)
    return Chem.MolToSmiles(mol)

def nitrosate_tertiary_n(mol, atom_id):
    """Nitrosate tertiary nitrogen, producing multiple products."""
    mol = Chem.RWMol(mol)
    products = []
    bonded_ids = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(atom_id).GetNeighbors()]
    for b_id in bonded_ids:
        temp_mol = Chem.RWMol(mol)
        try:
            temp_mol.RemoveBond(atom_id, b_id)
            new_n = Chem.Atom(7)
            new_o = Chem.Atom(8)
            new_n_idx = temp_mol.AddAtom(new_n)
            new_o_idx = temp_mol.AddAtom(new_o)
            temp_mol.AddBond(atom_id, new_n_idx, Chem.rdchem.BondType.SINGLE)
            temp_mol.AddBond(new_n_idx, new_o_idx, Chem.rdchem.BondType.DOUBLE)
            smi = Chem.MolToSmiles(temp_mol)
            products.append(smi)
        except:
            continue
    return products

def keep_nitrosated_fragment(smiles):
    """Keep only the fragment containing N=O."""
    fragments = smiles.split('.')
    for frag in fragments:
        mol = Chem.MolFromSmiles(frag)
        if mol:
            patt = Chem.MolFromSmarts("N(=O)")
            if mol.HasSubstructMatch(patt):
                return Chem.MolToSmiles(mol)
    return None

def run_nitrosation(input_csv, output_csv):
    """Run nitrosation reaction on molecules and save results."""
    logger.info(f"Running nitrosation reaction on {input_csv}")
    
    # Load the CSV
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise
    
    product_rows = []
    n_type_counter = {"Secondary": 0, "Tertiary": 0}
    processed = 0
    generated = 0
    failed = 0
    
    # Main loop
    for idx, row in df.iterrows():
        smiles = row['canonical_SMILES']
        atom_id = int(row['Atom_ID'])
        n_type = row['Type_of_nitrogen']
        
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            failed += 1
            continue
        
        processed += 1
        success = False
        
        try:
            if n_type == "Secondary":
                smi = nitrosate_secondary_n(mol, atom_id)
                clean_smi = keep_nitrosated_fragment(smi)
                if clean_smi:
                    new_row = row.copy()
                    new_row['Nitrosated_SMILES'] = clean_smi
                    product_rows.append(new_row)
                    generated += 1
                    n_type_counter["Secondary"] += 1
                    success = True
            
            elif n_type == "Tertiary":
                prod_list = nitrosate_tertiary_n(mol, atom_id)
                for p in prod_list:
                    clean_smi = keep_nitrosated_fragment(p)
                    if clean_smi:
                        new_row = row.copy()
                        new_row['Nitrosated_SMILES'] = clean_smi
                        product_rows.append(new_row)
                        generated += 1
                        n_type_counter["Tertiary"] += 1
                        success = True
            
            if not success:
                failed += 1
        
        except Exception as e:
            failed += 1
            logger.warning(f"Failed to nitrosate row {idx}: {e}")
            continue
    
    # Save Results
    product_df = pd.DataFrame(product_rows)
    product_df.to_csv(output_csv, index=False)
    logger.info(f"Saved nitrosated products to {output_csv}")
    
    # Print Report
    logger.info("==== Nitrosation Summary ====")
    logger.info(f"Total rows in CSV        : {len(df)}")
    logger.info(f"Total SMILES processed   : {processed}")
    logger.info(f"Total nitrosated products: {generated}")
    logger.info(f"Failed to nitrosate      : {failed}")
    logger.info(f"Secondary Nitrogens      : {n_type_counter['Secondary']}")
    logger.info(f"Tertiary Nitrogens       : {n_type_counter['Tertiary']}")
    
    # Plot Bar Chart
    plt.figure(figsize=(6, 4))
    bars = plt.bar(n_type_counter.keys(), n_type_counter.values(), color=['#619CFF', '#F8766D'])
    
    for bar, val in zip(bars, n_type_counter.values()):
        plt.text(bar.get_x() + bar.get_width()/2, val + 1, str(val),
                 ha='center', fontsize=12)
    
    plt.ylabel("Number of Nitrosated Products")
    plt.title("Nitrosated Product Count per Nitrogen Type")
    plt.tight_layout()
    plt.savefig(output_csv.parent / "nitrosated_distribution.png")
    plt.close()
    logger.info(f"Nitrosated product plot saved to {output_csv.parent / 'nitrosated_distribution.png'}")

if __name__ == "__main__":
    # For standalone testing (optional)
    import sys
    if len(sys.argv) != 3:
        print("Usage: python nitrosation_reaction.py <input_csv> <output_csv>")
        sys.exit(1)
    run_nitrosation(sys.argv[1], sys.argv[2])