import pandas as pd
from rdkit import Chem
import logging
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def denitrosate(smiles):
    """
    Denitrosates a molecule by manually removing the nitroso group.
    Returns the denitrosated SMILES and the index of the parent nitrogen atom.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # Create a mutable copy
    rw_mol = Chem.RWMol(mol)

    # Case 1: N-N=O (e.g., secondary/tertiary nitrosamines)
    # Pattern to find a nitrogen attached to a nitroso group
    pattern1 = Chem.MolFromSmarts("[#7]-[#7](=O)")
    matches1 = rw_mol.GetSubstructMatches(pattern1)
    if matches1:
        # Sort matches to handle multiple nitroso groups if necessary, processing from highest index to avoid re-indexing issues
        for match in sorted(matches1, key=lambda x: x[0], reverse=True):
            parent_n_idx, nitroso_n_idx = match[0], match[1]
            
            # Find the oxygen atom of the nitroso group
            o_idx = -1
            for nbr in rw_mol.GetAtomWithIdx(nitroso_n_idx).GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    o_idx = nbr.GetIdx()
                    break
            
            if o_idx != -1:
                try:
                    # Remove the N=O group and the bond to the parent nitrogen
                    rw_mol.RemoveBond(nitroso_n_idx, o_idx)
                    rw_mol.RemoveBond(parent_n_idx, nitroso_n_idx)
                    # Remove the atoms of the nitroso group, highest index first
                    for idx_to_remove in sorted([o_idx, nitroso_n_idx], reverse=True):
                        rw_mol.RemoveAtom(idx_to_remove)
                    
                    # The molecule is modified in place, but we return the parent N index
                    # from the first successful removal.
                    Chem.SanitizeMol(rw_mol)
                    return Chem.MolToSmiles(rw_mol), parent_n_idx
                except Exception as e:
                    logger.error(f"Error during N-N=O removal for SMILES {smiles}: {e}")
                    return None, None

    # Case 2: N=O (e.g., N-oxides on aromatic systems, treated as nitroso)
    # This pattern is more general, so it's checked after the more specific N-N=O
    pattern2 = Chem.MolFromSmarts("[#7](=O)")
    matches2 = rw_mol.GetSubstructMatches(pattern2)
    if matches2:
        for match in sorted(matches2, key=lambda x: x[0], reverse=True):
            n_idx = match[0]
            
            # Find the attached oxygen
            o_idx = -1
            for nbr in rw_mol.GetAtomWithIdx(n_idx).GetNeighbors():
                if nbr.GetAtomicNum() == 8 and rw_mol.GetBondBetweenAtoms(n_idx, nbr.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                    o_idx = nbr.GetIdx()
                    break

            if o_idx != -1:
                try:
                    # Remove the double bond to oxygen and the oxygen atom
                    rw_mol.RemoveBond(n_idx, o_idx)
                    rw_mol.RemoveAtom(o_idx)
                    
                    Chem.SanitizeMol(rw_mol)
                    return Chem.MolToSmiles(rw_mol), n_idx
                except Exception as e:
                    logger.error(f"Error during N=O removal for SMILES {smiles}: {e}")
                    return None, None

    return None, None

def run_denitrosation(input_csv: str, output_csv: str):
    """
    Runs the denitrosation process on a CSV file and saves the cleaned output.
    """
    logger.info(f"Running denitrosation on {input_csv}")
    
    output_path = Path(output_csv)
    
    # Load CSV
    try:
        df = pd.read_csv(input_csv)
    except FileNotFoundError:
        logger.error(f"Error: Input file not found at {input_csv}")
        return
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise

    recovered_smiles_list = []
    recovered_atom_ids = []
    unrecovered_rows = []

    # === Process all SMILES ===
    for _, row in df.iterrows():
        smiles = row.get('Nitrosated_SMILES')
        if pd.isna(smiles):
            recovered_smiles_list.append(None)
            recovered_atom_ids.append(None)
            unrecovered_rows.append(row)
            continue

        new_smiles, parent_n_idx = denitrosate(smiles)
        if new_smiles is not None:
            recovered_smiles_list.append(new_smiles)
            recovered_atom_ids.append(parent_n_idx)
        else:
            recovered_smiles_list.append(None)
            recovered_atom_ids.append(None)
            unrecovered_rows.append(row)

    df['Recovered_SMILES'] = recovered_smiles_list
    df['Recovered_Atom_ID'] = recovered_atom_ids

    # === Print processing summary ===
    total = len(df)
    fail = len(unrecovered_rows)
    success = total - fail

    logger.info(f"Total SMILES processed: {total}")
    logger.info(f"Successfully denitrosated: {success}")
    logger.info(f"Failed to denitrosate: {fail}")

    if unrecovered_rows:
        unrecovered_csv_path = output_path.parent / "unrecovered_smiles.csv"
        pd.DataFrame(unrecovered_rows).to_csv(unrecovered_csv_path, index=False)
        logger.info(f"Unrecovered SMILES saved to '{unrecovered_csv_path}'.")

    # === Clean invalid rows ===
    
    # Keep track of original size
    original_rows = len(df)
    
    # Create a boolean mask for rows to keep
    rows_to_keep = [True] * len(df)
    invalid_reason_counts = {"bad_smiles_or_atom_id": 0, "not_nitrogen": 0}

    for idx, row in df.iterrows():
        smiles = row['Recovered_SMILES']
        atom_id = row['Recovered_Atom_ID']

        if pd.isna(smiles) or pd.isna(atom_id):
            continue  # Keep rows that were not processed, they will be filtered later

        mol = Chem.MolFromSmiles(smiles)
        atom_id = int(atom_id)

        if mol is None or atom_id >= mol.GetNumAtoms():
            rows_to_keep[idx] = False
            invalid_reason_counts["bad_smiles_or_atom_id"] += 1
            continue

        atom = mol.GetAtomWithIdx(atom_id)
        if atom.GetAtomicNum() != 7:  # Not a Nitrogen atom
            rows_to_keep[idx] = False
            invalid_reason_counts["not_nitrogen"] += 1
    
    # Apply the mask
    df = df[rows_to_keep]

    if invalid_reason_counts["bad_smiles_or_atom_id"] > 0:
        logger.info(f"{invalid_reason_counts['bad_smiles_or_atom_id']} rows removed due to invalid Recovered_SMILES or Atom_ID out of bounds.")
        
    if invalid_reason_counts["not_nitrogen"] > 0:
        logger.info(f"{invalid_reason_counts['not_nitrogen']} rows removed because Recovered_Atom_ID was not a nitrogen atom.")

    # Remove rows where Recovered_SMILES is missing (includes failed denitrosations)
    before_na_drop = len(df)
    df.dropna(subset=['Recovered_SMILES'], inplace=True)
    after_na_drop = len(df)
    
    missing_smiles_count = before_na_drop - after_na_drop
    if missing_smiles_count > 0:
        logger.info(f"{missing_smiles_count} rows removed due to missing Recovered_SMILES.")

    # === Save cleaned output ===
    df.to_csv(output_path, index=False)
    logger.info(f"Final output with {len(df)} rows saved to '{output_path}'")





import pandas as pd
from rdkit import Chem

import pandas as pd
from rdkit import Chem
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def denitrosate(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    rw_mol = Chem.RWMol(mol)

    # Case 1: N(N=O)
    pattern1 = Chem.MolFromSmarts("[#7]-[#7](=O)")
    matches1 = mol.GetSubstructMatches(pattern1)
    if matches1:
        for match in matches1:
            parent_n_idx, nitroso_n_idx = match[0], match[1]
            for nbr in mol.GetAtomWithIdx(nitroso_n_idx).GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    o_idx = nbr.GetIdx()
                    try:
                        rw_mol.RemoveBond(nitroso_n_idx, o_idx)
                        rw_mol.RemoveBond(parent_n_idx, nitroso_n_idx)
                        for idx_to_remove in sorted([o_idx, nitroso_n_idx], reverse=True):
                            rw_mol.RemoveAtom(idx_to_remove)
                        Chem.SanitizeMol(rw_mol)
                        return Chem.MolToSmiles(rw_mol), parent_n_idx
                    except Exception as e:
                        logger.error(f"Error denitrosating (N-N=O): {smiles} — {e}")
                        return None, None

    # Case 2: N(=O)
    pattern2 = Chem.MolFromSmarts("[#7](=O)")
    matches2 = mol.GetSubstructMatches(pattern2)
    if matches2:
        for match in matches2:
            n_idx = match[0]
            for nbr in mol.GetAtomWithIdx(n_idx).GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    o_idx = nbr.GetIdx()
                    try:
                        rw_mol.RemoveBond(n_idx, o_idx)
                        rw_mol.RemoveAtom(o_idx)
                        Chem.SanitizeMol(rw_mol)
                        return Chem.MolToSmiles(rw_mol), n_idx
                    except Exception as e:
                        logger.error(f"Error denitrosating (N=O): {smiles} — {e}")
                        return None, None

    return None, None

def run_denitrosation(input_csv, output_csv):
    logger.info(f"Running denitrosation on {input_csv}")
    
    # Load CSV
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise

    recovered_smiles_list = []
    recovered_atom_ids = []
    unrecovered_rows = []

    # === Process all SMILES ===
    for idx, row in df.iterrows():
        smiles = row['Nitrosated_SMILES']
        if pd.isna(smiles):
            recovered_smiles_list.append(None)
            recovered_atom_ids.append(None)
            unrecovered_rows.append(row)
            continue

        new_smiles, parent_n_idx = denitrosate(smiles)
        if new_smiles:
            recovered_smiles_list.append(new_smiles)
            recovered_atom_ids.append(parent_n_idx)
        else:
            recovered_smiles_list.append(None)
            recovered_atom_ids.append(None)
            unrecovered_rows.append(row)

    df['Recovered_SMILES'] = recovered_smiles_list
    df['Recovered_Atom_ID'] = recovered_atom_ids

    # === Print processing summary ===
    total = len(df)
    success = total - len(unrecovered_rows)
    fail = len(unrecovered_rows)

    logger.info(f"✅ Total SMILES processed: {total}")
    logger.info(f"✅ Successfully denitrosated: {success}")
    logger.info(f"❌ Failed to denitrosate: {fail}")

    if unrecovered_rows:
        unrecovered_csv_path = output_csv.parent / "unrecovered_smiles.csv"
        pd.DataFrame(unrecovered_rows).to_csv(unrecovered_csv_path, index=False)
        logger.info(f"📝 Unrecovered SMILES saved to '{unrecovered_csv_path}'.")

    # === Clean invalid rows ===

    # Step 1: Remove rows where Recovered_Atom_ID is not nitrogen
    valid_n_rows = []
    invalid_n_count = 0

    for idx, row in df.iterrows():
        smiles = row['Recovered_SMILES']
        atom_idx = row['Recovered_Atom_ID']
        if pd.isna(smiles) or pd.isna(atom_idx):
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or int(atom_idx) >= mol.GetNumAtoms():
            invalid_n_count += 1
            continue
        atom = mol.GetAtomWithIdx(int(atom_idx))
        if atom.GetAtomicNum() == 7:  # Nitrogen
            valid_n_rows.append(row)
        else:
            invalid_n_count += 1

    df = pd.DataFrame(valid_n_rows)

    if invalid_n_count > 0:
        logger.info(f"❌ {invalid_n_count} rows removed — Recovered_Atom_ID not nitrogen.")

    # Step 2: Remove rows where Recovered_SMILES is empty
    before = len(df)
    df = df[df['Recovered_SMILES'].notna()]
    after = len(df)
    missing_smiles_count = before - after

    if missing_smiles_count > 0:
        logger.info(f"❌ {missing_smiles_count} rows removed — Recovered_SMILES missing.")

    # === Save cleaned output ===
    df.to_csv(output_csv, index=False)
    logger.info(f"✅ Final output saved to '{output_csv}'")