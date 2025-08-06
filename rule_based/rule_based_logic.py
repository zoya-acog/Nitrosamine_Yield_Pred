# import pandas as pd
# from rdkit import Chem
# import logging

# # Set up logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger(__name__)

# def apply_rule_based_logic(input_csv, output_csv):
#     """Apply rule-based logic to score molecules for nitrosamine likelihood."""
#     logger.info(f"Applying rule-based logic on {input_csv}")
    
#     # Load input
#     try:
#         df = pd.read_csv(input_csv)
#     except Exception as e:
#         logger.error(f"Error reading input CSV {input_csv}: {e}")
#         raise
    
#     # Unified Rule Definitions
#     rule_definitions = [
#         # Deactivating Rules (row-aware for Recovered_Atom_ID)
#         ("Tertiary nitrogen (at atom ID)",
#          lambda mol, row: (
#              (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
#              atom.GetDegree() == 3
#          ), -4, True),
        
#         ("Secondary nitrogen (at atom ID)",
#          lambda mol, row: (
#              (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
#              atom.GetDegree() == 2
#          ), -1, True),
        
#         ("Tertiary amine with pKa > 9.5",
#          lambda mol, row: (
#              (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
#              atom.GetDegree() == 3 and
#              atom.GetTotalNumHs(includeNeighbors=True) == 0 and
#              float(row["pKa_25"]) > 9.5
#          ), -2, True),
        
#         ("Secondary amine with pKa > 9.5",
#          lambda mol, row: (
#              (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
#              atom.GetDegree() == 2 and
#              atom.GetTotalNumHs(includeNeighbors=True) == 1 and
#              float(row["pKa_25"]) > 9.5
#          ), -2, True),
        
#         # SMARTS-Based Rules
#         ("Bulky leaving group on tertiary amine", "[N;D3]([C;H0;R0])([C;H0;R0])[C;H2]", -1, False),
#         ("Aromatic tertiary amine with EWG", "[n,c][N;D3]([C,N])([C,N])[C,N]", -1, False),
#         ("Primary amine", "[NX3H2]", -1, False),
#         ("Hydroxy group", "[OX2H1]", -1, False),
#         ("Sulfonamide", "S(=O)(=O)N", -1, False),
#         ("Steric hindrance near N", "[N;D2,D3]([C](C)C)", -1, False),
#         ("Guanidine", "NC(=N)N", -1, False),
#         ("Amide nitrogen in a ring", "[NX3;$([r5,r6]C(=O)N)]", -2, False),
#         ("Amides, ureas, carbamates", "NC(=O)", -3, False),
#         ("Sulfonamides", "NS(=O)(=O)", -2, False),
#         ("Secondary NH flanked by CO or SO2", "[CX3](=O)N[CX3](=O)", -5, False),
#         ("Wrong cleavage on tertiary amine (N–CH2–Ar intact)", "[N;D3]CC1=CC=CC=C1", -8, False),
        
#         # Activating Rules
#         ("Tertiary N with CH2-aryl group", "[N;D3]CC1=CC=CC=C1", +2, False),
#         ("Aromatic dialkyl tertiary amine", "c[N;D3]([C,N])([C,N])", +1, False),
#         ("Tertiary N with nitro (NO2)", "[N;D3][CX4][NX3](=O)=O", +2, False),
#     ]
    
#     # Compile SMARTS and functions
#     compiled_rules = []
#     for desc, pattern, score, row_aware in rule_definitions:
#         if callable(pattern):
#             compiled_rules.append((desc, pattern, score, row_aware))
#         else:
#             try:
#                 smarts = Chem.MolFromSmarts(pattern)
#                 compiled_rules.append((desc, smarts, score, row_aware))
#             except:
#                 logger.warning(f"Failed to parse SMARTS for: {desc}")
    
#     # Likelihood classification
#     def assign_likelihood(score):
#         if score >= -1:
#             return "Very likely"
#         elif score in [-2, -3]:
#             return "Likely"
#         elif score in [-4, -5]:
#             return "Less likely"
#         else:
#             return "Unlikely"
    
#     # Unified rule application
#     def apply_all_rules(mol, row):
#         rule_hits = {}
#         matched_rules = []
#         total_score = 0
        
#         try:
#             atom_id = int(row.get("Recovered_Atom_ID", -1))
#         except:
#             atom_id = -1
        
#         for desc, pattern, score, row_aware in compiled_rules:
#             try:
#                 if callable(pattern):
#                     matched = pattern(mol, row) if row_aware else pattern(mol)
#                 else:
#                     matches = mol.GetSubstructMatches(pattern)
#                     matched = any(atom_id in match for match in matches)
                
#                 if matched:
#                     rule_hits[desc] = score
#                     total_score += score
#                     matched_rules.append(desc)
#                 else:
#                     rule_hits[desc] = 0
#             except:
#                 rule_hits[desc] = 0
        
#         return total_score, matched_rules, rule_hits
    
#     # Main processing
#     scored_rows = []
    
#     for _, row in df.iterrows():
#         name = row["Name"]
#         smiles = row["Recovered_SMILES"]
#         atom_id = row.get("Recovered_Atom_ID")
#         mol = Chem.MolFromSmiles(smiles)
        
#         if mol is None:
#             continue
        
#         score, matched_rules, rule_hits = apply_all_rules(mol, row)
#         likelihood = assign_likelihood(score)
        
#         output = row.copy()
#         output["Matched_Rules"] = ";".join(matched_rules)
#         output["Score"] = score
#         output["Likelihood"] = likelihood
        
#         for desc, _, _, _ in compiled_rules:
#             output[desc] = rule_hits.get(desc, 0)
        
#         scored_rows.append(output)
    
#     # Save output
#     output_df = pd.DataFrame(scored_rows)
#     output_df.to_csv(output_csv, index=False)
#     logger.info(f"Saved scored output to {output_csv}")
    
#     # Summary
#     summary = output_df["Likelihood"].value_counts().reindex(["Very likely", "Likely", "Less likely", "Unlikely"], fill_value=0)
    
#     logger.info("\n✅ Scoring complete. Breakdown by likelihood category:")
#     for label, count in summary.items():
#         logger.info(f"{label:>12}: {count} molecules")


import pandas as pd
from rdkit import Chem
import logging

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def apply_rule_based_logic(input_csv, output_csv):
    """Apply rule-based logic to score molecules for nitrosamine likelihood."""
    logger.info(f"Applying rule-based logic on {input_csv}")
    
    # Load input
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise
    
    # Unified Rule Definitions
    rule_definitions = [
        # Deactivating Rules (row-aware for Recovered_Atom_ID)
        ("Tertiary nitrogen (at atom ID)",
         lambda mol, row: (
             (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
             atom.GetDegree() == 3
         ), -4, True),
        
        ("Secondary nitrogen (at atom ID)",
         lambda mol, row: (
             (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
             atom.GetDegree() == 2
         ), -1, True),
        
        ("Tertiary amine with pKa > 9.5",
         lambda mol, row: (
             (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
             atom.GetDegree() == 3 and
             atom.GetTotalNumHs(includeNeighbors=True) == 0 and
             float(row["pKa_25"]) > 9.5
         ), -2, True),
        
        ("Secondary amine with pKa > 9.5",
         lambda mol, row: (
             (atom := mol.GetAtomWithIdx(int(row["Recovered_Atom_ID"]))).GetAtomicNum() == 7 and
             atom.GetDegree() == 2 and
             atom.GetTotalNumHs(includeNeighbors=True) == 1 and
             float(row["pKa_25"]) > 9.5
         ), -2, True),
        
        # SMARTS-Based Rules
        ("Bulky leaving group on tertiary amine", "[N;D3]([C;H0;R0])([C;H0;R0])[C;H2]", -1, False),
        ("Aromatic tertiary amine with EWG", "[n,c][N;D3]([C,N])([C,N])[C,N]", -1, False),
        ("Primary amine", "[NX3H2]", -1, False),
        ("Hydroxy group", "[OX2H1]", -1, False),
        ("Sulfonamide", "S(=O)(=O)N", -1, False),
        ("Steric hindrance near N", "[N;D2,D3]([C](C)C)", -1, False),
        ("Guanidine", "NC(=N)N", -1, False),
        ("Amide nitrogen in a ring", "[NX3;$([r5,r6]C(=O)N)]", -2, False),
        ("Amides, ureas, carbamates", "NC(=O)", -3, False),
        ("Sulfonamides", "NS(=O)(=O)", -2, False),
        ("Secondary NH flanked by CO or SO2", "[CX3](=O)N[CX3](=O)", -5, False),
        ("Wrong cleavage on tertiary amine (N–CH2–Ar intact)", "[N;D3]CC1=CC=CC=C1", -8, False),
        
        # Activating Rules
        ("Tertiary N with CH2-aryl group", "[N;D3]CC1=CC=CC=C1", +2, False),
        ("Aromatic dialkyl tertiary amine", "c[N;D3]([C,N])([C,N])", +1, False),
        ("Tertiary N with nitro (NO2)", "[N;D3][CX4][NX3](=O)=O", +2, False),
    ]
    
    # Compile SMARTS and functions
    compiled_rules = []
    for desc, pattern, score, row_aware in rule_definitions:
        if callable(pattern):
            compiled_rules.append((desc, pattern, score, row_aware))
        else:
            try:
                smarts = Chem.MolFromSmarts(pattern)
                compiled_rules.append((desc, smarts, score, row_aware))
            except:
                logger.warning(f"Failed to parse SMARTS for: {desc}")
    
    # Likelihood classification
    def assign_likelihood(score):
        if score >= -1:
            return "Very likely"
        elif score in [-2, -3]:
            return "Likely"
        elif score in [-4, -5]:
            return "Less likely"
        else:
            return "Unlikely"
    
    # Unified rule application
    def apply_all_rules(mol, row):
        rule_hits = {}
        matched_rules = []
        total_score = 0
        
        try:
            atom_id = int(row.get("Recovered_Atom_ID", -1))
        except:
            atom_id = -1
        
        for desc, pattern, score, row_aware in compiled_rules:
            try:
                if callable(pattern):
                    matched = pattern(mol, row) if row_aware else pattern(mol)
                else:
                    matches = mol.GetSubstructMatches(pattern)
                    matched = any(atom_id in match for match in matches)
                
                if matched:
                    rule_hits[desc] = score
                    total_score += score
                    matched_rules.append(desc)
                else:
                    rule_hits[desc] = 0
            except:
                rule_hits[desc] = 0
        
        return total_score, matched_rules, rule_hits
    
    # Main processing
    scored_rows = []
    
    for _, row in df.iterrows():
        name = row["Name"]
        smiles = row["Recovered_SMILES"]
        atom_id = row.get("Recovered_Atom_ID")
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            continue
        
        score, matched_rules, rule_hits = apply_all_rules(mol, row)
        likelihood = assign_likelihood(score)
        
        output = row.copy()
        output["Matched_Rules"] = ";".join(matched_rules)
        output["Score"] = score
        output["Likelihood"] = likelihood
        
        for desc, _, _, _ in compiled_rules:
            output[desc] = rule_hits.get(desc, 0)
        
        scored_rows.append(output)
    
    # Save output
    output_df = pd.DataFrame(scored_rows)
    output_df.to_csv(output_csv, index=False)
    logger.info(f"Saved scored output to {output_csv}")
    
    # Summary
    summary = output_df["Likelihood"].value_counts().reindex(["Very likely", "Likely", "Less likely", "Unlikely"], fill_value=0)
    
    logger.info("\n✅ Scoring complete. Breakdown by likelihood category:")
    for label, count in summary.items():
        logger.info(f"{label:>12}: {count} molecules")
    
    return output_df