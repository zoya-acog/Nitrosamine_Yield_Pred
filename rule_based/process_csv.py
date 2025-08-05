import pandas as pd
import matplotlib.pyplot as plt
import logging
from .utils import canonicalize_smiles, has_nitrogen

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def process_csv(input_csv, output_csv):
    """Process input CSV, canonicalize SMILES, check for nitrogen, and plot distribution."""
    logger.info(f"Processing CSV: {input_csv}")
    
    # Step 1: Read CSV and select columns
    columns_to_keep = ['Name', 'SMILES', 'pKa_25']
    columns_to_drop = []  # e.g. ['Comment', 'InChI']
    
    try:
        df = pd.read_csv(input_csv)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_csv}: {e}")
        raise

    # Drop optional columns
    if columns_to_drop:
        df.drop(columns=columns_to_drop, errors='ignore', inplace=True)
    
    df = df[columns_to_keep].copy()
    
    # Step 2: Canonicalize SMILES
    df['canonical_SMILES'] = df['SMILES'].apply(canonicalize_smiles)
    
    # Drop rows with invalid SMILES
    df = df.dropna(subset=['canonical_SMILES']).reset_index(drop=True)
    
    # Step 3: Nitrogen Check
    df['has_nitrogen'] = df['canonical_SMILES'].apply(has_nitrogen)
    
    # Step 4: Save Processed CSV
    df.to_csv(output_csv, index=False)
    logger.info(f"Saved processed CSV to {output_csv}")
    
    # Step 5: Nitrogen Distribution Plot
    counts = df['has_nitrogen'].value_counts()
    counts_dict = {False: 0, True: 0}
    counts_dict.update(counts.to_dict())
    
    labels = ['No Nitrogen', 'Contains Nitrogen']
    values = [counts_dict[False], counts_dict[True]]
    colors = ['#f8766d', '#00ba38']
    
    plt.figure(figsize=(6, 4))
    bars = plt.bar(labels, values, color=colors)
    
    for bar, value in zip(bars, values):
        plt.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                 str(value), ha='center', fontsize=12)
    
    plt.ylabel("Number of Molecules")
    plt.title("Nitrogen Distribution in Canonical SMILES")
    plt.tight_layout()
    plt.savefig(output_csv.parent / "nitrogen_distribution.png")
    plt.close()
    logger.info(f"Nitrogen distribution plot saved to {output_csv.parent / 'nitrogen_distribution.png'}")