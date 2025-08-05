import argparse
import os
import sys
import logging
from pathlib import Path
from GAT.inference import predict_and_visualize

# Add current directory to sys.path
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

from rule_based.process_csv import process_csv
from rule_based.identify_nitrogen import identify_nitrogen_centers
from rule_based.nitrosation_reaction import run_nitrosation
from rule_based.denitrosation import run_denitrosation
from rule_based.rule_based_logic import apply_rule_based_logic
from rule_based.visualize_nitrogen import visualize_nitrogen_centers
from rule_based.visualize_nitrosamine import visualize_nitrosamine_products
from rule_based.visualize_atom_id import visualize_denitrosated_atoms

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("rule_based_pipeline.log"),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)

def setup_output_dir(output_dir):
    """Create output directory if it doesn't exist."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    return output_path

def validate_input_file(input_file, model_type):
    """Validate if input CSV file exists and has required columns based on model type."""
    if not os.path.exists(input_file):
        logger.error(f"Input file {input_file} does not exist.")
        sys.exit(1)
    try:
        import pandas as pd
        df = pd.read_csv(input_file)
        
        if model_type == "rule":
            required_columns = ["Name", "SMILES", "pKa_25"]
        elif model_type == "gat":
            required_columns = ["Name", "SMILES"] # GAT only needs Name and SMILES
        else:
            logger.error(f"Unknown model type: {model_type}")
            sys.exit(1)

        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            logger.error(f"Input CSV missing required columns: {missing}")
            sys.exit(1)
    except Exception as e:
        logger.error(f"Error reading input CSV {input_file}: {e}")
        sys.exit(1)
    return True

def run_pipeline(input_file, output_dir, steps, visualize):
    """Run the specified steps of the rule-based pipeline."""
    output_path = setup_output_dir(output_dir)
    logger.info(f"Starting pipeline with input: {input_file}, output_dir: {output_dir}")

    # Define default output file names
    processed_csv = output_path / "processed_nitro_data.csv"
    nitrogen_csv = output_path / "nitrogens_with_atom_ids.csv"
    nitrosated_csv = output_path / "nitrosated_products.csv"
    denitrosated_csv = output_path / "denitrosated_products.csv"
    final_csv = output_path / "new_with_atom_id.csv"

    # Step 1: Process CSV
    if "all" in steps or "process_csv" in steps:
        logger.info("-" * 50)
        logger.info("Running Step 1: Process CSV")
        try:
            process_csv(input_file, processed_csv)
            logger.info(f"Processed CSV saved to {processed_csv}")
        except Exception as e:
            logger.error(f"Error in process_csv: {e}")
            sys.exit(1)

    # Step 2: Identify Nitrogen Centers
    if "all" in steps or "identify_nitrogen" in steps:
        if not processed_csv.exists():
            logger.error(f"Required file {processed_csv} not found for identify_nitrogen")
            sys.exit(1)
        logger.info("-" * 50)
        logger.info("Running Step 2: Identify Nitrogen Centers")
        try:
            identify_nitrogen_centers(processed_csv, nitrogen_csv)
            logger.info(f"Nitrogen centers CSV saved to {nitrogen_csv}")
        except Exception as e:
            logger.error(f"Error in identify_nitrogen: {e}")
            sys.exit(1)

    # Optional: Visualize Nitrogen Centers
    if visualize and ("all" in steps or "identify_nitrogen" in steps):
        logger.info("-" * 50)
        logger.info("Visualizing Nitrogen Centers")
        try:
            visualize_nitrogen_centers(nitrogen_csv, output_path)
            logger.info(f"Nitrogen center visualizations saved in {output_path}")
        except Exception as e:
            logger.error(f"Error in visualize_nitrogen: {e}")

    # Step 3: Nitrosation Reaction
    if "all" in steps or "nitrosation" in steps:
        if not nitrogen_csv.exists():
            logger.error(f"Required file {nitrogen_csv} not found for nitrosation")
            sys.exit(1)
        logger.info("-" * 50)
        logger.info("Running Step 3: Nitrosation Reaction")
        try:
            run_nitrosation(nitrogen_csv, nitrosated_csv)
            logger.info(f"Nitrosated products CSV saved to {nitrosated_csv}")
        except Exception as e:
            logger.error(f"Error in nitrosation: {e}")
            sys.exit(1)

    # Optional: Visualize Nitrosamine Products
    if visualize and ("all" in steps or "nitrosation" in steps):
        logger.info("-" * 50)
        logger.info("Visualizing Nitrosamine Products")
        try:
            visualize_nitrosamine_products(nitrosated_csv, output_path)
            logger.info(f"Nitrosamine visualizations saved in {output_path}")
        except Exception as e:
            logger.error(f"Error in visualize_nitrosamine: {e}")

    # Step 4: Denitrosation Reaction
    if "all" in steps or "denitrosation" in steps:
        if not nitrosated_csv.exists():
            logger.error(f"Required file {nitrosated_csv} not found for denitrosation")
            sys.exit(1)
        logger.info("-" * 50)
        logger.info("Running Step 4: Denitrosation Reaction")
        try:
            run_denitrosation(nitrosated_csv, denitrosated_csv)
            logger.info(f"Denitrosated products CSV saved to {denitrosated_csv}")
        except Exception as e:
            logger.error(f"Error in denitrosation: {e}")
            sys.exit(1)

    # Optional: Visualize Denitrosated Atoms
    if visualize and ("all" in steps or "denitrosation" in steps):
        logger.info("-" * 50)
        logger.info("Visualizing Denitrosated Atoms")
        try:
            visualize_denitrosated_atoms(denitrosated_csv, output_path)
            logger.info(f"Denitrosated atom visualizations saved in {output_path}")
        except Exception as e:
            logger.error(f"Error in visualize_denitrosated: {e}")

    # Step 5: Apply Rule-Based Logic
    if "all" in steps or "rule_based" in steps:
        if not denitrosated_csv.exists():
            logger.error(f"Required file {denitrosated_csv} not found for rule_based")
            sys.exit(1)
        logger.info("-" * 50)
        logger.info("Running Step 5: Rule-Based Logic")
        try:
            apply_rule_based_logic(denitrosated_csv, final_csv)
            logger.info(f"Final output CSV saved to {final_csv}")
        except Exception as e:
            logger.error(f"Error in rule_based_logic: {e}")
            sys.exit(1)

    logger.info("Pipeline completed successfully!")

def run_gat_predict(input_file, model_path, output_dir, save_html=False):
    """Run GAT prediction on input CSV."""
    output_path = setup_output_dir(output_dir)
    csv_path = output_path / "gat_predictions.csv" if not save_html else None
    logger.info(f"Running GAT prediction with input: {input_file}, model: {model_path}, output: {csv_path}")
    try:
        predict_and_visualize(input_file, model_path, output_csv=csv_path)
        logger.info(f"GAT prediction {'saved to' if not save_html else 'displayed'} successfully")
    except Exception as e:
        logger.error(f"Error in GAT prediction: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="CLI for Rule-Based and GAT Chemical Pipeline")
    subparsers = parser.add_subparsers(dest="command")

    # Rule-based pipeline parser
    rule_parser = subparsers.add_parser("rule", help="Run rule-based pipeline")
    rule_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input CSV file (e.g., multiple_n.csv)"
    )
    rule_parser.add_argument(
        "-o", "--output-dir",
        default="data",
        help="Directory to save output files (default: data/)"
    )
    rule_parser.add_argument(
        "-s", "--steps",
        nargs="+",
        choices=["all", "process_csv", "identify_nitrogen", "nitrosation", "denitrosation", "rule_based"],
        default=["all"],
        help="Pipeline steps to run (default: all)"
    )
    rule_parser.add_argument(
        "--visualize",
        action="store_true",
        help="Generate visualizations for applicable steps"
    )

    # GAT prediction parser
    gat_parser = subparsers.add_parser("gat", help="Run GAT prediction")
    gat_parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input CSV file"
    )
    gat_parser.add_argument(
        "-m", "--model",
        required=True,
        help="Path to trained model file"
    )
    gat_parser.add_argument(
        "-o", "--output-dir",
        default="data",
        help="Directory to save output files (default: data/)"
    )
    gat_parser.add_argument(
        "--save-html",
        action="store_true",
        help="Save predictions as HTML file (deprecated, use for compatibility)"
    )

    args = parser.parse_args()

    if args.command == "rule":
        validate_input_file(args.input, args.command)
        run_pipeline(args.input, args.output_dir, args.steps, args.visualize)
    elif args.command == "gat":
        validate_input_file(args.input, args.command)
        run_gat_predict(args.input, args.model, args.output_dir, args.save_html)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()