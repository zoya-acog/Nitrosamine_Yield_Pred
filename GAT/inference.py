# import pandas as pd
# import random
# import numpy as np
# import torch
# import pytorch_lightning as pl
# from torch_geometric.nn import GATConv
# import torch.nn.functional as F
# from rdkit import Chem
# from rdkit.Chem import Draw
# from rdkit.Chem.Draw import rdMolDraw2D
# from PIL import Image
# import base64
# import io, os
# try:
#     from IPython.display import display, HTML
# except ImportError:
#     display, HTML = None, None
# import ast
# from torch_geometric.data import Data

# def preprocess_smiles(smi: str):
#     if not smi or smi == '-' or pd.isna(smi):
#         return None
#     m = Chem.MolFromSmiles(smi.strip())
#     if m is None:
#         return None
#     m = Chem.RemoveHs(m)
#     Chem.RemoveStereochemistry(m)
#     return m

# def identify_nitrosated_nitrogen(parent: Chem.Mol, product: Chem.Mol):
#     mapping = product.GetSubstructMatch(parent)
#     if not mapping:
#         return None
#     p2p = {i: mapping[i] for i in range(len(mapping))}
#     for smarts in ('[N]-[N]=[O]', '[N]=[O]'):
#         patt = Chem.MolFromSmarts(smarts)
#         matches = product.GetSubstructMatches(patt)
#         if matches:
#             prod_idx = matches[0][0]
#             for p_idx, pr_idx in p2p.items():
#                 if pr_idx == prod_idx:
#                     return p_idx
#     for p_idx, pr_idx in p2p.items():
#         atom     = parent.GetAtomWithIdx(p_idx)
#         prod_atom= product.GetAtomWithIdx(pr_idx)
#         if atom.GetAtomicNum()==7 and prod_atom.GetDegree()>atom.GetDegree():
#             return p_idx
#     return None

# def row_to_data(row, yield_threshold=1.0):
#     parent  = preprocess_smiles(row['SMILES'])
#     if parent is None:
#         return None
#     product = preprocess_smiles(row['NA SMILES*'])
#     is_pos  = row['4hr Yield (%)'] >= yield_threshold

#     # --- y labels ---
#     N = parent.GetNumAtoms()
#     y = torch.zeros(N, dtype=torch.long)

#     if is_pos:
#         # 1) explicit indices from output_with_index.csv
#         if row['nitro_indices']:
#             for idx in row['nitro_indices']:
#                 atom = parent.GetAtomWithIdx(idx)
#                 if atom.GetAtomicNum() == 7:  # Check if nitrogen
#                     y[idx] = 1
#         else:
#             # 2) try substructure‐based detect
#             nit = identify_nitrosated_nitrogen(parent, product) if product else None
#             if nit is not None:
#                 y[nit] = 1
#             else:
#                 # 3) fallback: mark all nitrogens
#                 for atom in parent.GetAtoms():
#                     if atom.GetAtomicNum() == 7:
#                         y[atom.GetIdx()] = 1

#     # --- mask: only evaluate nitrogens ---
#     mask = torch.tensor(
#         [atom.GetAtomicNum()==7 for atom in parent.GetAtoms()],
#         dtype=torch.bool
#     )

#     # --- features ---
#     x = torch.tensor([
#         [
#             atom.GetAtomicNum(),
#             atom.GetFormalCharge(),
#             int(atom.GetHybridization()),
#             atom.GetTotalNumHs()
#         ]
#         for atom in parent.GetAtoms()
#     ], dtype=torch.float)

#     # --- edges ---
#     edge_index, edge_attr = [], []
#     for b in parent.GetBonds():
#         i,j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
#         o   = b.GetBondTypeAsDouble()
#         edge_index += [[i,j],[j,i]]
#         edge_attr  += [[o],[o]]
#     edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
#     edge_attr  = torch.tensor(edge_attr, dtype=torch.float)

#     data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr,
#                 y=y, mask=mask)
#     data.smiles = row['SMILES']
#     return data

# def set_seed(seed=42):
#     """
#     Set all seeds to make results reproducible
#     """
#     random.seed(seed)         # for Python's random module
#     np.random.seed(seed)      # for NumPy
#     torch.manual_seed(seed)   # for PyTorch on CPU    
#     torch.cuda.manual_seed(seed) # for PyTorch on GPU
#     torch.cuda.manual_seed_all(seed)  # if using multi-GPU
    
#     # Set deterministic behavior for CUDA
#     if torch.cuda.is_available():
#         torch.backends.cudnn.deterministic = True
#         torch.backends.cudnn.benchmark = False
    
#     os.environ['PYTHONHASHSEED'] = str(seed) # for any other libraries that might use random numbers

#     pl.seed_everything(seed) # for PyTorch Lightning
    
#     print(f"Seed set to {seed}")
    
#     return seed


# class GATNet(torch.nn.Module):
    
#     def __init__(self, in_dim=4, hd=16, heads=4):
#         super().__init__()
#         self.c1 = GATConv(in_dim, hd, heads=heads)
#         self.c2 = GATConv(hd * heads, hd, heads=1)
#         self.lin = torch.nn.Linear(hd, 2)
    
#     def forward(self, data):
#         x, ei = data.x, data.edge_index
#         x = F.elu(self.c1(x, ei))
#         x = F.elu(self.c2(x, ei))
#         return self.lin(x)    
    

# def mol_image_with_labels(mol, highlight_atoms, nitrogen_indices):
#     #mol = Chem.AddHs(mol)  # Add hydrogens explicitly
#     drawer = rdMolDraw2D.MolDraw2DCairo(250, 250)
#     draw_options = drawer.drawOptions()
#     draw_options.includeAtomTags = True
#     for idx in nitrogen_indices:
#         draw_options.atomLabels[idx] = f'N{idx}'
#     drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
#     drawer.FinishDrawing()
#     png_data = drawer.GetDrawingText()
#     img = Image.open(io.BytesIO(png_data))
#     buffered = io.BytesIO()
#     img.save(buffered, format="PNG")
#     img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
#     return f'<img src="data:image/png;base64,{img_str}" width="120"/>'


# def predict_and_visualize(input_csv, model_path, output_csv=None, save_html_path=None):
#     seed = 42
#     set_seed(seed)
#     device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

#     model = GATNet(in_dim=4, hd=16, heads=4).to(device)
#     model.load_state_dict(torch.load(model_path, map_location=device))
#     model.eval()

#     df = pd.read_csv(input_csv)
#     records = []

#     for _, row in df.iterrows():
#         smi = row['SMILES']
#         mol = Chem.MolFromSmiles(smi)
#         if mol is None:
#             records.append({"Name": row.get("Name", ""), "SMILES": smi, "Error": "Invalid SMILES"})
#             continue

#         row = row.copy()
#         has_labels = False
#         if '4hr Yield (%)' in row and 'NA SMILES*' in row and 'nitro_indices' in row:
#             try:
#                 row['4hr Yield (%)'] = float(row.get('4hr Yield (%)', 0) or 0)
#                 if isinstance(row.get('nitro_indices'), str):
#                     row['nitro_indices'] = eval(row['nitro_indices'])
#                 has_labels = True
#             except:
#                 has_labels = False

#         try:
#             data = row_to_data(row) if has_labels else None
#         except:
#             data = None

#         if data is None:
#             parent = preprocess_smiles(smi)
#             if parent is None:
#                 records.append({"Name": row.get("Name", ""), "SMILES": smi, "Error": "Failed to parse SMILES"})
#                 continue

#             mask = torch.tensor([atom.GetAtomicNum() == 7 for atom in parent.GetAtoms()], dtype=torch.bool)
#             x = torch.tensor([
#                 [atom.GetAtomicNum(), atom.GetFormalCharge(), int(atom.GetHybridization()), atom.GetTotalNumHs()]
#                 for atom in parent.GetAtoms()
#             ], dtype=torch.float)
#             edge_index, edge_attr = [], []
#             for b in parent.GetBonds():
#                 i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
#                 o = b.GetBondTypeAsDouble()
#                 edge_index += [[i, j], [j, i]]
#                 edge_attr += [[o], [o]]
#             edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
#             edge_attr = torch.tensor(edge_attr, dtype=torch.float)
#             data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, mask=mask)
#             data.smiles = smi

#         data = data.to(device)
#         with torch.no_grad():
#             out = model(data)
#             probs = torch.softmax(out, dim=1)
#             preds = probs.argmax(dim=1)

#         n_mask = data.mask if hasattr(data, 'mask') else (data.x[:, 0] == 7)
#         n_indices = n_mask.nonzero(as_tuple=True)[0]
#         n_preds = preds[n_indices]
#         n_probs = probs[n_indices, n_preds]
#         n_true = data.y[n_indices] if hasattr(data, 'y') and data.y is not None and len(data.y) == data.x.size(0) else None

#         highlight_atoms = [i for i, p in zip(n_indices.tolist(), n_preds.tolist()) if p == 1]
#         img_html = mol_image_with_labels(mol, highlight_atoms=highlight_atoms, nitrogen_indices=n_indices.tolist())

#         for j, (i, pred, prob) in enumerate(zip(n_indices.tolist(), n_preds.tolist(), n_probs.tolist())):
#             record = {
#                 "Name": row.get("Name", ""),
#                 "SMILES": smi,
#                 "N_Atom_Index": i,
#                 "Pred_Label": int(pred),
#                 "Pred_Prob": round(float(prob), 4),
#                 "Structure": img_html
#             }
#             if n_true is not None:
#                 record["True_Label"] = int(n_true[j].item())
#             records.append(record)

#     result_df = pd.DataFrame(records)

#     def label_class(val):
#         return "Yes" if val == 1 else "No"

#     def format_nitrogen_info(subdf):
#         return "<br>".join([
#             f"N@{row.N_Atom_Index} : " +
#             (f"True={label_class(row.True_Label)}, " if 'True_Label' in row else "") +
#             f"Pred={label_class(row.Pred_Label)}, P={row.Pred_Prob}"
#             for _, row in subdf.iterrows()
#         ])

#     grouped_df = result_df.groupby(["Name", "SMILES", "Structure"]).apply(format_nitrogen_info).reset_index(name="Is N-nitrosatable?")
    
#     final_df = grouped_df.rename(columns={
#         "Name": "name",
#         "SMILES": "smile",
#         "Structure": "structure",
#         "Is N-nitrosatable?": "is n nitrostable?"
#     })

#     if output_csv:
#         if os.path.dirname(output_csv):
#             os.makedirs(os.path.dirname(output_csv), exist_ok=True)
#         final_df.to_csv(output_csv, index=False, quoting=1) # QUOTE_ALL
#         print(f"Saved CSV to {output_csv}")

#     if save_html_path:
#         if os.path.dirname(save_html_path):
#             os.makedirs(os.path.dirname(save_html_path), exist_ok=True)
#         with open(save_html_path, "w") as f:
#             f.write(grouped_df.to_html())
#         print(f"✔️ Saved HTML to {save_html_path}")

#     if display and HTML:
#         try:
#             get_ipython()
#             display(HTML(grouped_df.to_html(escape=False)))
#         except NameError:
#             print("--- Prediction Results ---")
#             print(final_df)
#             print("--------------------------")
#     else:
#         print("--- Prediction Results ---")
#         print(final_df)
#         print("--------------------------")

#     return grouped_df


import pandas as pd
import random
import numpy as np
import torch
import pytorch_lightning as pl
from torch_geometric.nn import GATConv
import torch.nn.functional as F
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image
import base64
import io, os
try:
    from IPython.display import display, HTML
except ImportError:
    display, HTML = None, None
import ast
from torch_geometric.data import Data
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, filename='rule_based_pipeline.log')
logger = logging.getLogger(__name__)

def preprocess_smiles(smi: str):
    if not smi or smi == '-' or pd.isna(smi):
        return None
    m = Chem.MolFromSmiles(smi.strip())
    if m is None:
        return None
    m = Chem.RemoveHs(m)
    Chem.RemoveStereochemistry(m)
    return m

def identify_nitrosated_nitrogen(parent: Chem.Mol, product: Chem.Mol):
    mapping = product.GetSubstructMatch(parent)
    if not mapping:
        return None
    p2p = {i: mapping[i] for i in range(len(mapping))}
    for smarts in ('[N]-[N]=[O]', '[N]=[O]'):
        patt = Chem.MolFromSmarts(smarts)
        matches = product.GetSubstructMatches(patt)
        if matches:
            prod_idx = matches[0][0]
            for p_idx, pr_idx in p2p.items():
                if pr_idx == prod_idx:
                    return p_idx
    for p_idx, pr_idx in p2p.items():
        atom     = parent.GetAtomWithIdx(p_idx)
        prod_atom= product.GetAtomWithIdx(pr_idx)
        if atom.GetAtomicNum()==7 and prod_atom.GetDegree()>atom.GetDegree():
            return p_idx
    return None

def row_to_data(row, yield_threshold=1.0):
    parent  = preprocess_smiles(row['SMILES'])
    if parent is None:
        return None
    product = preprocess_smiles(row['NA SMILES*'])
    is_pos  = row['4hr Yield (%)'] >= yield_threshold

    # --- y labels ---
    N = parent.GetNumAtoms()
    y = torch.zeros(N, dtype=torch.long)

    if is_pos:
        # 1) explicit indices from output_with_index.csv
        if row['nitro_indices']:
            for idx in row['nitro_indices']:
                atom = parent.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:  # Check if nitrogen
                    y[idx] = 1
        else:
            # 2) try substructure‐based detect
            nit = identify_nitrosated_nitrogen(parent, product) if product else None
            if nit is not None:
                y[nit] = 1
            else:
                # 3) fallback: mark all nitrogens
                for atom in parent.GetAtoms():
                    if atom.GetAtomicNum() == 7:
                        y[atom.GetIdx()] = 1

    # --- mask: only evaluate nitrogens ---
    mask = torch.tensor(
        [atom.GetAtomicNum()==7 for atom in parent.GetAtoms()],
        dtype=torch.bool
    )

    # --- features ---
    x = torch.tensor([
        [
            atom.GetAtomicNum(),
            atom.GetFormalCharge(),
            int(atom.GetHybridization()),
            atom.GetTotalNumHs()
        ]
        for atom in parent.GetAtoms()
    ], dtype=torch.float)

    # --- edges ---
    edge_index, edge_attr = [], []
    for b in parent.GetBonds():
        i,j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        o   = b.GetBondTypeAsDouble()
        edge_index += [[i,j],[j,i]]
        edge_attr  += [[o],[o]]
    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    edge_attr  = torch.tensor(edge_attr, dtype=torch.float)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr,
                y=y, mask=mask)
    data.smiles = row['SMILES']
    return data

def set_seed(seed=42):
    """
    Set all seeds to make results reproducible
    """
    random.seed(seed)         # for Python's random module
    np.random.seed(seed)      # for NumPy
    torch.manual_seed(seed)   # for PyTorch on CPU    
    torch.cuda.manual_seed(seed) # for PyTorch on GPU
    torch.cuda.manual_seed_all(seed)  # if using multi-GPU
    
    # Set deterministic behavior for CUDA
    if torch.cuda.is_available():
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    
    os.environ['PYTHONHASHSEED'] = str(seed) # for any other libraries that might use random numbers

    pl.seed_everything(seed) # for PyTorch Lightning
    
    print(f"Seed set to {seed}")
    
    return seed


class GATNet(torch.nn.Module):
    
    def __init__(self, in_dim=4, hd=16, heads=4):
        super().__init__()
        self.c1 = GATConv(in_dim, hd, heads=heads)
        self.c2 = GATConv(hd * heads, hd, heads=1)
        self.lin = torch.nn.Linear(hd, 2)
    
    def forward(self, data):
        x, ei = data.x, data.edge_index
        x = F.elu(self.c1(x, ei))
        x = F.elu(self.c2(x, ei))
        return self.lin(x)    

def mol_image_with_labels(mol, highlight_atoms, nitrogen_indices):
    """Generate SVG image for molecule with highlighted atoms and labeled nitrogens."""
    try:
        # Ensure molecule is valid
        if mol is None:
            logger.warning(f"Invalid molecule provided to mol_image_with_labels")
            return None
        
        # Compute 2D coordinates
        Chem.rdDepictor.Compute2DCoords(mol)
        
        # Create SVG drawer
        drawer = rdMolDraw2D.MolDraw2DSVG(200, 150)  # Match size with rule-based
        draw_options = drawer.drawOptions()
        draw_options.addAtomIndices = True
        draw_options.padding = 0.2
        
        # Add custom labels for nitrogen atoms
        for idx in nitrogen_indices:
            draw_options.atomLabels[idx] = f'N{idx}'
        
        # Draw molecule with highlighted atoms
        drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
        drawer.FinishDrawing()
        
        # Get SVG text
        svg = drawer.GetDrawingText().replace('svg:', '')
        
        # Encode to base64
        base64_svg = base64.b64encode(svg.encode('utf-8')).decode('utf-8')
        logger.info(f"Generated SVG for molecule with {len(highlight_atoms)} highlighted atoms")
        return base64_svg
    
    except Exception as e:
        logger.error(f"Failed to generate SVG for molecule: {e}")
        return None

def predict_and_visualize(input_csv, model_path, output_csv=None, save_html_path=None):
    seed = 42
    set_seed(seed)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = GATNet(in_dim=4, hd=16, heads=4).to(device)
    model.load_state_dict(torch.load(model_path, map_location=device))
    model.eval()

    df = pd.read_csv(input_csv)
    records = []

    for _, row in df.iterrows():
        smi = row['SMILES']
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            records.append({"Name": row.get("Name", ""), "SMILES": smi, "Error": "Invalid SMILES"})
            continue

        row = row.copy()
        has_labels = False
        if '4hr Yield (%)' in row and 'NA SMILES*' in row and 'nitro_indices' in row:
            try:
                row['4hr Yield (%)'] = float(row.get('4hr Yield (%)', 0) or 0)
                if isinstance(row.get('nitro_indices'), str):
                    row['nitro_indices'] = eval(row['nitro_indices'])
                has_labels = True
            except:
                has_labels = False

        try:
            data = row_to_data(row) if has_labels else None
        except:
            data = None

        if data is None:
            parent = preprocess_smiles(smi)
            if parent is None:
                records.append({"Name": row.get("Name", ""), "SMILES": smi, "Error": "Failed to parse SMILES"})
                continue

            mask = torch.tensor([atom.GetAtomicNum() == 7 for atom in parent.GetAtoms()], dtype=torch.bool)
            x = torch.tensor([
                [atom.GetAtomicNum(), atom.GetFormalCharge(), int(atom.GetHybridization()), atom.GetTotalNumHs()]
                for atom in parent.GetAtoms()
            ], dtype=torch.float)
            edge_index, edge_attr = [], []
            for b in parent.GetBonds():
                i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
                o = b.GetBondTypeAsDouble()
                edge_index += [[i, j], [j, i]]
                edge_attr += [[o], [o]]
            edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
            edge_attr = torch.tensor(edge_attr, dtype=torch.float)
            data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr, mask=mask)
            data.smiles = smi

        data = data.to(device)
        with torch.no_grad():
            out = model(data)
            probs = torch.softmax(out, dim=1)
            preds = probs.argmax(dim=1)

        n_mask = data.mask if hasattr(data, 'mask') else (data.x[:, 0] == 7)
        n_indices = n_mask.nonzero(as_tuple=True)[0]
        n_preds = preds[n_indices]
        n_probs = probs[n_indices, n_preds]
        n_true = data.y[n_indices] if hasattr(data, 'y') and data.y is not None and len(data.y) == data.x.size(0) else None

        highlight_atoms = [i for i, p in zip(n_indices.tolist(), n_preds.tolist()) if p == 1]
        img_html = mol_image_with_labels(mol, highlight_atoms=highlight_atoms, nitrogen_indices=n_indices.tolist())

        for j, (i, pred, prob) in enumerate(zip(n_indices.tolist(), n_preds.tolist(), n_probs.tolist())):
            record = {
                "Name": row.get("Name", ""),
                "SMILES": smi,
                "N_Atom_Index": i,
                "Pred_Label": int(pred),
                "Pred_Prob": round(float(prob), 4),
                "Structure": img_html
            }
            if n_true is not None:
                record["True_Label"] = int(n_true[j].item())
            records.append(record)

    result_df = pd.DataFrame(records)

    def label_class(val):
        return "Yes" if val == 1 else "No"

    def format_nitrogen_info(subdf):
        return "<br>".join([
            f"N-{row.N_Atom_Index} : " +
            (f"True={label_class(row.True_Label)}, " if 'True_Label' in row else "") +
            f"Pred={label_class(row.Pred_Label)}, P={row.Pred_Prob * 100:.2f}%"
            for _, row in subdf.iterrows()
        ])

    grouped_df = result_df.groupby(["Name", "SMILES", "Structure"]).apply(format_nitrogen_info).reset_index(name="Is N-nitrosatable?")
    
    final_df = grouped_df.rename(columns={
        "Name": "name",
        "SMILES": "SMILES",
        "Structure": "structure",
        "Is N-nitrosatable?": "is n nitrostable?"
    })

    if output_csv:
        if os.path.dirname(output_csv):
            os.makedirs(os.path.dirname(output_csv), exist_ok=True)
        final_df.to_csv(output_csv, index=False, quoting=1) # QUOTE_ALL
        print(f"Saved CSV to {output_csv}")

    if save_html_path:
        if os.path.dirname(save_html_path):
            os.makedirs(os.path.dirname(save_html_path), exist_ok=True)
        with open(save_html_path, "w") as f:
            f.write(grouped_df.to_html())
        print(f"✔️ Saved HTML to {save_html_path}")

    if display and HTML:
        try:
            get_ipython()
            display(HTML(grouped_df.to_html(escape=False)))
        except NameError:
            print("--- Prediction Results ---")
            print(final_df)
            print("--------------------------")
    else:
        print("--- Prediction Results ---")
        print(final_df)
        print("--------------------------")

    return grouped_df