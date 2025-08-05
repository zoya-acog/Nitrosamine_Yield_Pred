from rdkit import Chem

def canonicalize_smiles(smiles):
    """Canonicalize SMILES string."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        else:
            return None
    except:
        return None

def has_nitrogen(smiles):
    """Check if molecule contains nitrogen."""
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    return False
    