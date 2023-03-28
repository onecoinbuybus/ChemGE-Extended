from rdkit import Chem
from rdkit.AllChem import Descriptors

def get_logp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return Descriptors.MolLogP(mol)