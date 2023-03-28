import rdkit
from rdkit.Chem import Descriptors
import sys
import os
from rdkit.Chem import RDConfig
from pathlib import Path
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer

def get_sa(mol):
    return sascorer.calculateScore(mol)

def get_logp(mol):
    if mol is None:
        return 0
    return Descriptors.MolLogP(mol)

def get_tpa(mol):
    if mol is None:
        return 0
    return Descriptors.TPSA(mol)

def get_qed(mol):
    if mol is None:
        return 0
    return Descriptors.qed(mol)

