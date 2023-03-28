import numpy as np
from rdkit import Chem
from rdkit import rdBase
import sys
import os
from rdkit.Chem import RDConfig
from pathlib import Path
sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer

def calc_sa(mol):
    return sascorer.calculateScore(mol)