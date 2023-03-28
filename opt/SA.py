import numpy as np
from rdkit import Chem
import copy
import sys
import os
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
sys.path.append(parent_dir)
from gram_utils import *
from smiles_utils import canonicalize

def mutation(gene):
    idx = np.random.choice(len(gene))
    gene_mutant = copy.deepcopy(gene)
    gene_mutant[idx] = np.random.randint(0, 256)
    return gene_mutant

def simulated_annealing(ini_smi, reward, n_iters, temp, beta):
    init = np.ones(X.shape[1])
    best = init
    best_logp = logp_calc(best)
    curr, curr_eval = best, best_logp 
    route = []
    for i in range(n_iters):
        molecule_new = mutation(curr,0.8)
        molecule_new_c = molecule_new.copy()
        logp_new = calculator(molecule_new)
        route.append([molecule_new_c,logp_new])
        if logp_new > best_logp:
            best, best_logp = molecule_new, logp_new
            print('n_iter:',i, 'best_logp: ', best_logp)
        diff = logp_new - curr_eval
        t = temp / float(i + 1)
        metropolis = np.exp(diff / t)
        if diff > 0 or rand() < metropolis:
            curr, curr_eval = molecule_new, logp_new 
            curr_c = curr.copy()
    return [best, best_logp],route