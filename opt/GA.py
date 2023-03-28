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

def genetic_algorithm(population, N_lambda, N_mu, reward, generations=100):
    all_smiles = []

    for generation in range(generations):
        scores = [p[0] for p in population]
        mean_score = np.mean(scores)
        min_score = np.min(scores)
        std_score = np.std(scores)
        best_score = np.max(scores)
        idx = np.argmax(scores)
        best_smiles = population[idx][1]
        print("%{},{},{},{},{}".format(generation, best_score, mean_score, min_score, std_score))

        new_population = []
        for _ in range(N_lambda):
            p = population[np.random.randint(len(population))]
            p_gene = p[2]
            c_gene = mutation(p_gene)

            c_smiles = canonicalize(decode(GenetoCFG(c_gene)))
            if c_smiles not in all_smiles:
                c_score = reward(Chem.MolFromSmiles(c_smiles))
                c = (c_score, c_smiles, c_gene)
                new_population.append(c)
                all_smiles.append(c_smiles)

        population.extend(new_population)
        population = sorted(population, key=lambda x: x[0], reverse=True)[:N_mu]

    return population
