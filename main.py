import argparse
from rdkit import Chem
from rdkit import rdBase
rdBase.DisableLog('rdApp.error')
from gram_utils import *
from smiles_utils import *
from reward import *
from opt import GA, SA

def main():
    parser = argparse.ArgumentParser()
    # Load grammar
    parser.add_argument('--smifile', default='250k_rndm_zinc_drugs_clean.smi')
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--reward', type=callable, default=get_logp, help='reward function')
    parser.add_argument('--gen', type=int, default=10)

    args = parser.parse_args()
    
    np.random.seed(args.seed)

    reward = args.reward
    seed_smiles = []
    with open(args.smifile) as f:
        for line in f:
            smiles = line.rstrip()
            seed_smiles.append(smiles)

    gene_length = 300

    N_mu = 100
    N_lambda = 200

    initial_smiles = np.random.choice(seed_smiles, N_mu+N_lambda)
    initial_smiles = [canonicalize(s) for s in initial_smiles]
    initial_genes = [CFGtoGene(encode(s), max_len=gene_length) for s in initial_smiles]
    initial_scores = [reward(Chem.MolFromSmiles(s)) for s in initial_smiles]

    population = []
    for score, gene, smiles in zip(initial_scores, initial_genes,
                                   initial_smiles):
        population.append((score, smiles, gene))

    population = sorted(population, key=lambda x: x[0], reverse=True)[:N_mu]

    pop = GA.genetic_algorithm(population, N_lambda, N_mu, get_logp, generations=args.gen)
    pop = sorted(pop,key=lambda x: x[0], reverse=True)[:N_mu]

if __name__ == "__main__":
    main()