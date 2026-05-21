import argparse
import os
import numpy as np
import random
import scipy.sparse as sp
from Models import WCVAE
import torch
import torch_geometric.transforms as T

from config import DATA_DIR, DEFAULT_CANCERS, K
from utils.DataProcessing import CancerDataset

import warnings

# ignore FutureWarning & UserWarning
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)



def print_setting(args):
    print('\n===========================')
    for k, v, in args.__dict__.items():
        print('%s: %s' % (k, v))
    print('===========================\n')


if __name__ == '__main__':
    cancers = [name for name in DEFAULT_CANCERS if (DATA_DIR / name).exists()]

    parser = argparse.ArgumentParser()

    parser.add_argument("--conditional", action='store_true', default=True)
    parser.add_argument("--latent_size", type=int, default=10)
    parser.add_argument("--pretrain_lr", type=float, default=0.001)
    parser.add_argument("--total_iterations", type=int, default=10000)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument('--seed', type=int, default=42, help='Random seed.')

    args = parser.parse_args()
    print_setting(args)

    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(args.seed)
    np.random.seed(args.seed)
    random.seed(args.seed)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    for cancer_name in cancers:
        print('loading dataset--------{}'.format(cancer_name))
        dataset = CancerDataset(root=str(DATA_DIR / cancer_name), transform=T.ToSparseTensor())
        print('dataset:',dataset)
        print(len(dataset))
        data = dataset[0]
 
    
        print("===================Graph===================")
        adj_scipy = sp.csr_matrix(data.adj_t.to_scipy())
        # adj_scipy = build_directed_topk_adj(data.x, K)
        data = data.to(device)
        types = ['KEGG', 'KNN']
        for type in types:
            cvae_model = WCVAE.generated_generator(args, device, adj_scipy, data.x, cancer_name, type)
