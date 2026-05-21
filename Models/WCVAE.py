import os
import sys
import gc
import random
import numpy as np
import torch
import torch.optim as optim
from tqdm import trange
from config import LOG_DIR, PRETRAIN_DIR, K
from Models.VAE import VAE
from utils.EarlyStopping import EarlyStopping
from utils.Tensorboard import SummaryWriter


def cosine_similarity(a, b):
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    a_norm = np.linalg.norm(a, axis=1, keepdims=True)
    b_norm = np.linalg.norm(b, axis=1, keepdims=True).T
    return a @ b.T / np.maximum(a_norm * b_norm, 1e-12)



def loss_fn(recon_x, x, z, means_encoder, log_var_encoder, means_decoder, log_var_decoder):
    BCE = torch.nn.functional.mse_loss(recon_x, x, reduction='sum') 
    KLD = gaussian_loss(z, means_decoder, log_var_decoder, means_encoder, log_var_encoder)

    return (BCE + KLD) / x.size(0), BCE, KLD  

def log_normal(x, mu, var):
    var = torch.exp(var)
    eps = 1e-8
    if eps > 0.0:
        var = var + eps 
    log_two_pi = torch.log(x.new_tensor(2.0 * np.pi))
    return 0.5 * torch.mean(
        log_two_pi + torch.log(var) + torch.pow(x - mu, 2) / var, dim=-1)

def gaussian_loss(z, z_mu, z_var, z_mu_prior, z_var_prior):
    loss = log_normal(z, z_mu, z_var) - log_normal(z, z_mu_prior, z_var_prior)
    return loss.mean() 

# top-k
def generated_generator(args, device, adj_scipy, features, cancer_name, type):
    features = features.to(torch.device('cpu'))  
    x_list, c_list = [], []  
    kk = K
    for i in trange(adj_scipy.shape[0]):
        c = features[i].unsqueeze(0)
        # ------ KEGG ------
        if type == 'KEGG':
            neighbors_index = list(adj_scipy[i].nonzero()[1])
            if len(neighbors_index) == 0:
                continue

            x = features[neighbors_index]
            sim = np.abs(cosine_similarity(c.numpy(), x.numpy()).flatten())
            weights = abs(sim) /(np.abs(sim).sum() + 1e-8)

            if len(neighbors_index) <= kk:
                weights = weights / np.sum(weights) if np.sum(weights) != 0 else np.ones_like(weights) / len(weights)

            else:
                top_indices = np.argsort(weights)[-kk:]
                x = x[top_indices]
                weights = weights[top_indices]
                weights = weights / np.sum(weights) if np.sum(weights) != 0 else np.ones_like(weights) / len(weights)

        # ------ KNN ------
        elif type == 'KNN':

            neighbors_index = list(adj_scipy[i].nonzero()[1])
            x_real = features[neighbors_index] if neighbors_index else torch.empty((0, features.shape[1]))
            sim_real = np.abs(cosine_similarity(c.numpy(), x_real.numpy()).flatten()) if len(neighbors_index) > 0 else np.array([])

            if len(neighbors_index) < kk:

                candidate_index = np.setdiff1d(np.arange(adj_scipy.shape[0]), np.array([i] + neighbors_index))
                x_cand = features[candidate_index]
                sim_cand = np.abs(cosine_similarity(c.numpy(), x_cand.numpy()).flatten())

                topk_cand_idx = np.argsort(sim_cand)[-(kk - len(neighbors_index)):] if len(sim_cand) > 0 else []
                cand_index_final = candidate_index[topk_cand_idx] if len(sim_cand) > 0 else []
                x_fill = features[cand_index_final] if len(cand_index_final) > 0 else torch.empty((0, features.shape[1]))
                sim_fill = sim_cand[topk_cand_idx] if len(sim_cand) > 0 else np.array([])


                x = torch.cat([x_real, x_fill], dim=0)
                sim = np.concatenate([sim_real, sim_fill]) if sim_real.size > 0 or sim_fill.size > 0 else np.ones(kk) # 如果都没有就均匀
            else:
   
                x = features[neighbors_index]
                sim = np.abs(cosine_similarity(c.numpy(), x.numpy()).flatten())

        
            if x.shape[0] > kk:
                top_indices = np.argsort(sim)[-kk:]
                x = x[top_indices]
                sim = sim[top_indices]
            weights = sim / (np.sum(sim) + 1e-8)
            weights = weights / np.sum(weights) if np.sum(weights) != 0 else np.ones_like(weights) / len(weights)

        else:
            raise ValueError(f"Unknown type: {type}")

        weighted_x = x * torch.tensor(weights).unsqueeze(1).to(torch.float32)
        x_list.append(weighted_x)
        c_list.append(c.repeat(weighted_x.shape[0], 1))
    
    if not x_list or not c_list:
        raise ValueError("no avaliable feature")

    features_x = np.vstack([x.cpu().numpy() for x in x_list])
    features_c = np.vstack([c.cpu().numpy() for c in c_list])


    del x_list, c_list
    gc.collect() 

    hidden_features = 128
    wcvae = VAE(
        encoder_layer_sizes=[features.shape[1], hidden_features],
        latent_size=args.latent_size,
        decoder_layer_sizes=[hidden_features, features.shape[1]],
        conditional=args.conditional,
        conditional_size=features.shape[1]
    ).to(device)
    wcvae_optimizer = optim.Adam(wcvae.parameters(), lr=args.pretrain_lr)
    writer = SummaryWriter(log_dir=str(LOG_DIR))
    early_stopping = EarlyStopping(patience=50, verbose=False)

    for epoch in range(args.total_iterations):
        if len(features_x) == 0 or len(features_c) == 0:
            print("end this epoch, no valid features")
            continue
        
        replace = features_c.shape[0] < args.batch_size
        index = np.random.choice(features_c.shape[0], size=args.batch_size, replace=replace)
        x, c = features_x[index], features_c[index]
        x = torch.tensor(x, dtype=torch.float32).to(device)
        c = torch.tensor(c, dtype=torch.float32).to(device)

        wcvae.train()
        if args.conditional:
            recon_x, z, means_encoder, log_var_encoder, means_decoder, log_var_decoder = wcvae(x, c)
        else:
            recon_x, z, means_encoder, log_var_encoder, means_decoder, log_var_decoder = wcvae(x)

        wcvae_loss, BCE, KLD = loss_fn(recon_x, x, z, means_encoder, log_var_encoder, means_decoder, log_var_decoder)
        print(f"Epoch: {epoch}, BCE: {BCE}, KLD: {KLD}, Loss: {wcvae_loss}")

        writer.add_scalar('wcvae_loss', wcvae_loss, epoch)
        writer.add_scalar('BCE', BCE, epoch)
        writer.add_scalar('KLD', KLD, epoch)

        wcvae_optimizer.zero_grad()
        wcvae_loss.backward()
        wcvae_optimizer.step()

        model_path = PRETRAIN_DIR / cancer_name
        if not os.path.exists(model_path):
            os.makedirs(model_path)
        early_stopping(wcvae_loss.item(), wcvae, path=str(model_path / f"{cancer_name}_{type}_pretrain.pth"))

    
    writer.close()
    del features_x, features_c
    gc.collect()   
    return wcvae
    

