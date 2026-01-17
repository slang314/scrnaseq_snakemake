#!/usr/bin/env python
"""
Variational Autoencoder for scRNA-seq embedding.
Uses AnnLoader for efficient batching and Pyro for variational inference.

UNSUPERVISED VERSION - no cell type labels used during training.
This avoids circular reasoning where labels derived from the data
are used to train the model that produces embeddings.
"""

import torch
import torch.nn as nn
import pyro
import pyro.distributions as dist
from pyro.infer import SVI, TraceMeanField_ELBO
from pyro.optim import Adam
import numpy as np
import random
import scanpy as sc
from anndata.experimental.pytorch import AnnLoader
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# =============================================================================
# REPRODUCIBILITY
# =============================================================================
RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)
torch.manual_seed(RANDOM_SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(RANDOM_SEED)

# Settings - use v3 output
OUTPUT_DIR = "results/downstream_v3"
os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"Using device: {device}")
print(f"Random seed: {RANDOM_SEED}")

# =============================================================================
# Model Architecture
# =============================================================================

class MLP(nn.Module):
    """Multi-layer perceptron with LayerNorm and Dropout."""
    def __init__(self, input_dim, hidden_dims, out_dim):
        super().__init__()
        modules = []
        dims = [input_dim] + hidden_dims
        for in_size, out_size in zip(dims[:-1], dims[1:]):
            modules.append(nn.Linear(in_size, out_size))
            modules.append(nn.LayerNorm(out_size))
            modules.append(nn.ReLU())
            modules.append(nn.Dropout(p=0.1))
        modules.append(nn.Linear(hidden_dims[-1], out_dim))
        self.fc = nn.Sequential(*modules)

    def forward(self, x):
        return self.fc(x)


class VAE(nn.Module):
    """
    Variational Autoencoder for scRNA-seq data.
    Uses negative binomial likelihood for count data.

    UNSUPERVISED - no cell type classification component.
    """
    def __init__(self, input_dim, hidden_dims, latent_dim):
        super().__init__()
        self.latent_dim = latent_dim
        self.input_dim = input_dim

        # Encoder: x -> (z_mu, z_logvar)
        self.encoder = MLP(input_dim, hidden_dims, 2 * latent_dim)

        # Decoder: z -> reconstructed rates
        self.decoder = MLP(latent_dim, hidden_dims[::-1], input_dim)

        # Dispersion parameter (gene-specific)
        self.theta = nn.Parameter(torch.randn(input_dim))

    def encode(self, x):
        """Encode input to latent distribution parameters."""
        h = self.encoder(x)
        z_mu = h[:, :self.latent_dim]
        z_logvar = h[:, self.latent_dim:]
        return z_mu, z_logvar

    def reparameterize(self, mu, logvar):
        """Reparameterization trick for sampling."""
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z):
        """Decode latent to reconstruction."""
        return self.decoder(z)

    def model(self, x, size_factors):
        """Generative model p(x|z)p(z) - UNSUPERVISED."""
        pyro.module("vae", self)
        batch_size = x.shape[0]

        with pyro.plate("data", batch_size):
            # Prior p(z) = N(0, I)
            z_loc = x.new_zeros((batch_size, self.latent_dim))
            z_scale = x.new_ones((batch_size, self.latent_dim))
            z = pyro.sample("latent", dist.Normal(z_loc, z_scale).to_event(1))

            # Decode to get rate parameters
            dec_out = self.decode(z)
            dec_mu = torch.softmax(dec_out, dim=-1) * size_factors[:, None]

            # Gene-specific dispersion
            theta = torch.exp(self.theta).unsqueeze(0).expand(batch_size, -1)

            # Negative binomial likelihood
            # Using rate parameterization: mean = mu, var = mu + mu^2/theta
            logits = (dec_mu + 1e-6).log() - (theta + 1e-6).log()
            pyro.sample("obs", dist.NegativeBinomial(total_count=theta, logits=logits).to_event(1),
                       obs=x.int())

    def guide(self, x, size_factors):
        """Variational posterior q(z|x) - UNSUPERVISED."""
        batch_size = x.shape[0]

        with pyro.plate("data", batch_size):
            z_mu, z_logvar = self.encode(x)
            z_scale = torch.sqrt(torch.exp(z_logvar) + 1e-4)
            pyro.sample("latent", dist.Normal(z_mu, z_scale).to_event(1))

    def get_latent(self, x):
        """Get latent representation (mean of posterior)."""
        z_mu, _ = self.encode(x)
        return z_mu


# =============================================================================
# Training
# =============================================================================

def train_epoch(svi, dataloader):
    """Train for one epoch - UNSUPERVISED."""
    epoch_loss = 0.0
    for batch in dataloader:
        x = batch.X.float()  # Ensure float32
        size_factors = batch.obs['size_factors'].float()

        loss = svi.step(x, size_factors)
        epoch_loss += loss

    return epoch_loss / len(dataloader.dataset)


def evaluate(vae, dataloader):
    """Evaluate model and get latent embeddings."""
    vae.eval()
    latents = []

    with torch.no_grad():
        for batch in dataloader:
            x = batch.X.float()  # Ensure float32
            z = vae.get_latent(x)
            latents.append(z.cpu().numpy())

    return np.vstack(latents)


# =============================================================================
# Main
# =============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Loading data (UNSUPERVISED VAE)")
    print("=" * 60)

    # Load annotated data from v3 pipeline
    adata = sc.read_h5ad(f"{OUTPUT_DIR}/pbmc_10k_annotated_v3.h5ad")
    print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Use raw counts from layers
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()
        print("Using raw counts from layers['counts']")

    # Compute size factors
    if 'size_factors' not in adata.obs:
        adata.obs['size_factors'] = np.array(adata.X.sum(axis=1)).flatten().astype(np.float32)

    # NOTE: No cell type labels used - this is UNSUPERVISED
    print("Running UNSUPERVISED VAE (no cell type labels used)")

    # Subset to highly variable genes for efficiency
    if 'highly_variable' in adata.var:
        adata_hvg = adata[:, adata.var['highly_variable']].copy()
        print(f"Using {adata_hvg.n_vars} HVGs")
    else:
        adata_hvg = adata

    # Setup converters for AnnLoader
    def to_float32(x):
        arr = x.to_numpy() if hasattr(x, 'to_numpy') else np.array(x)
        return arr.astype(np.float32)

    encoders = {
        'obs': {
            'size_factors': to_float32
        }
    }

    # Create dataloader
    use_cuda = torch.cuda.is_available()
    dataloader = AnnLoader(
        adata_hvg,
        batch_size=256,
        shuffle=True,
        convert=encoders,
        use_cuda=use_cuda
    )

    print(f"\nDataloader created: {len(dataloader)} batches")

    # =============================================================================
    # Initialize model
    # =============================================================================
    print("\n" + "=" * 60)
    print("Initializing VAE (UNSUPERVISED)")
    print("=" * 60)

    input_dim = adata_hvg.n_vars
    hidden_dims = [256, 128]
    latent_dim = 20

    vae = VAE(
        input_dim=input_dim,
        hidden_dims=hidden_dims,
        latent_dim=latent_dim
    )

    if use_cuda:
        vae = vae.cuda()

    print(f"Model: {input_dim} -> {hidden_dims} -> {latent_dim}D latent")
    print(f"Total parameters: {sum(p.numel() for p in vae.parameters()):,}")

    # Setup optimizer and SVI
    pyro.clear_param_store()
    optimizer = Adam({"lr": 1e-3})
    svi = SVI(vae.model, vae.guide, optimizer, loss=TraceMeanField_ELBO())

    # =============================================================================
    # Training (UNSUPERVISED)
    # =============================================================================
    print("\n" + "=" * 60)
    print("Training VAE (UNSUPERVISED)")
    print("=" * 60)

    NUM_EPOCHS = 100
    losses = []

    for epoch in range(NUM_EPOCHS):
        vae.train()
        loss = train_epoch(svi, dataloader)
        losses.append(loss)

        if epoch % 10 == 0 or epoch == NUM_EPOCHS - 1:
            print(f"Epoch {epoch:3d}/{NUM_EPOCHS}: loss = {loss:.4f}")

    # Plot training curve
    plt.figure(figsize=(8, 4))
    plt.plot(losses)
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('VAE Training Loss')
    plt.savefig(f"{OUTPUT_DIR}/figures/vae_training_loss.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: vae_training_loss.png")

    # =============================================================================
    # Get embeddings
    # =============================================================================
    print("\n" + "=" * 60)
    print("Generating embeddings")
    print("=" * 60)

    # Create non-shuffled dataloader for embedding extraction
    dataloader_eval = AnnLoader(
        adata_hvg,
        batch_size=256,
        shuffle=False,
        convert=encoders,
        use_cuda=use_cuda
    )

    latent_embeddings = evaluate(vae, dataloader_eval)
    print(f"Latent embeddings shape: {latent_embeddings.shape}")

    # Add to adata
    adata.obsm['X_vae'] = latent_embeddings

    # =============================================================================
    # UMAP on VAE latent space
    # =============================================================================
    print("\n" + "=" * 60)
    print("Computing UMAP from VAE embeddings")
    print("=" * 60)

    # Store the original PCA-UMAP before overwriting
    if 'X_umap' in adata.obsm:
        adata.obsm['X_umap_pca'] = adata.obsm['X_umap'].copy()
        print("Stored original PCA-UMAP in X_umap_pca")

    # Compute neighbors and UMAP using VAE latent space
    sc.pp.neighbors(adata, use_rep='X_vae', n_neighbors=15, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)

    # Store as separate UMAP
    adata.obsm['X_umap_vae'] = adata.obsm['X_umap'].copy()

    # Re-cluster on VAE space
    sc.tl.leiden(adata, resolution=0.4, key_added='leiden_vae', random_state=RANDOM_SEED)
    n_clusters = adata.obs['leiden_vae'].nunique()
    print(f"Leiden clustering on VAE space: {n_clusters} clusters")

    # =============================================================================
    # Visualizations
    # =============================================================================
    print("\n" + "=" * 60)
    print("Creating visualizations")
    print("=" * 60)

    # Compare PCA-UMAP vs VAE-UMAP
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Left column: PCA-UMAP (stored in X_umap_pca)
    # Right column: VAE-UMAP (stored in X_umap_vae, also current X_umap)
    sc.pl.embedding(adata, basis='umap_pca', color='leiden', ax=axes[0, 0],
                   show=False, title='PCA-UMAP: Leiden clusters (original)')
    sc.pl.embedding(adata, basis='umap_vae', color='leiden_vae', ax=axes[0, 1],
                   show=False, title='VAE-UMAP: Leiden clusters')

    if 'celltypist_label' in adata.obs:
        sc.pl.embedding(adata, basis='umap_pca', color='celltypist_label', ax=axes[1, 0],
                       show=False, title='PCA-UMAP: Celltypist (post-hoc)', legend_fontsize=6)
        sc.pl.embedding(adata, basis='umap_vae', color='celltypist_label', ax=axes[1, 1],
                       show=False, title='VAE-UMAP: Celltypist (post-hoc)', legend_fontsize=6)

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/figures/umap_pca_vs_vae.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: umap_pca_vs_vae.png")

    # VAE UMAP with key markers - need to lookup by gene symbol
    # Create symbol to ID mapping
    if 'gene_symbol' in adata.var.columns:
        symbol_to_id = adata.var.reset_index().set_index('gene_symbol')['index'].to_dict()
        key_markers_symbols = ['CD3D', 'CD14', 'MS4A1', 'GNLY', 'FCGR3A', 'FCER1A']
        key_markers = [symbol_to_id.get(m) for m in key_markers_symbols if m in symbol_to_id]
        key_markers = [m for m in key_markers if m is not None and m in adata.raw.var_names]
    else:
        key_markers = ['CD3D', 'CD14', 'MS4A1', 'GNLY', 'FCGR3A', 'FCER1A']
        key_markers = [m for m in key_markers if m in adata.raw.var_names]

    if key_markers:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        for idx, marker in enumerate(key_markers[:6]):
            ax = axes.flatten()[idx]
            # Get symbol for title
            symbol = adata.var.loc[marker, 'gene_symbol'] if 'gene_symbol' in adata.var.columns else marker
            sc.pl.embedding(adata, basis='umap', color=marker, ax=ax, show=False, title=f'VAE: {symbol}')
        plt.tight_layout()
        plt.savefig(f"{OUTPUT_DIR}/figures/umap_vae_markers.png", dpi=150, bbox_inches='tight')
        plt.close()
        print("Saved: umap_vae_markers.png")

    # =============================================================================
    # Save results
    # =============================================================================
    print("\n" + "=" * 60)
    print("Saving results")
    print("=" * 60)

    # Save updated adata with VAE embeddings
    adata.write_h5ad(f"{OUTPUT_DIR}/pbmc_10k_annotated_v3_vae.h5ad")
    print(f"Saved: pbmc_10k_annotated_v3_vae.h5ad")

    # Save model (UNSUPERVISED - no n_classes)
    torch.save({
        'model_state_dict': vae.state_dict(),
        'latent_dim': latent_dim,
        'hidden_dims': hidden_dims,
        'input_dim': input_dim,
        'random_seed': RANDOM_SEED,
    }, f"{OUTPUT_DIR}/vae_model.pt")
    print(f"Saved: vae_model.pt")

    print("\n" + "=" * 60)
    print("VAE embedding complete (UNSUPERVISED)!")
    print("=" * 60)
    print(f"\nNew embeddings stored in adata.obsm['X_vae'] ({latent_dim}D)")
    print(f"New UMAP stored in adata.obsm['X_umap_vae']")
    print(f"New clustering stored in adata.obs['leiden_vae'] ({n_clusters} clusters)")
    print(f"\nNOTE: VAE was trained WITHOUT cell type labels - purely unsupervised")
