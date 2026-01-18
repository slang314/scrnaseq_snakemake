#!/usr/bin/env python
"""
Downstream scRNA-seq analysis pipeline - v3 (FIXED)

Fixes implemented:
1. Cell calling: Knee-based filtering to remove empty droplets
2. Gene IDs: Keep Ensembl IDs (no version) as var_names, symbols in var['gene_symbol']
3. HVG → PCA: Consistent path using scanpy's built-in approach
4. Reproducibility: Fixed random seeds, version logging
5. VAE: Will be run separately as unsupervised (no circular semi-supervision)
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
import random
warnings.filterwarnings('ignore')

# =============================================================================
# REPRODUCIBILITY: Set all random seeds
# =============================================================================
RANDOM_SEED = 42
random.seed(RANDOM_SEED)
np.random.seed(RANDOM_SEED)

# Scanpy settings
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

# =============================================================================
# VERSION LOGGING
# =============================================================================
print("=" * 60)
print("SOFTWARE VERSIONS")
print("=" * 60)
import sys
print(f"Python: {sys.version}")
print(f"Scanpy: {sc.__version__}")
print(f"NumPy: {np.__version__}")
print(f"Pandas: {pd.__version__}")
try:
    import anndata
    print(f"AnnData: {anndata.__version__}")
except:
    pass
try:
    import scipy
    print(f"SciPy: {scipy.__version__}")
except:
    pass

# Paths
INPUT_H5AD = "results/anndata/qc/pbmc_10k_v3.h5ad"
GENE_MAP = "resources/reference/gene_id_to_name.tsv"
OUTPUT_DIR = "results/downstream_v3"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
print("\n" + "=" * 60)
print("1. Loading data")
print("=" * 60)

adata = sc.read_h5ad(INPUT_H5AD)
print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

# =============================================================================
# 2. CELL CALLING: KNEE-BASED FILTERING
# =============================================================================
print("\n" + "=" * 60)
print("2. Cell calling: Knee-based filtering")
print("=" * 60)

# The problem: We have ~30k barcodes from a "10k" dataset
# This means we're including empty droplets with ambient RNA

# Calculate total UMI counts per cell
if 'total_counts' not in adata.obs:
    adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()

# Sort cells by total counts (descending) for knee plot
sorted_counts = np.sort(adata.obs['total_counts'].values)[::-1]
ranks = np.arange(1, len(sorted_counts) + 1)

# Find knee using the "kneedle" algorithm approach
# The knee is where cumulative distance from diagonal is maximized
log_ranks = np.log10(ranks)
log_counts = np.log10(sorted_counts + 1)

# Normalize to [0, 1] range
x_norm = (log_ranks - log_ranks.min()) / (log_ranks.max() - log_ranks.min())
y_norm = (log_counts - log_counts.min()) / (log_counts.max() - log_counts.min())

# Distance from diagonal (line from first to last point)
# For a convex curve, the knee is where this distance is maximum
distances = y_norm - x_norm

# Focus on reasonable range for 10x data (rank 5000-15000 for 10k dataset)
# Real cells are typically in the first ~10k barcodes
search_start = 5000
search_end = min(15000, len(distances) - 1)
knee_idx = search_start + np.argmax(distances[search_start:search_end])
knee_threshold = sorted_counts[knee_idx]

# Sanity check: threshold should be at least 500 UMI for 10x v3 data
MIN_UMI_THRESHOLD = 500
if knee_threshold < MIN_UMI_THRESHOLD:
    print(f"Knee threshold ({knee_threshold:.0f}) too low, using minimum of {MIN_UMI_THRESHOLD}")
    # Find the rank where UMI drops below MIN_UMI_THRESHOLD
    knee_idx = np.searchsorted(-sorted_counts, -MIN_UMI_THRESHOLD)
    knee_threshold = MIN_UMI_THRESHOLD

print(f"Knee detected at rank {knee_idx}, UMI threshold = {knee_threshold:.0f}")

# Plot barcode rank plot with knee
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Log-log plot
ax1 = axes[0]
ax1.loglog(ranks, sorted_counts, 'b-', linewidth=0.5)
ax1.axhline(y=knee_threshold, color='r', linestyle='--', label=f'Knee threshold ({knee_threshold:.0f} UMI)')
ax1.axvline(x=knee_idx, color='r', linestyle=':', alpha=0.5)
ax1.set_xlabel('Barcode Rank')
ax1.set_ylabel('Total UMI Counts')
ax1.set_title('Barcode Rank Plot (log-log)')
ax1.legend()

# Linear zoom on knee region
ax2 = axes[1]
ax2.plot(ranks[:20000], sorted_counts[:20000], 'b-', linewidth=0.5)
ax2.axhline(y=knee_threshold, color='r', linestyle='--', label=f'Knee threshold ({knee_threshold:.0f} UMI)')
ax2.axvline(x=knee_idx, color='r', linestyle=':', alpha=0.5)
ax2.fill_between(ranks[:knee_idx], 0, sorted_counts[:knee_idx], alpha=0.3, color='green', label='Real cells')
ax2.set_xlabel('Barcode Rank')
ax2.set_ylabel('Total UMI Counts')
ax2.set_title('Barcode Rank Plot (linear, first 20k)')
ax2.legend()

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/barcode_rank_plot.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: barcode_rank_plot.png")

# Apply knee-based filtering
n_before = adata.n_obs
adata = adata[adata.obs['total_counts'] >= knee_threshold].copy()
n_after = adata.n_obs

print(f"\nCell filtering with UMI >= {knee_threshold:.0f}:")
print(f"  Before: {n_before} barcodes")
print(f"  After:  {n_after} cells")
print(f"  Removed: {n_before - n_after} empty droplets ({100*(n_before-n_after)/n_before:.1f}%)")

# =============================================================================
# 3. FIX GENE ID HANDLING
# =============================================================================
print("\n" + "=" * 60)
print("3. Fixing gene ID handling")
print("=" * 60)

# Current var_names are Ensembl IDs with version (e.g., ENSG00000141510.18)
# We need to:
# 1. Strip version suffix for stable IDs
# 2. Keep Ensembl ID as var_names (stable identifier)
# 3. Store gene symbols in var['gene_symbol']

# Strip version from Ensembl IDs
original_ids = adata.var_names.tolist()
stripped_ids = [g.split('.')[0] for g in original_ids]

# Check for duplicates after stripping
unique_stripped = len(set(stripped_ids))
print(f"Original IDs: {len(original_ids)}")
print(f"After version stripping: {unique_stripped} unique")

if unique_stripped < len(original_ids):
    print(f"WARNING: {len(original_ids) - unique_stripped} duplicate IDs after stripping versions")
    print("Keeping original versioned IDs to avoid ambiguity")
    adata.var['gene_id_no_version'] = stripped_ids
else:
    # Safe to use stripped IDs
    adata.var['gene_id_versioned'] = original_ids
    adata.var_names = stripped_ids
    print("Using Ensembl IDs without version as var_names")

# Load gene symbol mapping
gene_map_df = pd.read_csv(GENE_MAP, sep='\t', header=None, names=['gene_id', 'gene_symbol'])

# Create mapping dict (handle both versioned and unversioned)
symbol_map = {}
for _, row in gene_map_df.iterrows():
    gene_id = row['gene_id']
    symbol = row['gene_symbol']
    symbol_map[gene_id] = symbol
    # Also map without version
    symbol_map[gene_id.split('.')[0]] = symbol

# Map to symbols (store in var, NOT as var_names)
adata.var['gene_symbol'] = [symbol_map.get(g, None) for g in adata.var_names]
n_mapped = adata.var['gene_symbol'].notna().sum()
print(f"Mapped {n_mapped}/{adata.n_vars} genes to symbols")
print(f"Gene symbols stored in adata.var['gene_symbol']")
print(f"var_names remain as Ensembl IDs (stable identifiers)")

# Show example
print("\nExample var DataFrame:")
print(adata.var.head())

# =============================================================================
# 4. PRESERVE RAW COUNTS
# =============================================================================
print("\n" + "=" * 60)
print("4. Preserving raw counts")
print("=" * 60)

# Ensure total_counts is computed (needed for size_factors)
if 'total_counts' not in adata.obs:
    adata.obs['total_counts'] = np.array(adata.X.sum(axis=1)).flatten()

adata.layers['counts'] = adata.X.copy()
print(f"Raw counts preserved in adata.layers['counts']")

# =============================================================================
# 5. NORMALIZATION
# =============================================================================
print("\n" + "=" * 60)
print("5. Normalization")
print("=" * 60)

# Store size factors
adata.obs['size_factors'] = adata.obs['total_counts'].copy()

# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)
print("Normalized to 10k counts per cell")

# Log transform
sc.pp.log1p(adata)
print("Log1p transformed")

# =============================================================================
# 6. HVG SELECTION (FIXED)
# =============================================================================
print("\n" + "=" * 60)
print("6. HVG selection")
print("=" * 60)

# Find highly variable genes
# Using flavor='seurat' which doesn't require skmisc
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    subset=False  # Don't subset yet, just mark
)
n_hvg = adata.var['highly_variable'].sum()
print(f"Identified {n_hvg} highly variable genes")

# Store normalized data in .raw (all genes, for DE later)
adata.raw = adata
print("Stored normalized data in adata.raw (all genes)")

# =============================================================================
# 7. SCALING AND PCA (FIXED - CONSISTENT PATH)
# =============================================================================
print("\n" + "=" * 60)
print("7. Scaling and PCA (consistent HVG path)")
print("=" * 60)

# CRITICAL FIX: Use scanpy's built-in approach
# Scale only HVGs, run PCA on same data

# Subset to HVGs for dimensionality reduction
adata_hvg = adata[:, adata.var['highly_variable']].copy()
print(f"Subset to {adata_hvg.n_vars} HVGs")

# Scale the HVG data
sc.pp.scale(adata_hvg, max_value=10)
print("Scaled HVGs to unit variance (clipped at 10)")

# Run PCA on the SAME scaled HVG data
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack', random_state=RANDOM_SEED)
print("Computed 50 PCs on scaled HVGs")

# Copy PCA results back to main adata
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']

# Store which genes were used for PCA
adata.uns['pca']['highly_variable_genes'] = adata.var_names[adata.var['highly_variable']].tolist()

print(f"PCA variance explained (first 10): {adata.uns['pca']['variance_ratio'][:10].round(3)}")

# =============================================================================
# 8. NEIGHBORS AND UMAP
# =============================================================================
print("\n" + "=" * 60)
print("8. Neighbors and UMAP")
print("=" * 60)

# Build neighbor graph
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=RANDOM_SEED)
print("Built neighbor graph (k=15, 30 PCs)")

# UMAP
sc.tl.umap(adata, random_state=RANDOM_SEED)
print("Computed UMAP")

# =============================================================================
# 9. LEIDEN CLUSTERING
# =============================================================================
print("\n" + "=" * 60)
print("9. Leiden clustering")
print("=" * 60)

# Try multiple resolutions
for res in [0.2, 0.4, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}', random_state=RANDOM_SEED)
    n_clusters = adata.obs[f'leiden_{res}'].nunique()
    print(f"Leiden res={res}: {n_clusters} clusters")

# Use resolution 0.4 as default
adata.obs['leiden'] = adata.obs['leiden_0.4']
print(f"\nDefault: leiden_0.4 ({adata.obs['leiden'].nunique()} clusters)")

# =============================================================================
# 10. DOUBLET DETECTION
# =============================================================================
print("\n" + "=" * 60)
print("10. Doublet detection (Scrublet)")
print("=" * 60)

import scrublet as scr

# Run on raw counts
scrub = scr.Scrublet(
    adata.layers['counts'],
    expected_doublet_rate=0.06,
    random_state=RANDOM_SEED
)
doublet_scores, _ = scrub.scrub_doublets(
    min_counts=2,
    min_cells=3,
    min_gene_variability_pctl=85,
    n_prin_comps=30
)

# Use manual threshold - auto-threshold often too aggressive
# Standard 10x doublet rate is ~6% for 10k cells
DOUBLET_THRESHOLD = 0.25  # Manual threshold (auto-threshold is often too low)

adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = doublet_scores > DOUBLET_THRESHOLD

n_doublets = adata.obs['predicted_doublet'].sum()
auto_threshold = scrub.threshold_
print(f"Scrublet auto threshold: {auto_threshold:.3f} (not used)")
print(f"Manual threshold: {DOUBLET_THRESHOLD}")
print(f"Predicted doublets: {n_doublets} ({100*n_doublets/adata.n_obs:.1f}%)")

# Plot doublet scores
fig, ax = plt.subplots(figsize=(8, 4))
ax.hist(doublet_scores, bins=100, density=True, alpha=0.7)
ax.axvline(DOUBLET_THRESHOLD, color='red', linestyle='--', label=f'manual threshold={DOUBLET_THRESHOLD}')
ax.axvline(auto_threshold, color='orange', linestyle=':', alpha=0.5, label=f'auto threshold={auto_threshold:.2f}')
ax.set_xlabel('Doublet score')
ax.set_ylabel('Density')
ax.set_title('Scrublet doublet score distribution')
ax.legend()
plt.savefig(f"{OUTPUT_DIR}/figures/doublet_scores.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: doublet_scores.png")

# =============================================================================
# 11. VISUALIZATIONS
# =============================================================================
print("\n" + "=" * 60)
print("11. Visualizations")
print("=" * 60)

# UMAP with clusters
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

sc.pl.umap(adata, color='leiden', ax=axes[0, 0], show=False,
           title=f'Leiden clusters (n={adata.obs["leiden"].nunique()})', legend_loc='on data')
sc.pl.umap(adata, color='leiden_0.6', ax=axes[0, 1], show=False,
           title=f'Leiden res=0.6 (n={adata.obs["leiden_0.6"].nunique()})', legend_loc='on data')
sc.pl.umap(adata, color='doublet_score', ax=axes[1, 0], show=False,
           title='Doublet score', cmap='Reds')
sc.pl.umap(adata, color='total_counts', ax=axes[1, 1], show=False,
           title='Total UMI counts', cmap='viridis')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/umap_overview.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: umap_overview.png")

# =============================================================================
# 12. MARKER GENE DETECTION
# =============================================================================
print("\n" + "=" * 60)
print("12. Marker gene detection")
print("=" * 60)

# Find markers using raw (normalized, all genes)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='markers_leiden')
print("Computed marker genes per cluster (Wilcoxon)")

# Get marker dataframe and add gene symbols
markers_df = sc.get.rank_genes_groups_df(adata, group=None, key='markers_leiden')

# Map gene IDs to symbols in the markers table
markers_df['gene_symbol'] = markers_df['names'].map(
    lambda x: adata.var.loc[x, 'gene_symbol'] if x in adata.var.index else None
)

markers_df.to_csv(f"{OUTPUT_DIR}/marker_genes_all_clusters.tsv", sep='\t', index=False)
print("Saved: marker_genes_all_clusters.tsv")

# Top 5 markers per cluster (showing symbols)
print("\nTop 5 markers per cluster:")
top_markers = markers_df.groupby('group').head(5)
for cluster in sorted(adata.obs['leiden'].unique(), key=int):
    cluster_markers = top_markers[top_markers['group'] == cluster]
    symbols = cluster_markers['gene_symbol'].fillna(cluster_markers['names']).tolist()
    print(f"  Cluster {cluster}: {', '.join(symbols)}")

# =============================================================================
# 12b. CELL TYPE DIFFERENTIAL EXPRESSION (CellTypist annotations)
# =============================================================================
print("\n" + "=" * 60)
print("12b. Cell type differential expression")
print("=" * 60)

# Find DE genes between cell types using CellTypist annotations
sc.tl.rank_genes_groups(adata, groupby='celltypist_label', method='wilcoxon', key_added='markers_celltype')
print("Computed DE genes per cell type (Wilcoxon, one-vs-rest)")

# Get cell type marker dataframe
celltype_markers_df = sc.get.rank_genes_groups_df(adata, group=None, key='markers_celltype')

# Map gene IDs to symbols
celltype_markers_df['gene_symbol'] = celltype_markers_df['names'].map(
    lambda x: adata.var.loc[x, 'gene_symbol'] if x in adata.var.index else None
)

celltype_markers_df.to_csv(f"{OUTPUT_DIR}/de_genes_by_celltype.tsv", sep='\t', index=False)
print("Saved: de_genes_by_celltype.tsv")

# Top 5 DE genes per cell type
print("\nTop 5 DE genes per cell type:")
celltype_top = celltype_markers_df.groupby('group').head(5)
for celltype in sorted(adata.obs['celltypist_label'].unique()):
    ct_markers = celltype_top[celltype_top['group'] == celltype]
    symbols = ct_markers['gene_symbol'].fillna(ct_markers['names']).tolist()
    print(f"  {celltype}: {', '.join(symbols[:5])}")

# =============================================================================
# 13. PBMC MARKER VISUALIZATION
# =============================================================================
print("\n" + "=" * 60)
print("13. PBMC marker visualization")
print("=" * 60)

# Define canonical PBMC markers (gene symbols)
pbmc_markers = {
    'T cells': ['CD3D', 'CD3E', 'CD3G'],
    'CD4 T': ['CD4', 'IL7R'],
    'CD8 T': ['CD8A', 'CD8B'],
    'NK cells': ['GNLY', 'NKG7', 'KLRD1'],
    'B cells': ['CD19', 'MS4A1', 'CD79A'],
    'Monocytes': ['CD14', 'LYZ', 'S100A8', 'S100A9'],
    'CD16 Mono': ['FCGR3A', 'MS4A7'],
    'Dendritic': ['FCER1A', 'CST3'],
    'Platelets': ['PPBP', 'PF4'],
}

# Create symbol to ID mapping for lookup
symbol_to_id = adata.var.reset_index().set_index('gene_symbol')['index'].to_dict()

# Find available markers (need to look up by symbol)
available_markers = []
marker_ids = []
for cell_type, symbols in pbmc_markers.items():
    for symbol in symbols:
        if symbol in symbol_to_id:
            gene_id = symbol_to_id[symbol]
            if gene_id in adata.raw.var_names:
                available_markers.append(symbol)
                marker_ids.append(gene_id)

print(f"Found {len(available_markers)} PBMC markers: {available_markers}")

# Plot markers using gene IDs (with symbol labels)
if marker_ids:
    # Create a temporary view with symbol names for plotting
    fig, axes = plt.subplots(3, 3, figsize=(15, 15))
    key_markers = ['CD3D', 'CD14', 'MS4A1', 'GNLY', 'FCGR3A', 'CD4', 'CD8A', 'PPBP', 'FCER1A']

    for idx, symbol in enumerate(key_markers):
        if idx >= 9:
            break
        ax = axes.flatten()[idx]
        if symbol in symbol_to_id and symbol_to_id[symbol] in adata.raw.var_names:
            gene_id = symbol_to_id[symbol]
            sc.pl.umap(adata, color=gene_id, ax=ax, show=False, title=symbol)
        else:
            ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/figures/umap_pbmc_markers.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: umap_pbmc_markers.png")

# =============================================================================
# 14. CELLTYPIST ANNOTATION (POST-HOC ONLY)
# =============================================================================
print("\n" + "=" * 60)
print("14. Celltypist annotation (post-hoc)")
print("=" * 60)

try:
    import celltypist
    from celltypist import models

    # Create a copy with gene symbols as var_names for celltypist
    # Celltypist requires gene symbols
    adata_ct = adata.raw.to_adata().copy()

    # Map var_names to symbols
    adata_ct.var['original_id'] = adata_ct.var_names.tolist()
    new_names = []
    for gene_id in adata_ct.var_names:
        symbol = adata.var.loc[gene_id, 'gene_symbol'] if gene_id in adata.var.index else None
        new_names.append(symbol if pd.notna(symbol) else gene_id)
    adata_ct.var_names = new_names
    adata_ct.var_names_make_unique()

    # Download and load model
    models.download_models(force_update=False, model='Immune_All_Low.pkl')
    model = models.Model.load(model='Immune_All_Low.pkl')

    # Annotate
    predictions = celltypist.annotate(adata_ct, model=model, majority_voting=True)
    adata_pred = predictions.to_adata()

    # Transfer labels back to main adata
    adata.obs['celltypist_label'] = adata_pred.obs['majority_voting'].values
    adata.obs['celltypist_conf'] = adata_pred.obs['conf_score'].values

    print("Celltypist annotations:")
    print(adata.obs['celltypist_label'].value_counts())

    # Plot
    fig, ax = plt.subplots(figsize=(12, 10))
    sc.pl.umap(adata, color='celltypist_label', ax=ax, show=False, legend_fontsize=8)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/figures/umap_celltypist.png", dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved: umap_celltypist.png")

except Exception as e:
    print(f"Celltypist error: {e}")
    print("Continuing without automated annotation")

# =============================================================================
# 15. SAVE RESULTS
# =============================================================================
print("\n" + "=" * 60)
print("15. Saving results")
print("=" * 60)

# Remove doublets for clean version
adata_clean = adata[~adata.obs['predicted_doublet']].copy()
print(f"Clean dataset (doublets removed): {adata_clean.n_obs} cells")

# Save
adata.write_h5ad(f"{OUTPUT_DIR}/pbmc_10k_annotated_v3.h5ad")
print(f"Saved: pbmc_10k_annotated_v3.h5ad ({adata.n_obs} cells)")

adata_clean.write_h5ad(f"{OUTPUT_DIR}/pbmc_10k_annotated_v3_clean.h5ad")
print(f"Saved: pbmc_10k_annotated_v3_clean.h5ad ({adata_clean.n_obs} cells)")

# Summary
print("\n" + "=" * 60)
print("ANALYSIS SUMMARY")
print("=" * 60)
summary = {
    'cells_after_knee_filter': adata.n_obs,
    'genes': adata.n_vars,
    'hvgs': int(adata.var['highly_variable'].sum()),
    'n_clusters_res0.4': adata.obs['leiden'].nunique(),
    'predicted_doublets': int(adata.obs['predicted_doublet'].sum()),
    'cells_after_doublet_removal': adata_clean.n_obs,
    'random_seed': RANDOM_SEED,
}

for k, v in summary.items():
    print(f"  {k}: {v}")

pd.Series(summary).to_csv(f"{OUTPUT_DIR}/analysis_summary.tsv", sep='\t')

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
