#!/usr/bin/env python
"""
Downstream scRNA-seq analysis pipeline for PBMC 10k dataset.

Steps:
1. Load QC'd data, preserve raw counts
2. Normalize, log1p, HVG selection
3. PCA, neighbors, doublet detection (Scrublet)
4. UMAP, Leiden clustering
5. Marker gene detection
6. Cell type annotation
7. Pseudobulk aggregation
"""

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')

# Settings
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

# Paths
INPUT_H5AD = "results/anndata/qc/pbmc_10k_v3.h5ad"
OUTPUT_DIR = "results/downstream"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(f"{OUTPUT_DIR}/figures", exist_ok=True)

# =============================================================================
# 1. Load data and preserve raw counts
# =============================================================================
print("=" * 60)
print("1. Loading QC'd data and preserving raw counts")
print("=" * 60)

adata = sc.read_h5ad(INPUT_H5AD)
print(f"Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Preserve raw counts in layers BEFORE any normalization
adata.layers["counts"] = adata.X.copy()
print("Raw counts preserved in adata.layers['counts']")

# =============================================================================
# 2. Normalization, log1p, HVG selection
# =============================================================================
print("\n" + "=" * 60)
print("2. Normalization and HVG selection")
print("=" * 60)

# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)
print("Normalized to 10k counts per cell")

# Log transform
sc.pp.log1p(adata)
print("Log1p transformed")

# Find highly variable genes (using 'seurat' flavor which doesn't require skmisc)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5,
                            subset=False)
n_hvg = adata.var['highly_variable'].sum()
print(f"Identified {n_hvg} highly variable genes")

# Store normalized data in .raw (for DE later - uses all genes)
adata.raw = adata
print("Normalized data stored in adata.raw (all genes, for DE)")

# =============================================================================
# 3. Scale, PCA, neighbors
# =============================================================================
print("\n" + "=" * 60)
print("3. PCA and neighbor graph")
print("=" * 60)

# Subset to HVGs for PCA
adata_hvg = adata[:, adata.var['highly_variable']].copy()
print(f"Subset to {adata_hvg.shape[1]} HVGs for PCA")

# Scale
sc.pp.scale(adata_hvg, max_value=10)
print("Scaled to unit variance (clipped at 10)")

# PCA
sc.tl.pca(adata_hvg, n_comps=50, svd_solver='arpack')
print("Computed 50 PCs")

# Copy PCA results back to main object
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']
adata.varm['PCs'] = np.zeros((adata.shape[1], 50))
adata.varm['PCs'][adata.var['highly_variable'], :] = adata_hvg.varm['PCs']

# Neighbors
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
print("Built neighbor graph (k=15, 30 PCs)")

# =============================================================================
# 4. Doublet detection with Scrublet
# =============================================================================
print("\n" + "=" * 60)
print("4. Doublet detection (Scrublet)")
print("=" * 60)

try:
    import scrublet as scr

    # Run Scrublet on raw counts
    scrub = scr.Scrublet(adata.layers['counts'], expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(
        min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
    )

    adata.obs['doublet_score'] = doublet_scores
    adata.obs['predicted_doublet'] = predicted_doublets

    n_doublets = predicted_doublets.sum()
    print(f"Scrublet detected {n_doublets} predicted doublets ({100*n_doublets/len(predicted_doublets):.1f}%)")
    print(f"Doublet score range: {doublet_scores.min():.3f} - {doublet_scores.max():.3f}")
    print(f"Median doublet score: {np.median(doublet_scores):.3f}")

except ImportError:
    print("Scrublet not installed, using scanpy's scrublet implementation")
    sc.pp.scrublet(adata, batch_key=None, expected_doublet_rate=0.06)
    n_doublets = adata.obs['predicted_doublet'].sum()
    print(f"Detected {n_doublets} predicted doublets ({100*n_doublets/adata.n_obs:.1f}%)")

# =============================================================================
# 5. UMAP and Leiden clustering
# =============================================================================
print("\n" + "=" * 60)
print("5. UMAP and Leiden clustering")
print("=" * 60)

# UMAP
sc.tl.umap(adata, min_dist=0.3)
print("Computed UMAP")

# Leiden clustering at multiple resolutions
for res in [0.3, 0.5, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')
    n_clusters = adata.obs[f'leiden_{res}'].nunique()
    print(f"Leiden res={res}: {n_clusters} clusters")

# Use resolution 0.5 as default
adata.obs['leiden'] = adata.obs['leiden_0.5']
print(f"\nDefault clustering: leiden_0.5 ({adata.obs['leiden'].nunique()} clusters)")

# =============================================================================
# 6. Visualizations
# =============================================================================
print("\n" + "=" * 60)
print("6. Creating visualizations")
print("=" * 60)

# UMAP colored by clusters
fig, axes = plt.subplots(2, 2, figsize=(14, 12))

sc.pl.umap(adata, color='leiden', ax=axes[0, 0], show=False,
           title='Leiden clusters (res=0.5)', legend_loc='on data')
sc.pl.umap(adata, color='leiden_0.8', ax=axes[0, 1], show=False,
           title='Leiden clusters (res=0.8)', legend_loc='on data')
sc.pl.umap(adata, color='doublet_score', ax=axes[1, 0], show=False,
           title='Doublet score', cmap='Reds')
sc.pl.umap(adata, color='predicted_doublet', ax=axes[1, 1], show=False,
           title='Predicted doublets')

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/umap_clusters_doublets.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: umap_clusters_doublets.png")

# QC metrics on UMAP
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
sc.pl.umap(adata, color='n_genes_by_counts', ax=axes[0], show=False, cmap='viridis')
sc.pl.umap(adata, color='total_counts', ax=axes[1], show=False, cmap='viridis')
if 'pct_counts_mt' in adata.obs.columns:
    sc.pl.umap(adata, color='pct_counts_mt', ax=axes[2], show=False, cmap='Reds')
else:
    axes[2].set_visible(False)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/figures/umap_qc_metrics.png", dpi=150, bbox_inches='tight')
plt.close()
print("Saved: umap_qc_metrics.png")

# Doublet scores per cluster
doublet_by_cluster = adata.obs.groupby('leiden').agg({
    'doublet_score': ['mean', 'median'],
    'predicted_doublet': 'sum'
}).round(3)
doublet_by_cluster.columns = ['mean_doublet_score', 'median_doublet_score', 'n_predicted_doublets']
doublet_by_cluster['n_cells'] = adata.obs['leiden'].value_counts().sort_index()
doublet_by_cluster['pct_doublets'] = (100 * doublet_by_cluster['n_predicted_doublets'] /
                                       doublet_by_cluster['n_cells']).round(1)
print("\nDoublet scores by cluster:")
print(doublet_by_cluster.to_string())

doublet_by_cluster.to_csv(f"{OUTPUT_DIR}/doublet_scores_by_cluster.tsv", sep='\t')

# Flag clusters with high doublet enrichment (>20%)
high_doublet_clusters = doublet_by_cluster[doublet_by_cluster['pct_doublets'] > 20].index.tolist()
if high_doublet_clusters:
    print(f"\nWARNING: Clusters with >20% doublets: {high_doublet_clusters}")
    print("Consider removing these before annotation.")

# =============================================================================
# 7. Marker gene detection
# =============================================================================
print("\n" + "=" * 60)
print("7. Marker gene detection")
print("=" * 60)

# Find markers for each cluster (uses .raw which has all genes, normalized)
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon',
                        key_added='markers_leiden')
print("Computed marker genes per cluster (Wilcoxon)")

# Extract top markers
markers_df = sc.get.rank_genes_groups_df(adata, group=None, key='markers_leiden')
markers_df.to_csv(f"{OUTPUT_DIR}/marker_genes_all_clusters.tsv", sep='\t', index=False)
print(f"Saved all markers to marker_genes_all_clusters.tsv")

# Top 5 markers per cluster
top_markers = markers_df.groupby('group').head(5)
print("\nTop 5 markers per cluster:")
for cluster in sorted(adata.obs['leiden'].unique(), key=int):
    genes = top_markers[top_markers['group'] == cluster]['names'].tolist()
    print(f"  Cluster {cluster}: {', '.join(genes)}")

# Plot top markers
sc.pl.rank_genes_groups(adata, n_genes=10, key='markers_leiden',
                        save='_top10_per_cluster.png', show=False)
# Move scanpy's auto-saved figure
import shutil
if os.path.exists("figures/rank_genes_groups_leiden_top10_per_cluster.png"):
    shutil.move("figures/rank_genes_groups_leiden_top10_per_cluster.png",
                f"{OUTPUT_DIR}/figures/rank_genes_groups_top10.png")
    print("Saved: rank_genes_groups_top10.png")

# =============================================================================
# 8. PBMC marker genes visualization
# =============================================================================
print("\n" + "=" * 60)
print("8. PBMC canonical marker genes")
print("=" * 60)

# Define canonical PBMC markers
pbmc_markers = {
    'T cells': ['CD3D', 'CD3E', 'CD3G'],
    'CD4 T': ['CD4', 'IL7R'],
    'CD8 T': ['CD8A', 'CD8B'],
    'NK cells': ['GNLY', 'NKG7', 'KLRD1', 'NCAM1'],
    'B cells': ['CD19', 'MS4A1', 'CD79A'],
    'Monocytes': ['CD14', 'LYZ', 'S100A8', 'S100A9'],
    'CD16 Mono': ['FCGR3A', 'MS4A7'],
    'Dendritic': ['FCER1A', 'CST3', 'CLEC10A'],
    'Platelets': ['PPBP', 'PF4'],
}

# Check which markers are present
available_markers = []
marker_to_type = {}
for cell_type, genes in pbmc_markers.items():
    for gene in genes:
        if gene in adata.var_names:
            available_markers.append(gene)
            marker_to_type[gene] = cell_type
        elif gene in adata.raw.var_names:
            available_markers.append(gene)
            marker_to_type[gene] = cell_type

print(f"Found {len(available_markers)} PBMC markers in dataset")

# Dotplot of markers
if available_markers:
    sc.pl.dotplot(adata, var_names=available_markers, groupby='leiden',
                  standard_scale='var', save='_pbmc_markers.png', show=False)
    if os.path.exists("figures/dotplot_pbmc_markers.png"):
        shutil.move("figures/dotplot_pbmc_markers.png",
                    f"{OUTPUT_DIR}/figures/dotplot_pbmc_markers.png")
        print("Saved: dotplot_pbmc_markers.png")

    # UMAP with key markers
    key_markers = ['CD3D', 'CD14', 'MS4A1', 'GNLY', 'FCGR3A', 'FCER1A']
    key_markers = [m for m in key_markers if m in adata.raw.var_names]

    if key_markers:
        sc.pl.umap(adata, color=key_markers, ncols=3, save='_key_markers.png', show=False)
        if os.path.exists("figures/umap_key_markers.png"):
            shutil.move("figures/umap_key_markers.png",
                        f"{OUTPUT_DIR}/figures/umap_key_markers.png")
            print("Saved: umap_key_markers.png")

# =============================================================================
# 9. Manual cell type annotation based on markers
# =============================================================================
print("\n" + "=" * 60)
print("9. Cell type annotation")
print("=" * 60)

# Score cells for each marker set
for cell_type, genes in pbmc_markers.items():
    genes_present = [g for g in genes if g in adata.raw.var_names]
    if genes_present:
        sc.tl.score_genes(adata, gene_list=genes_present,
                         score_name=f'score_{cell_type.replace(" ", "_")}')

# Create annotation based on marker expression
# We'll use a simple approach: check mean expression of key markers per cluster
print("\nCluster annotation based on marker expression:")

cluster_annotations = {}
for cluster in sorted(adata.obs['leiden'].unique(), key=int):
    mask = adata.obs['leiden'] == cluster
    cluster_data = adata[mask]

    # Calculate mean expression of key markers
    scores = {}
    for cell_type, genes in pbmc_markers.items():
        genes_present = [g for g in genes if g in adata.raw.var_names]
        if genes_present:
            # Get expression from raw (normalized, log1p)
            expr = adata.raw[mask, genes_present].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray()
            scores[cell_type] = np.mean(expr)

    # Assign most likely cell type
    if scores:
        best_type = max(scores, key=scores.get)
        cluster_annotations[cluster] = best_type
        print(f"  Cluster {cluster}: {best_type} (score: {scores[best_type]:.2f})")

# Add annotations to adata
adata.obs['cell_type_manual'] = adata.obs['leiden'].map(cluster_annotations)
print(f"\nCell type distribution:")
print(adata.obs['cell_type_manual'].value_counts())

# =============================================================================
# 10. Automated annotation with celltypist (if available)
# =============================================================================
print("\n" + "=" * 60)
print("10. Automated annotation (celltypist)")
print("=" * 60)

try:
    import celltypist
    from celltypist import models

    # Download immune cell model if not present
    models.download_models(force_update=False, model='Immune_All_Low.pkl')
    model = models.Model.load(model='Immune_All_Low.pkl')

    # Run prediction
    predictions = celltypist.annotate(adata, model=model, majority_voting=True)
    adata_ct = predictions.to_adata()

    # Add predictions to our adata
    adata.obs['celltypist_label'] = adata_ct.obs['majority_voting']
    adata.obs['celltypist_conf'] = adata_ct.obs['conf_score']

    print("Celltypist annotations added:")
    print(adata.obs['celltypist_label'].value_counts().head(15))

    # Plot celltypist results
    sc.pl.umap(adata, color='celltypist_label', legend_loc='on data',
               save='_celltypist.png', show=False)
    if os.path.exists("figures/umap_celltypist.png"):
        shutil.move("figures/umap_celltypist.png",
                    f"{OUTPUT_DIR}/figures/umap_celltypist.png")
        print("Saved: umap_celltypist.png")

except ImportError:
    print("Celltypist not installed. Skipping automated annotation.")
    print("Install with: pip install celltypist")
except Exception as e:
    print(f"Celltypist error: {e}")
    print("Continuing with manual annotations only.")

# =============================================================================
# 11. Remove/flag doublets
# =============================================================================
print("\n" + "=" * 60)
print("11. Doublet handling")
print("=" * 60)

n_doublets = adata.obs['predicted_doublet'].sum()
print(f"Total predicted doublets: {n_doublets} ({100*n_doublets/adata.n_obs:.1f}%)")

# Create a filtered version without doublets
adata_clean = adata[~adata.obs['predicted_doublet']].copy()
print(f"After removing doublets: {adata_clean.n_obs} cells")

# =============================================================================
# 12. Pseudobulk aggregation
# =============================================================================
print("\n" + "=" * 60)
print("12. Pseudobulk aggregation")
print("=" * 60)

# Aggregate counts by cell type (using clean data)
def pseudobulk_by_group(adata, groupby, layer='counts'):
    """Aggregate raw counts by group."""
    groups = adata.obs[groupby].unique()
    pseudobulk = {}

    for group in groups:
        mask = adata.obs[groupby] == group
        counts = adata[mask].layers[layer]
        if hasattr(counts, 'toarray'):
            counts = counts.toarray()
        pseudobulk[group] = counts.sum(axis=0)

    pb_df = pd.DataFrame(pseudobulk, index=adata.var_names)
    return pb_df

# Pseudobulk by manual cell type annotation
pb_celltype = pseudobulk_by_group(adata_clean, 'cell_type_manual')
pb_celltype.to_csv(f"{OUTPUT_DIR}/pseudobulk_by_celltype.tsv", sep='\t')
print(f"Saved pseudobulk by cell type: {pb_celltype.shape}")

# Pseudobulk by cluster
pb_cluster = pseudobulk_by_group(adata_clean, 'leiden')
pb_cluster.to_csv(f"{OUTPUT_DIR}/pseudobulk_by_cluster.tsv", sep='\t')
print(f"Saved pseudobulk by cluster: {pb_cluster.shape}")

# Cell counts per group
cell_counts = adata_clean.obs['cell_type_manual'].value_counts()
cell_counts.to_csv(f"{OUTPUT_DIR}/cell_counts_by_celltype.tsv", sep='\t')
print("\nCell counts by cell type (after doublet removal):")
print(cell_counts)

# =============================================================================
# 13. Save final annotated object
# =============================================================================
print("\n" + "=" * 60)
print("13. Saving results")
print("=" * 60)

# Save full annotated object (with doublets flagged)
adata.write_h5ad(f"{OUTPUT_DIR}/pbmc_10k_annotated.h5ad")
print(f"Saved: pbmc_10k_annotated.h5ad ({adata.n_obs} cells)")

# Save clean object (doublets removed)
adata_clean.write_h5ad(f"{OUTPUT_DIR}/pbmc_10k_annotated_clean.h5ad")
print(f"Saved: pbmc_10k_annotated_clean.h5ad ({adata_clean.n_obs} cells)")

# Summary statistics
summary = {
    'total_cells_after_qc': adata.n_obs,
    'total_genes': adata.n_vars,
    'n_hvgs': adata.var['highly_variable'].sum(),
    'n_clusters': adata.obs['leiden'].nunique(),
    'n_predicted_doublets': int(adata.obs['predicted_doublet'].sum()),
    'cells_after_doublet_removal': adata_clean.n_obs,
}
summary_df = pd.Series(summary)
summary_df.to_csv(f"{OUTPUT_DIR}/analysis_summary.tsv", sep='\t', header=False)
print("\nAnalysis Summary:")
for k, v in summary.items():
    print(f"  {k}: {v}")

print("\n" + "=" * 60)
print("Analysis complete!")
print("=" * 60)
print(f"\nOutputs saved to: {OUTPUT_DIR}/")
print("  - pbmc_10k_annotated.h5ad (full annotated object)")
print("  - pbmc_10k_annotated_clean.h5ad (doublets removed)")
print("  - pseudobulk_by_celltype.tsv")
print("  - pseudobulk_by_cluster.tsv")
print("  - marker_genes_all_clusters.tsv")
print("  - figures/")
