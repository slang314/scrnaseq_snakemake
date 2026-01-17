#!/usr/bin/env python
"""
Scanpy QC script for scRNA-seq data.

Performs standard quality control:
1. Calculate QC metrics (gene counts, UMI counts, mitochondrial %)
2. Generate QC plots
3. Filter cells and genes based on thresholds
4. Save filtered AnnData and metrics
"""
import argparse
from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for headless systems
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser(description="Scanpy QC for scRNA-seq")
    ap.add_argument("--input", dest="inp", required=True, help="Input h5ad file")
    ap.add_argument("--out", required=True, help="Output QC'd h5ad file")
    ap.add_argument("--metrics", required=True, help="Output metrics TSV file")
    ap.add_argument("--mito-prefix", required=True, help="Mitochondrial gene prefix (e.g., MT- for human)")
    ap.add_argument("--min-genes", type=int, required=True, help="Minimum genes per cell")
    ap.add_argument("--min-cells-per-gene", type=int, required=True, help="Minimum cells per gene")
    ap.add_argument("--max-pct-mt", type=float, required=True, help="Maximum mitochondrial percentage")
    ap.add_argument("--max-counts", type=float, default=None, help="Maximum total counts (optional)")
    ap.add_argument("--max-genes", type=float, default=None, help="Maximum genes per cell (optional)")
    args = ap.parse_args()

    print(f"Reading {args.inp}")
    adata = sc.read_h5ad(args.inp)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Get sample name from obs if available
    sample_name = "unknown"
    if 'sample' in adata.obs.columns:
        sample_name = adata.obs['sample'].iloc[0]

    # Determine gene symbols column
    if 'gene_name' in adata.var.columns:
        gene_symbols = adata.var['gene_name'].astype(str)
    elif 'gene_names' in adata.var.columns:
        gene_symbols = adata.var['gene_names'].astype(str)
    else:
        gene_symbols = adata.var.index.astype(str)

    # Annotate mitochondrial genes
    adata.var['mt'] = gene_symbols.str.upper().str.startswith(args.mito_prefix.upper())
    n_mt_genes = adata.var['mt'].sum()
    print(f"Found {n_mt_genes} mitochondrial genes with prefix '{args.mito_prefix}'")

    # Calculate QC metrics
    print("Calculating QC metrics...")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Create output directories
    out_dir = Path(args.metrics).parent
    out_dir.mkdir(parents=True, exist_ok=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)

    # Configure scanpy plotting
    sc.settings.figdir = str(out_dir)
    sc.settings.set_figure_params(dpi=100, facecolor='white')

    # Generate QC plots
    print("Generating QC plots...")

    # Violin plots
    try:
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        sc.pl.violin(adata, 'n_genes_by_counts', ax=axes[0], show=False)
        sc.pl.violin(adata, 'total_counts', ax=axes[1], show=False)
        sc.pl.violin(adata, 'pct_counts_mt', ax=axes[2], show=False)
        plt.tight_layout()
        fig.savefig(out_dir / f"{sample_name}_qc_violin.png", dpi=100, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"Warning: Could not generate violin plot: {e}")

    # Scatter plot
    try:
        fig, ax = plt.subplots(figsize=(6, 5))
        sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', color='pct_counts_mt',
                      ax=ax, show=False)
        fig.savefig(out_dir / f"{sample_name}_qc_scatter.png", dpi=100, bbox_inches='tight')
        plt.close(fig)
    except Exception as e:
        print(f"Warning: Could not generate scatter plot: {e}")

    # Record pre-filter stats
    n_cells_before = adata.n_obs
    n_genes_before = adata.n_vars

    # Filtering
    print("Applying filters...")

    # 1. Filter genes: keep genes expressed in at least min_cells_per_gene cells
    sc.pp.filter_genes(adata, min_cells=args.min_cells_per_gene)
    n_genes_after_gene_filter = adata.n_vars
    print(f"  Genes after min_cells filter: {n_genes_after_gene_filter}")

    # 2. Filter cells: keep cells with at least min_genes genes
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    n_cells_after_min_genes = adata.n_obs
    print(f"  Cells after min_genes filter: {n_cells_after_min_genes}")

    # 3. Apply mitochondrial threshold
    keep = adata.obs['pct_counts_mt'] < args.max_pct_mt

    # 4. Optional upper bounds
    if args.max_counts is not None:
        keep = keep & (adata.obs['total_counts'] <= args.max_counts)
        print(f"  Applying max_counts filter: {args.max_counts}")

    if args.max_genes is not None:
        keep = keep & (adata.obs['n_genes_by_counts'] <= args.max_genes)
        print(f"  Applying max_genes filter: {args.max_genes}")

    adata = adata[keep].copy()

    n_cells_after = adata.n_obs
    n_genes_after = adata.n_vars

    print(f"After QC: {n_cells_after} cells x {n_genes_after} genes")
    print(f"Removed {n_cells_before - n_cells_after} cells ({100*(n_cells_before - n_cells_after)/n_cells_before:.1f}%)")

    # Write metrics
    metrics = {
        'sample': sample_name,
        'input_h5ad': args.inp,
        'output_h5ad': args.out,
        'min_genes': args.min_genes,
        'min_cells_per_gene': args.min_cells_per_gene,
        'max_pct_mt': args.max_pct_mt,
        'max_counts': args.max_counts if args.max_counts else 'NA',
        'max_genes': args.max_genes if args.max_genes else 'NA',
        'cells_before': n_cells_before,
        'genes_before': n_genes_before,
        'cells_after': n_cells_after,
        'genes_after': n_genes_after,
        'cells_removed': n_cells_before - n_cells_after,
        'pct_cells_removed': round(100 * (n_cells_before - n_cells_after) / n_cells_before, 2) if n_cells_before > 0 else 0,
        'median_genes_per_cell': int(adata.obs['n_genes_by_counts'].median()) if n_cells_after > 0 else 0,
        'median_counts_per_cell': int(adata.obs['total_counts'].median()) if n_cells_after > 0 else 0,
        'median_pct_mt': round(adata.obs['pct_counts_mt'].median(), 2) if n_cells_after > 0 else 0,
    }

    pd.DataFrame([metrics]).to_csv(args.metrics, sep='\t', index=False)
    print(f"Metrics saved to {args.metrics}")

    # Store QC info in adata
    adata.uns['qc'] = metrics

    # Write output
    adata.write_h5ad(args.out)
    print(f"QC'd AnnData saved to {args.out}")
    print("Done!")


if __name__ == "__main__":
    main()
