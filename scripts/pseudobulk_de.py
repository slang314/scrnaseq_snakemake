#!/usr/bin/env python
"""
Pseudobulk Differential Expression Analysis

This script:
1. Loads annotated scRNA-seq data (from downstream_analysis_v3.py or Cell Ranger)
2. Aggregates cells into pseudobulk samples by cell type (+ sample if multi-sample)
3. Performs pairwise differential expression between all cell types
4. Outputs DE results and visualizations

Requires: pydeseq2, scanpy, pandas, numpy, matplotlib, seaborn
"""

import argparse
import os
import sys
import warnings
from itertools import combinations
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse

warnings.filterwarnings('ignore')

# =============================================================================
# CONFIGURATION
# =============================================================================

DEFAULT_INPUT = "results/downstream_v3/pbmc_10k_annotated_v3_clean.h5ad"
DEFAULT_OUTPUT_DIR = "results/pseudobulk"
DEFAULT_ANNOTATION_COL = "celltypist_label"  # Column containing cell type annotations
DEFAULT_SAMPLE_COL = "sample"  # Column containing sample IDs (for multi-sample)

# DE thresholds
PADJ_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0  # |log2FC| > 1 means 2-fold change
MIN_CELLS_PER_GROUP = 10  # Minimum cells to form a pseudobulk sample


# =============================================================================
# PSEUDOBULK AGGREGATION
# =============================================================================

def create_pseudobulk(adata, groupby_cols, layer='counts', min_cells=10):
    """
    Aggregate single-cell counts into pseudobulk samples.

    Parameters
    ----------
    adata : AnnData
        Annotated single-cell data with raw counts
    groupby_cols : list
        Columns to group by (e.g., ['sample', 'celltypist_label'])
    layer : str
        Layer containing raw counts (default: 'counts')
    min_cells : int
        Minimum cells required to form a pseudobulk sample

    Returns
    -------
    pseudobulk_df : DataFrame
        Pseudobulk count matrix (samples × genes)
    metadata_df : DataFrame
        Sample metadata
    """
    print(f"Creating pseudobulk by grouping: {groupby_cols}")

    # Get raw counts
    if layer and layer in adata.layers:
        X = adata.layers[layer]
        print(f"  Using counts from layer '{layer}'")
    else:
        X = adata.X
        print(f"  Using adata.X (assuming raw counts)")

    # Convert to dense if sparse
    if sparse.issparse(X):
        X = X.toarray()

    # Create group labels
    if len(groupby_cols) == 1:
        groups = adata.obs[groupby_cols[0]].astype(str)
    else:
        groups = adata.obs[groupby_cols].astype(str).agg('_'.join, axis=1)

    # Get unique groups
    unique_groups = groups.unique()
    print(f"  Found {len(unique_groups)} pseudobulk groups")

    # Aggregate counts per group
    pseudobulk_data = []
    metadata = []

    for group in unique_groups:
        mask = groups == group
        n_cells = mask.sum()

        if n_cells < min_cells:
            print(f"  Skipping '{group}': only {n_cells} cells (min={min_cells})")
            continue

        # Sum counts across cells in this group
        group_counts = X[mask].sum(axis=0)
        pseudobulk_data.append(group_counts)

        # Extract metadata
        meta = {'pseudobulk_id': group, 'n_cells': n_cells}
        for col in groupby_cols:
            meta[col] = adata.obs.loc[mask, col].iloc[0]
        metadata.append(meta)

    # Create DataFrames
    pseudobulk_df = pd.DataFrame(
        np.array(pseudobulk_data),
        index=[m['pseudobulk_id'] for m in metadata],
        columns=adata.var_names
    )
    metadata_df = pd.DataFrame(metadata).set_index('pseudobulk_id')

    print(f"  Created {len(pseudobulk_df)} pseudobulk samples")
    print(f"  Total cells included: {metadata_df['n_cells'].sum()}")

    return pseudobulk_df, metadata_df


# =============================================================================
# DIFFERENTIAL EXPRESSION
# =============================================================================

def run_deseq2(counts_df, metadata_df, group_col, group1, group2):
    """
    Run DESeq2-style differential expression between two groups.

    Uses pyDESeq2 for proper negative binomial modeling.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("ERROR: pydeseq2 not installed. Install with: pip install pydeseq2")
        sys.exit(1)

    # Subset to the two groups
    mask = metadata_df[group_col].isin([group1, group2])
    counts_subset = counts_df.loc[mask].T  # DESeq2 expects genes × samples
    meta_subset = metadata_df.loc[mask].copy()

    # Ensure counts are integers
    counts_subset = counts_subset.astype(int)

    # Create condition column for DE
    meta_subset['condition'] = meta_subset[group_col].map({group1: 'A', group2: 'B'})

    # Filter low-count genes
    min_counts = 10
    keep_genes = counts_subset.sum(axis=1) >= min_counts
    counts_subset = counts_subset.loc[keep_genes]

    if counts_subset.shape[0] < 100:
        print(f"  Warning: Only {counts_subset.shape[0]} genes passed filter")

    # Run DESeq2
    dds = DeseqDataSet(
        counts=counts_subset,
        metadata=meta_subset,
        design_factors="condition"
    )

    dds.deseq2()

    # Get results (B vs A, so positive log2FC means higher in group2)
    stat_res = DeseqStats(dds, contrast=["condition", "B", "A"])
    stat_res.summary()

    results_df = stat_res.results_df.copy()
    results_df['gene'] = results_df.index
    results_df['comparison'] = f"{group2}_vs_{group1}"
    results_df['group1'] = group1
    results_df['group2'] = group2

    return results_df


def run_simple_de(counts_df, metadata_df, group_col, group1, group2):
    """
    Simple DE using Wilcoxon rank-sum test (fallback if pyDESeq2 unavailable).
    Less statistically rigorous but works for basic comparisons.
    """
    from scipy.stats import ranksums

    mask1 = metadata_df[group_col] == group1
    mask2 = metadata_df[group_col] == group2

    counts1 = counts_df.loc[mask1]
    counts2 = counts_df.loc[mask2]

    # Normalize by library size for comparison
    norm1 = counts1.div(counts1.sum(axis=1), axis=0) * 1e6  # CPM
    norm2 = counts2.div(counts2.sum(axis=1), axis=0) * 1e6

    results = []
    for gene in counts_df.columns:
        vals1 = norm1[gene].values
        vals2 = norm2[gene].values

        # Skip genes with no expression
        if vals1.sum() == 0 and vals2.sum() == 0:
            continue

        # Wilcoxon test
        stat, pval = ranksums(vals1, vals2)

        # Log2 fold change (mean CPM)
        mean1 = vals1.mean() + 1  # Add pseudocount
        mean2 = vals2.mean() + 1
        log2fc = np.log2(mean2 / mean1)

        results.append({
            'gene': gene,
            'log2FoldChange': log2fc,
            'pvalue': pval,
            'baseMean': (mean1 + mean2) / 2,
            'group1': group1,
            'group2': group2,
            'comparison': f"{group2}_vs_{group1}"
        })

    results_df = pd.DataFrame(results)

    # FDR correction
    from scipy.stats import false_discovery_control
    if len(results_df) > 0:
        results_df['padj'] = false_discovery_control(results_df['pvalue'].values, method='bh')

    return results_df


def run_pairwise_de(counts_df, metadata_df, group_col, method='deseq2'):
    """
    Run DE for all pairwise combinations of groups.
    """
    groups = metadata_df[group_col].unique()
    pairs = list(combinations(sorted(groups), 2))

    print(f"\nRunning {len(pairs)} pairwise comparisons...")

    all_results = []

    for i, (group1, group2) in enumerate(pairs):
        print(f"\n[{i+1}/{len(pairs)}] {group2} vs {group1}")

        # Check sample sizes
        n1 = (metadata_df[group_col] == group1).sum()
        n2 = (metadata_df[group_col] == group2).sum()
        print(f"  n={n1} vs n={n2}")

        if n1 < 2 or n2 < 2:
            print(f"  Skipping: need at least 2 samples per group")
            continue

        try:
            if method == 'deseq2':
                results = run_deseq2(counts_df, metadata_df, group_col, group1, group2)
            else:
                results = run_simple_de(counts_df, metadata_df, group_col, group1, group2)

            # Count significant genes
            sig = results[(results['padj'] < PADJ_THRESHOLD) &
                         (results['log2FoldChange'].abs() > LOG2FC_THRESHOLD)]
            print(f"  Significant genes: {len(sig)} (padj<{PADJ_THRESHOLD}, |log2FC|>{LOG2FC_THRESHOLD})")

            all_results.append(results)

        except Exception as e:
            print(f"  Error: {e}")
            continue

    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        return combined_df
    else:
        return pd.DataFrame()


# =============================================================================
# VISUALIZATION
# =============================================================================

def plot_volcano(results_df, comparison, output_dir, gene_symbols=None):
    """Create volcano plot for a single comparison."""
    df = results_df[results_df['comparison'] == comparison].copy()

    if len(df) == 0:
        return

    # Add significance category
    df['significant'] = 'NS'
    df.loc[(df['padj'] < PADJ_THRESHOLD) & (df['log2FoldChange'] > LOG2FC_THRESHOLD), 'significant'] = 'Up'
    df.loc[(df['padj'] < PADJ_THRESHOLD) & (df['log2FoldChange'] < -LOG2FC_THRESHOLD), 'significant'] = 'Down'

    # Plot
    fig, ax = plt.subplots(figsize=(8, 6))

    colors = {'NS': 'lightgrey', 'Up': 'red', 'Down': 'blue'}

    for sig_type in ['NS', 'Down', 'Up']:
        subset = df[df['significant'] == sig_type]
        ax.scatter(
            subset['log2FoldChange'],
            -np.log10(subset['padj'] + 1e-300),
            c=colors[sig_type],
            alpha=0.5,
            s=10,
            label=f"{sig_type} ({len(subset)})"
        )

    # Add threshold lines
    ax.axhline(-np.log10(PADJ_THRESHOLD), color='grey', linestyle='--', linewidth=0.5)
    ax.axvline(LOG2FC_THRESHOLD, color='grey', linestyle='--', linewidth=0.5)
    ax.axvline(-LOG2FC_THRESHOLD, color='grey', linestyle='--', linewidth=0.5)

    # Label top genes
    top_genes = df.nsmallest(10, 'padj')
    for _, row in top_genes.iterrows():
        gene_label = row['gene']
        if gene_symbols is not None and row['gene'] in gene_symbols:
            gene_label = gene_symbols[row['gene']]
        ax.annotate(
            gene_label,
            (row['log2FoldChange'], -np.log10(row['padj'] + 1e-300)),
            fontsize=6,
            alpha=0.7
        )

    ax.set_xlabel('log2 Fold Change')
    ax.set_ylabel('-log10(adjusted p-value)')
    ax.set_title(comparison.replace('_', ' '))
    ax.legend(loc='upper right')

    plt.tight_layout()
    safe_name = comparison.replace('/', '_').replace(' ', '_')
    plt.savefig(f"{output_dir}/volcano_{safe_name}.png", dpi=150, bbox_inches='tight')
    plt.close()


def plot_de_summary(results_df, output_dir):
    """Create summary heatmap of DE gene counts."""
    # Count significant genes per comparison
    sig_counts = results_df[
        (results_df['padj'] < PADJ_THRESHOLD) &
        (results_df['log2FoldChange'].abs() > LOG2FC_THRESHOLD)
    ].groupby('comparison').size()

    up_counts = results_df[
        (results_df['padj'] < PADJ_THRESHOLD) &
        (results_df['log2FoldChange'] > LOG2FC_THRESHOLD)
    ].groupby('comparison').size()

    down_counts = results_df[
        (results_df['padj'] < PADJ_THRESHOLD) &
        (results_df['log2FoldChange'] < -LOG2FC_THRESHOLD)
    ].groupby('comparison').size()

    summary_df = pd.DataFrame({
        'Up': up_counts,
        'Down': down_counts,
        'Total': sig_counts
    }).fillna(0).astype(int)

    # Plot
    fig, ax = plt.subplots(figsize=(10, max(6, len(summary_df) * 0.4)))

    summary_df[['Up', 'Down']].plot(
        kind='barh',
        stacked=True,
        color=['red', 'blue'],
        ax=ax,
        alpha=0.7
    )

    ax.set_xlabel('Number of DE genes')
    ax.set_title(f'Differential Expression Summary\n(padj<{PADJ_THRESHOLD}, |log2FC|>{LOG2FC_THRESHOLD})')
    ax.legend(title='Direction')

    plt.tight_layout()
    plt.savefig(f"{output_dir}/de_summary.png", dpi=150, bbox_inches='tight')
    plt.close()

    return summary_df


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Pseudobulk differential expression analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with defaults
  python pseudobulk_de.py

  # Custom input and annotation column
  python pseudobulk_de.py --input my_data.h5ad --annotation leiden

  # Multi-sample analysis (group by sample + cell type)
  python pseudobulk_de.py --multi-sample
        """
    )

    parser.add_argument('--input', '-i', default=DEFAULT_INPUT,
                        help=f'Input h5ad file (default: {DEFAULT_INPUT})')
    parser.add_argument('--output', '-o', default=DEFAULT_OUTPUT_DIR,
                        help=f'Output directory (default: {DEFAULT_OUTPUT_DIR})')
    parser.add_argument('--annotation', '-a', default=DEFAULT_ANNOTATION_COL,
                        help=f'Column with cell type annotations (default: {DEFAULT_ANNOTATION_COL})')
    parser.add_argument('--sample-col', default=DEFAULT_SAMPLE_COL,
                        help=f'Column with sample IDs (default: {DEFAULT_SAMPLE_COL})')
    parser.add_argument('--multi-sample', action='store_true',
                        help='Group by sample + cell type (for multi-sample experiments)')
    parser.add_argument('--method', choices=['deseq2', 'wilcoxon'], default='deseq2',
                        help='DE method (default: deseq2)')
    parser.add_argument('--min-cells', type=int, default=MIN_CELLS_PER_GROUP,
                        help=f'Minimum cells per pseudobulk group (default: {MIN_CELLS_PER_GROUP})')

    args = parser.parse_args()

    # Use args.min_cells directly instead of global

    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    os.makedirs(f"{args.output}/volcano_plots", exist_ok=True)

    print("=" * 60)
    print("PSEUDOBULK DIFFERENTIAL EXPRESSION")
    print("=" * 60)
    print(f"Input: {args.input}")
    print(f"Output: {args.output}")
    print(f"Annotation column: {args.annotation}")
    print(f"Method: {args.method}")
    print(f"Min cells per group: {args.min_cells}")

    # Load data
    print("\n" + "=" * 60)
    print("Loading data")
    print("=" * 60)

    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

    # Check annotation column exists
    if args.annotation not in adata.obs.columns:
        print(f"\nERROR: Annotation column '{args.annotation}' not found.")
        print(f"Available columns: {list(adata.obs.columns)}")
        sys.exit(1)

    # Show cell type distribution
    print(f"\nCell type distribution ({args.annotation}):")
    print(adata.obs[args.annotation].value_counts())

    # Determine grouping columns
    if args.multi_sample and args.sample_col in adata.obs.columns:
        groupby_cols = [args.sample_col, args.annotation]
        de_group_col = args.annotation  # DE still by cell type
    else:
        groupby_cols = [args.annotation]
        de_group_col = args.annotation

    # Create pseudobulk
    print("\n" + "=" * 60)
    print("Creating pseudobulk samples")
    print("=" * 60)

    counts_df, metadata_df = create_pseudobulk(
        adata,
        groupby_cols=groupby_cols,
        layer='counts',
        min_cells=args.min_cells
    )

    # Save pseudobulk data
    counts_df.to_csv(f"{args.output}/pseudobulk_counts.tsv", sep='\t')
    metadata_df.to_csv(f"{args.output}/pseudobulk_metadata.tsv", sep='\t')
    print(f"\nSaved: pseudobulk_counts.tsv, pseudobulk_metadata.tsv")

    # Get gene symbols for labeling
    gene_symbols = None
    if 'gene_symbol' in adata.var.columns:
        gene_symbols = adata.var['gene_symbol'].to_dict()

    # Run pairwise DE
    print("\n" + "=" * 60)
    print("Running differential expression")
    print("=" * 60)

    results_df = run_pairwise_de(
        counts_df,
        metadata_df,
        group_col=de_group_col,
        method=args.method
    )

    if len(results_df) == 0:
        print("\nNo DE results generated. Check your data and parameters.")
        sys.exit(1)

    # Add gene symbols to results
    if gene_symbols is not None:
        results_df['gene_symbol'] = results_df['gene'].map(gene_symbols)

    # Save results
    results_df.to_csv(f"{args.output}/de_results_all.tsv", sep='\t', index=False)
    print(f"\nSaved: de_results_all.tsv ({len(results_df)} total tests)")

    # Filter significant results
    sig_results = results_df[
        (results_df['padj'] < PADJ_THRESHOLD) &
        (results_df['log2FoldChange'].abs() > LOG2FC_THRESHOLD)
    ]
    sig_results.to_csv(f"{args.output}/de_results_significant.tsv", sep='\t', index=False)
    print(f"Saved: de_results_significant.tsv ({len(sig_results)} significant)")

    # Create visualizations
    print("\n" + "=" * 60)
    print("Creating visualizations")
    print("=" * 60)

    # Summary plot
    summary_df = plot_de_summary(results_df, args.output)
    summary_df.to_csv(f"{args.output}/de_summary.tsv", sep='\t')
    print("Saved: de_summary.png, de_summary.tsv")

    # Volcano plots for each comparison
    comparisons = results_df['comparison'].unique()
    print(f"Creating {len(comparisons)} volcano plots...")

    for comparison in comparisons:
        plot_volcano(results_df, comparison, f"{args.output}/volcano_plots", gene_symbols)

    print(f"Saved volcano plots to {args.output}/volcano_plots/")

    # Final summary
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nResults saved to: {args.output}/")
    print(f"  - pseudobulk_counts.tsv: Aggregated count matrix")
    print(f"  - pseudobulk_metadata.tsv: Sample metadata")
    print(f"  - de_results_all.tsv: All DE results")
    print(f"  - de_results_significant.tsv: Significant DE genes only")
    print(f"  - de_summary.tsv/png: Summary of DE gene counts")
    print(f"  - volcano_plots/: Volcano plot for each comparison")

    print(f"\nTotal comparisons: {len(comparisons)}")
    print(f"Total significant DE genes: {len(sig_results)}")


if __name__ == "__main__":
    main()
