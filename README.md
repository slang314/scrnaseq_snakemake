# scRNA-seq Snakemake Pipeline

A reproducible single-cell RNA-seq analysis pipeline using alevin-fry for quantification and Scanpy for QC and downstream analysis.

## Overview

This pipeline processes 10x Chromium scRNA-seq data through:

1. **Quantification** - salmon alevin + alevin-fry (pseudoalignment, UMI deduplication, cell calling)
2. **QC** - Scanpy quality control with filtering and visualization
3. **Downstream analysis** - Clustering, cell type annotation, differential expression, marker genes, VAE embedding

## Quick Start

```bash
# 1. Clone the repository
git clone https://github.com/slang314/scrnaseq_snakemake.git
cd scrnaseq_snakemake

# 2. Create conda environments
conda env create -f envs/af.yaml
conda env create -f envs/scanpy.yaml

# 3. Build reference index (requires af environment)
conda activate af
./scripts/setup_reference.sh

# 4. Download demo data
./scripts/setup_demo_data.sh

# 5. Run the pipeline
conda deactivate
snakemake --cores 16 --use-conda
```

## Requirements

- [Conda](https://docs.conda.io/en/latest/miniconda.html) or [Mamba](https://mamba.readthedocs.io/)
- ~50 GB disk space for reference + demo data
- 16+ GB RAM recommended

## Project Structure

```
scrnaseq_snakemake/
├── Snakefile                 # Main pipeline definition
├── config/
│   ├── config.yaml           # Pipeline parameters
│   └── samples.tsv           # Sample sheet (FASTQ paths)
├── envs/
│   ├── af.yaml               # salmon/alevin-fry environment
│   └── scanpy.yaml           # Python analysis environment
├── scripts/
│   ├── setup_reference.sh    # Download & build reference
│   ├── setup_demo_data.sh    # Download demo FASTQ data
│   ├── af_to_anndata.py      # Convert alevin-fry → AnnData
│   ├── scanpy_qc.py          # QC filtering
│   ├── downstream_analysis_v3.py  # Clustering & annotation
│   └── vae_embedding.py      # VAE latent space analysis
├── resources/
│   └── reference/af_ref/     # Reference files (built by setup script)
├── data/                     # FASTQ files (user-provided)
└── results/                  # Pipeline outputs
```

## Setup

### 1. Reference Files

The pipeline requires a salmon index for the human transcriptome. Run the setup script:

```bash
conda activate af
./scripts/setup_reference.sh
```

This downloads from GENCODE and builds:
- `resources/reference/af_ref/index/` - Salmon index (~900 MB)
- `resources/reference/af_ref/t2g.tsv` - Transcript-to-gene mapping
- `resources/reference/af_ref/3M-february-2018.txt` - 10x barcode whitelist

### 2. Demo Data

Download 10x Genomics PBMC datasets for testing:

```bash
./scripts/setup_demo_data.sh
```

Options:
- **pbmc_1k_v3** - 1,000 cells (~5 GB) - quick testing
- **pbmc_10k_v3** - 10,000 cells (~48 GB) - full demo

### 3. Sample Sheet

Edit `config/samples.tsv` to point to your FASTQ files:

```
sample          fq1                                     fq2
pbmc_10k_v3     /path/to/sample_L001_R1_001.fastq.gz,/path/to/sample_L002_R1_001.fastq.gz    /path/to/sample_L001_R2_001.fastq.gz,/path/to/sample_L002_R2_001.fastq.gz
```

- Multiple lanes: comma-separate paths (no spaces)
- R1 = barcode + UMI reads
- R2 = cDNA reads

## Running the Pipeline

### Main Pipeline (Snakemake)

```bash
# Dry run - see what will be executed
snakemake --dry-run --cores 16

# Run pipeline
snakemake --cores 16 --use-conda

# Force re-run all steps
snakemake --cores 16 --use-conda --forceall
```

### Downstream Analysis (Manual)

After the main pipeline completes:

```bash
conda activate scanpy

# Clustering, cell typing, marker genes
python scripts/downstream_analysis_v3.py

# VAE embedding (optional)
python scripts/vae_embedding.py
```

## Outputs

### Main Pipeline

```
results/
├── af/{sample}/quant/           # alevin-fry count matrices
├── anndata/
│   ├── raw/{sample}.h5ad        # Raw counts
│   ├── qc/{sample}.h5ad         # QC-filtered
│   └── merged_qc.h5ad           # Merged samples
└── qc/{sample}/
    ├── qc_metrics.tsv           # QC statistics
    ├── {sample}_qc_violin.png   # QC violin plots
    └── {sample}_qc_scatter.png  # QC scatter plots
```

### Downstream Analysis

```
results/downstream_v3/
├── figures/
│   ├── barcode_rank_plot.png    # Knee plot for cell calling
│   ├── umap_overview.png        # UMAP with clusters
│   ├── umap_celltypist.png      # Cell type annotations
│   ├── umap_pbmc_markers.png    # Marker gene expression
│   ├── doublet_scores.png       # Scrublet doublet detection
│   ├── volcano_plots_all_celltypes.png  # DE volcano plots (overview)
│   ├── volcano_individual/      # Individual volcano plots per cell type
│   ├── vae_training_loss.png    # VAE training curve
│   ├── umap_pca_vs_vae.png      # PCA vs VAE comparison
│   └── umap_vae_markers.png     # Markers on VAE UMAP
├── marker_genes_all_clusters.tsv    # Marker genes per Leiden cluster
├── de_genes_by_celltype.tsv         # DE genes per cell type (CellTypist)
├── pbmc_10k_annotated_v3.h5ad       # Full annotated dataset
└── pbmc_10k_annotated_v3_clean.h5ad # Doublets removed
```

## Configuration

Edit `config/config.yaml`:

```yaml
reference:
  af_dir: "resources/reference/af_ref"

alevin_fry:
  threads: 16
  mem_gb: 48
  chemistry: "3prime-v3"

scanpy_qc:
  species: "human"
  mito_prefix: "MT-"
  min_genes: 200
  min_cells_per_gene: 3
  max_pct_mt: 20
```

## Conda Environments

| Environment | Purpose | Create |
|-------------|---------|--------|
| `af` | salmon, alevin-fry | `conda env create -f envs/af.yaml` |
| `scanpy` | Scanpy, analysis | `conda env create -f envs/scanpy.yaml` |

Snakemake automatically uses the correct environment for each rule when run with `--use-conda`.

## Adding Your Own Data

1. Place FASTQ files in `data/fastq/your_sample/`
2. Add entry to `config/samples.tsv`
3. Run: `snakemake --cores 16 --use-conda`

## License

MIT

## Citation

If you use this pipeline, please cite:
- [salmon](https://doi.org/10.1038/nmeth.4197)
- [alevin-fry](https://doi.org/10.1038/s41592-022-01408-3)
- [Scanpy](https://doi.org/10.1186/s13059-017-1382-0)
