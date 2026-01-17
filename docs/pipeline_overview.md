# scRNA-seq Pipeline Overview

## Complete Data Flow: FASTQs → VAE-UMAP

This document traces every step of the pipeline, noting where all information originates.

---

## 1. Input Data

### Raw FASTQ Files
```
data/fastq/pbmc_10k_v3_fastqs/
├── pbmc_10k_v3_S1_L001_R1_001.fastq.gz  # Read 1: Cell barcode (16bp) + UMI (12bp)
├── pbmc_10k_v3_S1_L001_R2_001.fastq.gz  # Read 2: cDNA sequence (transcript)
├── pbmc_10k_v3_S1_L002_R1_001.fastq.gz  # Lane 2
└── pbmc_10k_v3_S1_L002_R2_001.fastq.gz
```

**Source**: 10x Genomics public dataset "10k PBMCs from a Healthy Donor (v3 chemistry)"
- ~10,000 peripheral blood mononuclear cells
- 10x Chromium 3' v3 chemistry
- ~682 million read pairs

### Reference Files
```
resources/reference/
├── af_ref/
│   ├── index/              # Salmon index (built from transcriptome)
│   ├── t2g.tsv             # Transcript → Gene mapping
│   └── 3M-february-2018.txt # 10x barcode whitelist (3.7M valid barcodes)
├── annotation.gtf.gz       # GENCODE v44 (GRCh38) gene annotations
├── transcripts.fa.gz       # Human transcriptome sequences
└── gene_id_to_name.tsv     # Ensembl ID → Gene symbol mapping (extracted from GTF)
```

**Sources**:
- **Transcriptome/GTF**: GENCODE v44 (Ensembl 110), GRCh38 human genome
- **Barcode whitelist**: 10x Genomics (defines valid cell barcodes for v3 chemistry)
- **Salmon index**: Pre-built from transcriptome FASTA
- **Gene mapping**: Extracted from GENCODE GTF using custom script

---

## 2. Alignment & Quantification (Snakemake Rule: `alevin_fry_quant`)

### Step 2a: Salmon Alevin (Alignment)
```bash
salmon alevin \
    -i resources/reference/af_ref/index \    # Salmon index
    -l ISR \                                  # Library type: Inward Stranded Read
    -1 {R1_fastqs} \                          # Cell barcode + UMI reads
    -2 {R2_fastqs} \                          # Transcript reads
    --chromiumV3 \                            # 10x v3 chemistry (16bp BC + 12bp UMI)
    -p 16 \                                   # Threads
    -o results/af/pbmc_10k_v3/alevin \
    --sketch                                  # Sketch mode (faster, writes RAD file)
```

**What happens**:
1. Reads R2 (cDNA) are pseudoaligned to transcriptome
2. R1 is parsed for 16bp cell barcode + 12bp UMI
3. Outputs RAD (Reduced Alignment Data) file with mappings

**Output**: `results/af/pbmc_10k_v3/alevin/map.rad`

### Step 2b: Alevin-fry Permit List
```bash
alevin-fry generate-permit-list \
    -i results/af/pbmc_10k_v3/alevin \
    -d fw \                                   # Direction: forward
    -o results/af/pbmc_10k_v3/permit \
    --unfiltered-pl resources/reference/af_ref/3M-february-2018.txt
```

**What happens**:
1. Compares observed barcodes against 10x whitelist (3.7M valid barcodes)
2. Corrects barcodes with 1 edit distance from whitelist
3. Identifies "true" cells vs empty droplets

**Key info source**: `3M-february-2018.txt` - 10x Genomics barcode whitelist defines which 16bp sequences are valid cell barcodes

### Step 2c: Alevin-fry Collate
```bash
alevin-fry collate \
    -i results/af/pbmc_10k_v3/permit \
    -r results/af/pbmc_10k_v3/alevin \
    -t 16
```

**What happens**: Groups reads by corrected cell barcode for efficient counting

### Step 2d: Alevin-fry Quant
```bash
alevin-fry quant \
    -i results/af/pbmc_10k_v3/permit \
    -o results/af/pbmc_10k_v3/quant \
    -t 16 \
    -m resources/reference/af_ref/t2g.tsv \  # Transcript → Gene mapping
    -r cr-like \                              # Resolution: Cell Ranger-like
    --use-mtx                                 # Output as MatrixMarket format
```

**What happens**:
1. Collapses UMIs to get unique molecule counts
2. Maps transcripts → genes using `t2g.tsv`
3. Resolves multi-mapping reads (cr-like = Cell Ranger algorithm)

**Output**:
```
results/af/pbmc_10k_v3/quant/alevin/
├── quants_mat.mtx          # Sparse count matrix (cells × genes)
├── quants_mat_rows.txt     # Cell barcodes
└── quants_mat_cols.txt     # Gene IDs (Ensembl IDs with version)
```

**Key info source**: `t2g.tsv` maps transcript IDs → gene IDs (from GENCODE GTF)

---

## 3. Count Matrix to AnnData (Snakemake Rule: `af_to_anndata`)

```python
# scripts/af_to_anndata.py
adata = anndata.AnnData(X=count_matrix)
adata.obs_names = barcodes      # Cell barcodes
adata.var_names = gene_ids      # Ensembl gene IDs (e.g., ENSG00000141510.18)
```

**Output**: `results/anndata/raw/pbmc_10k_v3.h5ad`
- Raw UMI counts (integers)

---

## 4. Quality Control (Snakemake Rule: `scanpy_qc`)

```python
# scripts/scanpy_qc.py
# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

# Filter cells
sc.pp.filter_cells(adata, min_genes=10)      # Cells with ≥10 genes detected
sc.pp.filter_genes(adata, min_cells=3)       # Genes detected in ≥3 cells

# Filter by mitochondrial content
adata = adata[adata.obs['pct_counts_mt'] < 20]
```

**Filters applied** (from `config/config.yaml`):
| Parameter | Value | Effect |
|-----------|-------|--------|
| min_genes | 10 | Remove low-quality cells |
| min_cells_per_gene | 3 | Remove rare genes |
| max_pct_mt | 20% | Remove dying cells |

**Output**: `results/anndata/qc/pbmc_10k_v3.h5ad`

---

## 5. Downstream Analysis (scripts/downstream_analysis_v3.py)

### 5a. Gene ID Handling

```python
# Strip version from Ensembl IDs for stability
# ENSG00000141510.18 → ENSG00000141510
stripped_ids = [g.split('.')[0] for g in adata.var_names]
adata.var['gene_id_versioned'] = original_ids  # Preserve original
adata.var_names = stripped_ids                  # Use stable IDs

# Map to gene symbols (stored in var column, NOT as var_names)
gene_map = pd.read_csv("resources/reference/gene_id_to_name.tsv")
adata.var['gene_symbol'] = [gene_map.get(g, None) for g in adata.var_names]
```

**Key design decision**: Ensembl IDs remain as `var_names` (stable identifiers), gene symbols stored in `var['gene_symbol']` column. This avoids issues with:
- Duplicate gene symbols
- Missing mappings
- Downstream tool compatibility

**Source**: `resources/reference/gene_id_to_name.tsv` extracted from GENCODE v44 GTF

### 5b. Preserve Raw Counts

```python
adata.layers['counts'] = adata.X.copy()  # Raw UMI counts
adata.obs['size_factors'] = adata.obs['total_counts'].copy()
```

### 5c. Normalization

```python
# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata, target_sum=1e4)

# Log transform: log(x + 1)
sc.pp.log1p(adata)

# Store normalized data for DE analysis (all genes)
adata.raw = adata
```

### 5d. HVG Selection

```python
sc.pp.highly_variable_genes(
    adata,
    min_mean=0.0125,
    max_mean=3,
    min_disp=0.5,
    subset=False  # Mark but don't subset
)
# Result: ~3,029 HVGs
```

### 5e. Scaling and PCA (Consistent Path)

```python
# CRITICAL: Subset to HVGs first, then scale, then PCA on SAME object
adata_hvg = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50, random_state=42)

# Copy results back to main object
adata.obsm['X_pca'] = adata_hvg.obsm['X_pca']
adata.uns['pca'] = adata_hvg.uns['pca']
```

**Why this matters**: Scaling and PCA must operate on the same HVG subset. Previous versions scaled all genes but ran PCA on HVGs only, creating inconsistency.

### 5f. Neighbors and UMAP

```python
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, random_state=42)
sc.tl.umap(adata, random_state=42)
```

### 5g. Leiden Clustering

```python
for res in [0.2, 0.4, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}', random_state=42)

adata.obs['leiden'] = adata.obs['leiden_0.4']  # Default
```

---

## 6. Doublet Detection

```python
import scrublet as scr

scrub = scr.Scrublet(adata.layers['counts'], expected_doublet_rate=0.06, random_state=42)
doublet_scores, _ = scrub.scrub_doublets()

# Manual threshold (auto-threshold is often unreliable)
DOUBLET_THRESHOLD = 0.25
adata.obs['doublet_score'] = doublet_scores
adata.obs['predicted_doublet'] = doublet_scores > DOUBLET_THRESHOLD
```

**How Scrublet works**:
1. Simulates doublets by averaging random cell pairs
2. Embeds real + simulated cells in PCA space
3. Scores cells by k-NN density of simulated doublets nearby

**Why manual threshold**: Scrublet's auto-threshold can be unreliable (e.g., 0.016 detecting 75% doublets). Manual threshold of 0.25 is more reasonable for 10x data.

**Expected doublet rate**: ~6% for 10k cells (10x formula: ~0.8% per 1000 cells)

---

## 7. Cell Type Annotation (Post-hoc)

### Celltypist (Automated)

```python
import celltypist
from celltypist import models

# Create temporary copy with gene symbols for celltypist
adata_ct = adata.raw.to_adata().copy()
adata_ct.var_names = [adata.var.loc[g, 'gene_symbol'] or g for g in adata_ct.var_names]
adata_ct.var_names_make_unique()

# Download and load pre-trained model
models.download_models(model='Immune_All_Low.pkl')
model = models.Model.load(model='Immune_All_Low.pkl')

# Predict cell types
predictions = celltypist.annotate(adata_ct, model=model, majority_voting=True)
adata.obs['celltypist_label'] = predictions.to_adata().obs['majority_voting']
```

**Model source**: `Immune_All_Low.pkl`
- From Celltypist model repository (https://www.celltypist.org/)
- Trained on: ~20 immune cell datasets, millions of cells
- Reference: Domínguez Conde et al., Science 2022
- Contains: Logistic regression weights for ~10,000 genes
- Outputs: 100+ immune cell subtypes

**How it works**:
1. Matches input genes to model's gene set (~4,692 genes overlap)
2. Applies pre-trained logistic regression
3. Over-clusters data, then majority votes per cluster

**Important**: Celltypist labels are used for **post-hoc annotation only**, not for training any models.

---

## 8. VAE-based Embedding (scripts/vae_embedding.py)

### Architecture

```python
class VAE(nn.Module):
    # Encoder: gene expression → latent mean & variance
    # 3,029 HVGs → [256, 128] → 20D latent (μ, σ)

    # Decoder: latent → reconstructed expression
    # 20D → [128, 256] → 3,029 genes

    # Loss: Negative binomial reconstruction + KL divergence
```

### Training Details

```python
# UNSUPERVISED - no cell type labels used
input_dim = 3029      # HVGs
hidden_dims = [256, 128]
latent_dim = 20
epochs = 100
batch_size = 256
learning_rate = 1e-3
random_seed = 42
```

**Key design decision**: The VAE is **purely unsupervised**. Using celltypist labels during VAE training would create circular reasoning (labels derived from the data used to train embedding that produces the data).

### VAE UMAP and Clustering

```python
# Store original PCA-UMAP
adata.obsm['X_umap_pca'] = adata.obsm['X_umap'].copy()

# Compute new embedding from VAE latent space
sc.pp.neighbors(adata, use_rep='X_vae', n_neighbors=15, random_state=42)
sc.tl.umap(adata, random_state=42)
adata.obsm['X_umap_vae'] = adata.obsm['X_umap'].copy()

# Cluster in VAE space
sc.tl.leiden(adata, resolution=0.4, key_added='leiden_vae', random_state=42)
```

---

## 9. Marker Gene Detection

```python
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', key_added='markers_leiden')

# Get results with gene symbols
markers_df = sc.get.rank_genes_groups_df(adata, group=None, key='markers_leiden')
markers_df['gene_symbol'] = markers_df['names'].map(lambda x: adata.var.loc[x, 'gene_symbol'])
```

### Canonical PBMC Markers (from literature)

```python
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
```

**Source**: Well-established in immunology literature and 10x Genomics tutorials.

---

## 10. Final Outputs

### Output Directory Structure

```
results/downstream_v3/
├── pbmc_10k_annotated_v3.h5ad           # Full dataset with all annotations
├── pbmc_10k_annotated_v3_clean.h5ad     # Doublets removed
├── pbmc_10k_annotated_v3_vae.h5ad       # With VAE embeddings
├── vae_model.pt                          # Saved VAE model weights
├── marker_genes_all_clusters.tsv        # DE genes per cluster
├── analysis_summary.tsv                  # Summary statistics
└── figures/
    ├── umap_overview.png
    ├── umap_pbmc_markers.png
    ├── umap_celltypist.png
    ├── umap_pca_vs_vae.png
    ├── umap_vae_markers.png
    ├── doublet_scores.png
    └── vae_training_loss.png
```

### AnnData Structure

```
AnnData object with n_obs × n_vars = [cells] × [genes]
    obs: 'n_genes_by_counts', 'total_counts', 'pct_counts_mt',
         'size_factors', 'leiden', 'leiden_0.2', ..., 'leiden_1.0',
         'leiden_vae', 'doublet_score', 'predicted_doublet',
         'celltypist_label', 'celltypist_conf'
    var: 'mt', 'n_cells_by_counts', 'gene_id_versioned', 'gene_symbol',
         'highly_variable', 'means', 'dispersions', 'dispersions_norm'
    uns: 'pca', 'neighbors', 'umap', 'leiden', 'markers_leiden'
    obsm: 'X_pca' (50D), 'X_umap', 'X_umap_pca', 'X_vae' (20D), 'X_umap_vae'
    layers: 'counts' (raw UMI counts)
    raw: normalized log-expression for all genes (for DE)
```

### Data States

| Location | Contents |
|----------|----------|
| `adata.layers['counts']` | Raw UMI counts (integers) |
| `adata.X` | Log-normalized expression |
| `adata.raw.X` | Log-normalized, all genes (for DE) |
| `adata.obsm['X_pca']` | 50 principal components |
| `adata.obsm['X_vae']` | 20D VAE latent space |
| `adata.obsm['X_umap_pca']` | UMAP from PCA |
| `adata.obsm['X_umap_vae']` | UMAP from VAE |

---

## Information Flow Summary

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           EXTERNAL RESOURCES                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│  GENCODE v44 GTF ──────────► Gene IDs, Gene Symbols, Transcript mappings   │
│  10x Barcode Whitelist ────► Valid cell barcodes (3.7M sequences)          │
│  Human Transcriptome ──────► Salmon index for alignment                     │
│  Celltypist Model ─────────► Pre-trained immune cell classifier            │
│  Literature ───────────────► Canonical marker genes (CD3D, CD14, etc.)     │
└─────────────────────────────────────────────────────────────────────────────┘
                                      │
                                      ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                              PIPELINE FLOW                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  FASTQ ──► Salmon Alevin ──► RAD file ──► Alevin-fry ──► Count Matrix      │
│                 │                              │                             │
│           (uses index)                   (uses t2g.tsv)                      │
│                                          (uses whitelist)                    │
│                                                │                             │
│                                                ▼                             │
│                                         QC Filtering                         │
│                                    (min_genes, min_cells)                    │
│                                                │                             │
│                                                ▼                             │
│                                         Gene ID Fix                          │
│                                    (strip version, map symbols)              │
│                                                │                             │
│                                                ▼                             │
│                                    Normalize + Log1p                         │
│                                    (preserve counts in layer)                │
│                                                │                             │
│                                                ▼                             │
│                                         HVG Selection                        │
│                                        (~3,029 genes)                        │
│                                                │                             │
│                               ┌────────────────┴────────────────┐            │
│                               │                                 │            │
│                               ▼                                 ▼            │
│                         Scale + PCA                       VAE Training       │
│                          (50 PCs)                    (UNSUPERVISED, 20D)     │
│                               │                                 │            │
│                               ▼                                 ▼            │
│                            k-NN ◄──────────────────────────► k-NN           │
│                               │                                 │            │
│                               ▼                                 ▼            │
│                            UMAP                              UMAP            │
│                               │                                 │            │
│                               ▼                                 ▼            │
│                      Leiden Clustering              Leiden Clustering        │
│                      (res=0.4 → ~97)                (res=0.4 → ~36)         │
│                               │                                 │            │
│                               └────────────┬────────────────────┘            │
│                                            │                                 │
│                                            ▼                                 │
│                                    Celltypist Annotation                     │
│                                  (POST-HOC labels only)                      │
│                                            │                                 │
│                                            ▼                                 │
│                                   Doublet Detection                          │
│                                  (Scrublet, manual threshold)                │
│                                            │                                 │
│                                            ▼                                 │
│                                    Marker Detection                          │
│                                      (Wilcoxon)                              │
│                                            │                                 │
│                                            ▼                                 │
│                                    Save .h5ad files                          │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Reproducibility

### Software Versions (logged at runtime)
```
Python: 3.11.x
Scanpy: 1.11.5
NumPy: 2.3.5
Pandas: 2.3.3
AnnData: 0.12.7
SciPy: 1.17.0
PyTorch: (with CUDA if available)
Pyro: (latest)
Celltypist: (latest)
```

### Random Seeds

All stochastic operations use `RANDOM_SEED = 42`:
- NumPy: `np.random.seed(42)`
- Python: `random.seed(42)`
- PyTorch: `torch.manual_seed(42)`
- Scanpy PCA: `random_state=42`
- Scanpy UMAP: `random_state=42`
- Scanpy Leiden: `random_state=42`
- Scrublet: `random_state=42`

### Key Parameters

```yaml
# config/config.yaml
alevin_fry:
  chemistry: "3prime-v3"
  threads: 16

scanpy_qc:
  min_genes: 10
  min_cells_per_gene: 3
  max_pct_mt: 20

# HVG selection
hvg:
  min_mean: 0.0125
  max_mean: 3
  min_disp: 0.5

# Dimensionality reduction
pca:
  n_comps: 50

neighbors:
  n_neighbors: 15
  n_pcs: 30

# Clustering
leiden:
  resolution: 0.4  # default

# VAE (in vae_embedding.py)
vae:
  hidden_dims: [256, 128]
  latent_dim: 20
  epochs: 100
  batch_size: 256
  learning_rate: 1e-3

# Doublet detection
scrublet:
  expected_doublet_rate: 0.06
  threshold: 0.25  # manual, not auto
```

---

## Scripts

| Script | Purpose |
|--------|---------|
| `scripts/af_to_anndata.py` | Convert alevin-fry output to AnnData |
| `scripts/scanpy_qc.py` | Quality control filtering |
| `scripts/downstream_analysis_v3.py` | Main analysis (normalize, HVG, PCA, clustering, annotation) |
| `scripts/vae_embedding.py` | Unsupervised VAE training and embedding |
