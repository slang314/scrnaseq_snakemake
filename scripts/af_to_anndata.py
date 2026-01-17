#!/usr/bin/env python
"""
Convert alevin-fry output to AnnData format.

Alevin-fry outputs a sparse matrix in MTX format along with barcode and gene files.
This script reads those and creates a properly oriented AnnData object (cells x genes).
"""
import argparse
from pathlib import Path
import gzip
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad


def read_lines(path: Path):
    """Read lines from a file, handling gzip if needed."""
    if path.suffix == '.gz':
        with gzip.open(path, 'rt') as f:
            return [line.strip() for line in f]
    else:
        with open(path, 'r') as f:
            return [line.strip() for line in f]


def find_quant_dir(af_dir: Path):
    """Find the quant directory within alevin-fry output."""
    candidates = [
        af_dir / "quant",
        af_dir / "alevin",
        af_dir,
    ]
    for c in candidates:
        if (c / "alevin").exists():
            return c / "alevin"
        if any(c.glob("*.mtx*")):
            return c
    raise FileNotFoundError(f"Could not find quant output in {af_dir}")


def find_matrix_files(quant_dir: Path):
    """
    Locate matrix + barcodes + features in alevin-fry output.
    Returns (mtx_path, barcodes_path, features_path).
    """
    # Matrix file candidates
    mtx_candidates = [
        "quants_mat.mtx",
        "quants_mat.mtx.gz",
        "matrix.mtx",
        "matrix.mtx.gz",
    ]
    mtx_path = None
    for name in mtx_candidates:
        p = quant_dir / name
        if p.exists():
            mtx_path = p
            break

    if mtx_path is None:
        raise FileNotFoundError(f"Could not find matrix file in {quant_dir}")

    # Barcode file candidates (rows = cells)
    bc_candidates = [
        "quants_mat_rows.txt",
        "barcodes.tsv",
        "barcodes.tsv.gz",
    ]
    bc_path = None
    for name in bc_candidates:
        p = quant_dir / name
        if p.exists():
            bc_path = p
            break

    if bc_path is None:
        raise FileNotFoundError(f"Could not find barcodes file in {quant_dir}")

    # Features/genes file candidates (cols = genes)
    feat_candidates = [
        "quants_mat_cols.txt",
        "features.tsv",
        "features.tsv.gz",
        "genes.tsv",
        "genes.tsv.gz",
    ]
    feat_path = None
    for name in feat_candidates:
        p = quant_dir / name
        if p.exists():
            feat_path = p
            break

    if feat_path is None:
        raise FileNotFoundError(f"Could not find features file in {quant_dir}")

    return mtx_path, bc_path, feat_path


def main():
    ap = argparse.ArgumentParser(description="Convert alevin-fry output to AnnData")
    ap.add_argument("--sample", required=True, help="Sample name")
    ap.add_argument("--af-dir", required=True, help="Alevin-fry output directory")
    ap.add_argument("--out", required=True, help="Output h5ad file path")
    args = ap.parse_args()

    af_dir = Path(args.af_dir)

    # Find the directory containing the matrix
    quant_dir = find_quant_dir(af_dir)
    print(f"Found quant directory: {quant_dir}")

    # Find matrix files
    mtx_path, bc_path, feat_path = find_matrix_files(quant_dir)
    print(f"Matrix: {mtx_path}")
    print(f"Barcodes: {bc_path}")
    print(f"Features: {feat_path}")

    # Read matrix - alevin-fry outputs cells x genes (rows x cols)
    print("Reading matrix...")
    X = scipy.io.mmread(str(mtx_path))

    # Convert to CSR for efficient row slicing
    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = sp.csr_matrix(X)

    # Read barcodes (cells)
    barcodes = read_lines(bc_path)

    # Read features (genes)
    feat_lines = read_lines(feat_path)

    # Parse features - could be 1, 2, or 3 columns
    if '\t' in feat_lines[0]:
        feat_split = [line.split('\t') for line in feat_lines]
        n_cols = len(feat_split[0])

        if n_cols >= 2:
            gene_ids = [row[0] for row in feat_split]
            gene_names = [row[1] if len(row) > 1 else row[0] for row in feat_split]
            var = pd.DataFrame({
                'gene_id': gene_ids,
                'gene_name': gene_names,
            }, index=gene_ids)
        else:
            gene_ids = [row[0] for row in feat_split]
            var = pd.DataFrame(index=gene_ids)
    else:
        # Single column - just gene names/IDs
        var = pd.DataFrame(index=feat_lines)

    var.index.name = 'gene'

    # Verify dimensions
    print(f"Matrix shape: {X.shape}")
    print(f"Barcodes: {len(barcodes)}")
    print(f"Features: {len(var)}")

    # Check if we need to transpose
    # AnnData expects cells x genes
    if X.shape[0] == len(var) and X.shape[1] == len(barcodes):
        print("Transposing matrix (was genes x cells)")
        X = X.T.tocsr()
    elif X.shape[0] != len(barcodes) or X.shape[1] != len(var):
        raise ValueError(
            f"Matrix shape {X.shape} doesn't match barcodes ({len(barcodes)}) "
            f"and features ({len(var)})"
        )

    # Create obs (cell metadata)
    obs = pd.DataFrame(index=barcodes)
    obs.index.name = 'barcode'
    obs['sample'] = args.sample

    # Create AnnData
    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Add provenance info
    adata.uns['provenance'] = {
        'sample': args.sample,
        'source': 'alevin-fry',
        'af_dir': str(af_dir),
        'quant_dir': str(quant_dir),
    }

    # Write output
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print(f"Writing AnnData to {out_path}")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Genes: {adata.n_vars}")

    adata.write_h5ad(out_path)
    print("Done!")


if __name__ == "__main__":
    main()
