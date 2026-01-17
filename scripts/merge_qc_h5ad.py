#!/usr/bin/env python
"""
Merge multiple QC'd AnnData objects into a single combined object.
"""
import argparse
from pathlib import Path
import scanpy as sc


def main():
    ap = argparse.ArgumentParser(description="Merge QC'd h5ad files")
    ap.add_argument("--samples", required=True, help="Path to samples.tsv")
    ap.add_argument("--qc-h5ads", nargs="+", required=True, help="QC'd h5ad files to merge")
    ap.add_argument("--out", required=True, help="Output merged h5ad file")
    args = ap.parse_args()

    print(f"Merging {len(args.qc_h5ads)} samples...")

    # Read all AnnData objects
    adatas = []
    sample_names = []

    for path in args.qc_h5ads:
        sample_name = Path(path).stem
        print(f"  Reading {sample_name} from {path}")
        adata = sc.read_h5ad(path)

        if 'sample' not in adata.obs.columns:
            adata.obs['sample'] = sample_name

        adatas.append(adata)
        sample_names.append(sample_name)

    if len(adatas) == 0:
        raise ValueError("No AnnData objects to merge")

    if len(adatas) == 1:
        print("Only one sample, copying directly...")
        merged = adatas[0].copy()
    else:
        print("Concatenating samples...")
        merged = sc.concat(
            adatas,
            label='sample',
            keys=sample_names,
            index_unique='-',
            join='outer',
        )

    print(f"Merged AnnData: {merged.n_obs} cells x {merged.n_vars} genes")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(out_path)
    print("Done!")


if __name__ == "__main__":
    main()
