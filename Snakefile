"""
scRNA-seq Pipeline: alevin-fry quantification + Scanpy QC

This pipeline:
1. Quantifies 10x scRNA-seq data using salmon alevin + alevin-fry
2. Converts output to AnnData format
3. Runs Scanpy QC (earliest possible point after count matrix exists)
4. Optionally merges QC'd samples
"""

import pandas as pd
from pathlib import Path

configfile: "config/config.yaml"

# Load sample sheet
SAMPLES_DF = pd.read_csv("config/samples.tsv", sep="\t")
SAMPLES = SAMPLES_DF["sample"].tolist()


def get_fq1(wildcards):
    """Get R1 FASTQ path(s) for a sample as a list."""
    fq1 = SAMPLES_DF.set_index("sample").loc[wildcards.sample, "fq1"]
    return fq1.split(",")


def get_fq2(wildcards):
    """Get R2 FASTQ path(s) for a sample as a list."""
    fq2 = SAMPLES_DF.set_index("sample").loc[wildcards.sample, "fq2"]
    return fq2.split(",")


def get_fq1_str(wildcards):
    """Get R1 FASTQ path(s) as space-separated string for shell."""
    return " ".join(get_fq1(wildcards))


def get_fq2_str(wildcards):
    """Get R2 FASTQ path(s) as space-separated string for shell."""
    return " ".join(get_fq2(wildcards))


# Define final outputs
rule all:
    input:
        expand("results/anndata/qc/{sample}.h5ad", sample=SAMPLES),
        expand("results/qc/{sample}/qc_metrics.tsv", sample=SAMPLES),
        "results/anndata/merged_qc.h5ad" if config["downstream"].get("merge_qc", True) else []


# =============================================================================
# RULE 1: Run salmon alevin + alevin-fry quantification
# =============================================================================
rule alevin_fry_quant:
    input:
        r1=get_fq1,
        r2=get_fq2,
    output:
        quant_dir=directory("results/af/{sample}/quant"),
        done="results/af/{sample}/DONE.txt",
    params:
        af_ref=config["reference"]["af_dir"],
        salmon_flags=config["alevin_fry"].get("salmon_flags", "--chromiumV3"),
        outdir="results/af/{sample}",
        r1_str=get_fq1_str,
        r2_str=get_fq2_str,
    threads: config["alevin_fry"]["threads"]
    resources:
        mem_gb=config["alevin_fry"]["mem_gb"]
    conda:
        "envs/af.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        echo "R1 files: {params.r1_str}"
        echo "R2 files: {params.r2_str}"

        # Step 1: Run salmon alevin
        salmon alevin \
            -i {params.af_ref}/index \
            -l ISR \
            -1 {params.r1_str} \
            -2 {params.r2_str} \
            {params.salmon_flags} \
            -p {threads} \
            -o {params.outdir}/alevin \
            --sketch

        # Step 2: alevin-fry generate-permit-list (knee-distance for cell calling)
        alevin-fry generate-permit-list \
            -i {params.outdir}/alevin \
            -d fw \
            -o {params.outdir}/permit \
            --knee-distance

        # Step 3: alevin-fry collate
        alevin-fry collate \
            -i {params.outdir}/permit \
            -r {params.outdir}/alevin \
            -t {threads}

        # Step 4: alevin-fry quant
        alevin-fry quant \
            -i {params.outdir}/permit \
            -o {output.quant_dir} \
            -t {threads} \
            -m {params.af_ref}/t2g.tsv \
            -r cr-like \
            --use-mtx

        echo "done" > {output.done}
        """


# =============================================================================
# RULE 2: Convert alevin-fry output to AnnData
# =============================================================================
rule af_to_anndata:
    input:
        quant_dir="results/af/{sample}/quant",
        done="results/af/{sample}/DONE.txt",
    output:
        h5ad="results/anndata/raw/{sample}.h5ad",
    conda:
        "envs/scanpy.yaml"
    shell:
        r"""
        python scripts/af_to_anndata.py \
            --sample {wildcards.sample} \
            --af-dir results/af/{wildcards.sample} \
            --out {output.h5ad}
        """


# =============================================================================
# RULE 3: Scanpy QC (earliest possible point after count matrix)
# =============================================================================
rule scanpy_qc:
    input:
        h5ad="results/anndata/raw/{sample}.h5ad",
    output:
        qc_h5ad="results/anndata/qc/{sample}.h5ad",
        metrics="results/qc/{sample}/qc_metrics.tsv",
    params:
        mito_prefix=config["scanpy_qc"]["mito_prefix"],
        min_genes=config["scanpy_qc"]["min_genes"],
        min_cells_per_gene=config["scanpy_qc"]["min_cells_per_gene"],
        max_pct_mt=config["scanpy_qc"]["max_pct_mt"],
        max_counts=config["scanpy_qc"].get("max_counts"),
        max_genes=config["scanpy_qc"].get("max_genes"),
    conda:
        "envs/scanpy.yaml"
    shell:
        r"""
        python scripts/scanpy_qc.py \
            --input {input.h5ad} \
            --out {output.qc_h5ad} \
            --metrics {output.metrics} \
            --mito-prefix "{params.mito_prefix}" \
            --min-genes {params.min_genes} \
            --min-cells-per-gene {params.min_cells_per_gene} \
            --max-pct-mt {params.max_pct_mt} \
            $(if [ "{params.max_counts}" != "None" ] && [ -n "{params.max_counts}" ]; then echo "--max-counts {params.max_counts}"; fi) \
            $(if [ "{params.max_genes}" != "None" ] && [ -n "{params.max_genes}" ]; then echo "--max-genes {params.max_genes}"; fi)
        """


# =============================================================================
# RULE 4: Merge QC'd samples (optional)
# =============================================================================
rule merge_qc:
    input:
        h5ads=expand("results/anndata/qc/{sample}.h5ad", sample=SAMPLES),
    output:
        "results/anndata/merged_qc.h5ad",
    conda:
        "envs/scanpy.yaml"
    shell:
        r"""
        python scripts/merge_qc_h5ad.py \
            --samples config/samples.tsv \
            --qc-h5ads {input.h5ads} \
            --out {output}
        """
