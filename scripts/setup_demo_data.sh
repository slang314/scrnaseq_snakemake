#!/bin/bash
set -euo pipefail

# ============================================================
# Download demo FASTQ data (10x PBMC datasets)
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$PROJECT_DIR/data/fastq"

# 10x Genomics PBMC datasets
PBMC_1K_URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar"
PBMC_10K_URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_fastqs.tar"

echo "============================================================"
echo "Demo data setup for scRNA-seq pipeline"
echo "============================================================"
echo ""
echo "Available datasets:"
echo "  1) pbmc_1k_v3  - 1,000 PBMCs (~5 GB, faster for testing)"
echo "  2) pbmc_10k_v3 - 10,000 PBMCs (~48 GB, full demo)"
echo "  3) both"
echo ""

read -p "Which dataset to download? [1/2/3]: " choice

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

download_and_extract() {
    local name=$1
    local url=$2

    if [ -d "${name}_fastqs" ]; then
        echo "$name already exists, skipping."
        return
    fi

    echo ""
    echo "Downloading $name..."
    wget -q --show-progress -O "${name}_fastqs.tar" "$url"

    echo "Extracting..."
    tar -xf "${name}_fastqs.tar"

    echo "$name ready."
}

case $choice in
    1)
        download_and_extract "pbmc_1k_v3" "$PBMC_1K_URL"
        ;;
    2)
        download_and_extract "pbmc_10k_v3" "$PBMC_10K_URL"
        ;;
    3)
        download_and_extract "pbmc_1k_v3" "$PBMC_1K_URL"
        download_and_extract "pbmc_10k_v3" "$PBMC_10K_URL"
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo ""
echo "============================================================"
echo "Demo data setup complete!"
echo "============================================================"
echo ""
echo "FASTQ files are in: $DATA_DIR"
echo ""
echo "Update config/samples.tsv to point to your data, e.g.:"
echo ""
echo "sample        fq1                                          fq2"
echo "pbmc_1k_v3    $DATA_DIR/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R1_001.fastq.gz,...    $DATA_DIR/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_R2_001.fastq.gz,..."
echo ""
