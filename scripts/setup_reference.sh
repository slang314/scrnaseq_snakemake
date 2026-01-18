#!/bin/bash
set -euo pipefail

# ============================================================
# Setup script for scRNA-seq pipeline reference files
# Downloads and builds salmon index for human transcriptome
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
REF_DIR="$PROJECT_DIR/resources/reference"
AF_REF="$REF_DIR/af_ref"

# GENCODE release (human)
GENCODE_RELEASE="44"
GENCODE_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"

# 10x barcode whitelist
BARCODE_URL="https://github.com/10XGenomics/cellranger/raw/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"

echo "============================================================"
echo "Setting up reference files for scRNA-seq pipeline"
echo "============================================================"
echo "Project directory: $PROJECT_DIR"
echo "Reference directory: $REF_DIR"
echo ""

# Create directories
mkdir -p "$AF_REF/index"
mkdir -p "$REF_DIR/download"

cd "$REF_DIR/download"

# ------------------------------------------------------------
# 1. Download GENCODE transcriptome and GTF
# ------------------------------------------------------------
echo "1. Downloading GENCODE v${GENCODE_RELEASE} files..."

FASTA="gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"
GTF="gencode.v${GENCODE_RELEASE}.annotation.gtf.gz"

if [ ! -f "$FASTA" ]; then
    echo "   Downloading transcriptome FASTA..."
    wget -q --show-progress "${GENCODE_BASE}/${FASTA}"
else
    echo "   FASTA already exists, skipping."
fi

if [ ! -f "$GTF" ]; then
    echo "   Downloading GTF annotation..."
    wget -q --show-progress "${GENCODE_BASE}/${GTF}"
else
    echo "   GTF already exists, skipping."
fi

# ------------------------------------------------------------
# 2. Download 10x barcode whitelist
# ------------------------------------------------------------
echo ""
echo "2. Downloading 10x Chromium v3 barcode whitelist..."

if [ ! -f "$AF_REF/3M-february-2018.txt" ]; then
    wget -q --show-progress -O 3M-february-2018.txt.gz "$BARCODE_URL"
    gunzip -c 3M-february-2018.txt.gz > "$AF_REF/3M-february-2018.txt"
    rm 3M-february-2018.txt.gz
    echo "   Saved to $AF_REF/3M-february-2018.txt"
else
    echo "   Barcode whitelist already exists, skipping."
fi

# ------------------------------------------------------------
# 3. Build transcript-to-gene mapping
# ------------------------------------------------------------
echo ""
echo "3. Building transcript-to-gene mapping..."

if [ ! -f "$AF_REF/t2g.tsv" ]; then
    python "$SCRIPT_DIR/build_t2g.py" "$GTF" "$AF_REF/t2g.tsv"
else
    echo "   t2g.tsv already exists, skipping."
fi

# ------------------------------------------------------------
# 4. Build gene ID to gene name mapping
# ------------------------------------------------------------
echo ""
echo "4. Building gene ID to name mapping..."

GENE_MAP="$REF_DIR/gene_id_to_name.tsv"
if [ ! -f "$GENE_MAP" ]; then
    echo "   Extracting gene IDs and names from GTF..."
    zcat "$GTF" | awk -F'\t' '
        $3 == "gene" {
            gene_id = ""; gene_name = ""
            n = split($9, attrs, ";")
            for (i = 1; i <= n; i++) {
                gsub(/^ +| +$/, "", attrs[i])
                if (attrs[i] ~ /^gene_id/) {
                    gsub(/gene_id "/, "", attrs[i])
                    gsub(/"/, "", attrs[i])
                    gene_id = attrs[i]
                }
                if (attrs[i] ~ /^gene_name/) {
                    gsub(/gene_name "/, "", attrs[i])
                    gsub(/"/, "", attrs[i])
                    gene_name = attrs[i]
                }
            }
            if (gene_id != "" && gene_name != "") {
                print gene_id "\t" gene_name
            }
        }
    ' > "$GENE_MAP"
    echo "   Created $GENE_MAP ($(wc -l < "$GENE_MAP") genes)"
else
    echo "   gene_id_to_name.tsv already exists, skipping."
fi

# ------------------------------------------------------------
# 5. Build salmon index
# ------------------------------------------------------------
echo ""
echo "5. Building salmon index (this may take 10-20 minutes)..."

if [ ! -f "$AF_REF/index/info.json" ]; then
    # Check if salmon is available
    if ! command -v salmon &> /dev/null; then
        echo "ERROR: salmon not found. Please activate the af conda environment:"
        echo "  conda env create -f envs/af.yaml"
        echo "  conda activate af"
        exit 1
    fi

    echo "   Decompressing FASTA..."
    gunzip -c "$FASTA" > transcripts.fa

    echo "   Running salmon index..."
    salmon index \
        -t transcripts.fa \
        -i "$AF_REF/index" \
        -k 31 \
        --gencode \
        -p 8

    rm transcripts.fa
    echo "   Index built successfully."
else
    echo "   Salmon index already exists, skipping."
fi

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
echo ""
echo "============================================================"
echo "Reference setup complete!"
echo "============================================================"
echo ""
echo "Files created:"
echo "  $AF_REF/index/           - Salmon index"
echo "  $AF_REF/t2g.tsv          - Transcript-to-gene mapping"
echo "  $AF_REF/3M-february-2018.txt - 10x barcode whitelist"
echo "  $REF_DIR/gene_id_to_name.tsv - Gene ID to symbol mapping"
echo ""
echo "Downloaded files cached in:"
echo "  $REF_DIR/download/"
echo ""
