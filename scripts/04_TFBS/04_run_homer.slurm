#!/bin/bash

#SBATCH --job-name=EDTA_TFBS
#SBATCH --partition=64c512g
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=logs/EDTA_TFBS_%j.out
#SBATCH --error=logs/EDTA_TFBS_%j.err

# --------------------------------------------------
# Step 04: HOMER motif enrichment analysis
#
# Target: ATAC-seq peaks overlapping target gene bodies
# Background: All ATAC-seq peaks
# Reference genome: danRer11 (GRCz11)
# --------------------------------------------------

set -eo pipefail


# --------------------------------------------------
# 1. Set project directory
# --------------------------------------------------

PROJECT_DIR="$HOME/projects/EDTA"

cd "$PROJECT_DIR"


# --------------------------------------------------
# 2. Activate HOMER environment
# --------------------------------------------------

source /dssg/home/acct-esezk/esezk-wangpeng6/.conda/etc/profile.d/conda.sh 2>/dev/null || true

module load miniconda3

source activate /dssg/home/acct-esezk/esezk-wangpeng6/.conda/envs/meme/


# --------------------------------------------------
# 3. Define analysis parameters
# --------------------------------------------------

GENOME="danRer11"

COMPOUND="${COMPOUND:?Please specify COMPOUND=BET or FLU}"

if [[ "$COMPOUND" != "BET" && "$COMPOUND" != "FLU" ]]; then
    echo "Error: COMPOUND must be BET or FLU."
    exit 1
fi

THREADS="${SLURM_CPUS_PER_TASK:-8}"


# --------------------------------------------------
# 4. Define input and output paths
# --------------------------------------------------

INPUT_DIR="$PROJECT_DIR/results/04_TFBS"

TARGET="$INPUT_DIR/${COMPOUND}_posi.homer.pos"

BACKGROUND="$INPUT_DIR/ATAC_background.homer.pos"

OUTPUT_DIR="$INPUT_DIR/${COMPOUND}_posi_motif"


# --------------------------------------------------
# 5. Check input files
# --------------------------------------------------

if [[ ! -s "$TARGET" ]]; then
    echo "Error: Target file not found or empty: $TARGET"
    exit 1
fi

if [[ ! -s "$BACKGROUND" ]]; then
    echo "Error: Background file not found or empty: $BACKGROUND"
    exit 1
fi

if ! command -v findMotifsGenome.pl >/dev/null 2>&1; then
    echo "Error: HOMER is not available in the current environment."
    exit 1
fi

mkdir -p "$OUTPUT_DIR"


# --------------------------------------------------
# 6. Display analysis information
# --------------------------------------------------

echo "=========================================="
echo "EDTA: HOMER Motif Enrichment"
echo "=========================================="
echo "Compound: $COMPOUND"
echo "Genome: $GENOME"
echo "Target: $TARGET"
echo "Background: $BACKGROUND"
echo "Motif sequence size: 200 bp"
echo "Threads: $THREADS"
echo "Output: $OUTPUT_DIR"
echo "=========================================="


# --------------------------------------------------
# 7. Run HOMER motif enrichment
# --------------------------------------------------

findMotifsGenome.pl \
    "$TARGET" \
    "$GENOME" \
    "$OUTPUT_DIR" \
    -size given \
    -p "$THREADS" \
    -mask


# --------------------------------------------------
# 8. Finish
# --------------------------------------------------

echo "=========================================="
echo "HOMER analysis completed successfully."
echo "Compound: $COMPOUND"
echo "=========================================="
