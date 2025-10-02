#!/bin/bash
# Run ALL Pipeline Steps Automatically
# =====================================
# This script runs the entire pipeline with the --all flag
# No prompts, runs everything sequentially
#
# Usage: bash run_all_steps.sh [conda_env_name]
#   conda_env_name: optional, name of conda environment (default: auto-detect)

set -e

# Get script directory (works regardless of where script is called from)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Working directory: $SCRIPT_DIR"

# Activate conda environment (try multiple approaches)
CONDA_ENV="${1:-BILL}"  # Use first argument or default to BILL

if command -v conda &> /dev/null; then
    echo "Conda found. Attempting to activate environment: $CONDA_ENV"
    
    # Try multiple conda initialization paths
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
    elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
        source "/opt/conda/etc/profile.d/conda.sh"
    fi
    
    # Try to activate the environment
    if conda activate "$CONDA_ENV" 2>/dev/null; then
        echo "✓ Conda environment '$CONDA_ENV' activated"
    else
        echo "⚠ Warning: Could not activate conda environment '$CONDA_ENV'"
        echo "  Continuing with system Python..."
    fi
else
    echo "⚠ Conda not found. Using system Python."
fi

# Verify Python is available
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found in PATH"
    exit 1
fi

echo "Using Python: $(which python3)"
echo ""

echo "========================================="
echo "RUNNING COMPLETE PIPELINE WITH --all"
echo "========================================="
echo ""
echo "This will run ALL pipeline steps:"
echo "  1. Pre-QC"
echo "  2. Kraken (if available)"
echo "  3. Post-QC"
echo "  4. QUAST Pre"
echo "  5. BUSCO Pre"
echo "  6. QUAST Post"
echo "  7. BUSCO Post"
echo "  8. QC Merge"
echo "  9. Quality Flags"
echo " 10. Variants"
echo " 11. Pangenome-Cactus"
echo " 12. Pangenome-PGGB"
echo " 13. Phylogeny"
echo " 14. Synteny"
echo " 15. Functional"
echo " 16. Summary"
echo ""
echo "WARNING: This may take several hours!"
echo ""

# Give user 5 seconds to cancel
echo "Starting in 5 seconds... (Ctrl+C to cancel)"
sleep 5

echo ""
echo "Starting pipeline execution..."
echo "Start time: $(date)"
echo ""

# Configuration - EDIT THESE PATHS FOR YOUR SYSTEM
# ==================================================

# Input genomes (cleaned/assembled)
CLEANED_GENOMES=(
    "p15_1:clean_fasta/p15_1_clean.fasta"
    "p15_6:clean_fasta/p15_6_clean.fasta"
    "p90_2:clean_fasta/p90_2_clean.fasta"
    "p90_6:clean_fasta/p90_6_clean.fasta"
)

# Reference genome for pangenome
REFERENCE_FASTA="clean_fasta/ref/KHV-U_trunc.fasta"

# Optional: Raw genomes for pre-QC (comment out if not needed)
# RAW_GENOMES=(
#     "p15:data/p15_khv.fasta"
#     "p90:data/p90_khv.fasta"
# )

# Number of threads (auto-detect or set manually)
THREADS=${SLURM_CPUS_PER_TASK:-$(nproc 2>/dev/null || echo 4)}
echo "Using $THREADS threads"

# BUSCO lineage (adjust for your organism)
BUSCO_LINEAGE="viruses_odb10"

# Expected genome characteristics (adjust for your organism)
EXPECTED_SIZE=295146
EXPECTED_GC=59.2

# PGGB parameters
PGGB_SEGMENT=3000
PGGB_IDENTITY=95
PGGB_PASSES=10

# ==================================================
# End of configuration
# ==================================================

echo ""
echo "Checking input files..."
echo ""

# Check if reference exists
if [ ! -f "$REFERENCE_FASTA" ]; then
    echo "ERROR: Reference genome not found: $REFERENCE_FASTA"
    echo "Please edit the REFERENCE_FASTA variable in this script."
    exit 1
fi
echo "✓ Reference genome: $REFERENCE_FASTA"

# Check cleaned genomes
for genome_entry in "${CLEANED_GENOMES[@]}"; do
    genome_file="${genome_entry#*:}"  # Extract path after colon
    if [ ! -f "$genome_file" ]; then
        echo "ERROR: Genome file not found: $genome_file"
        echo "Please edit the CLEANED_GENOMES array in this script."
        exit 1
    fi
    echo "✓ Genome: $genome_file"
done

echo ""
echo "All input files found. Starting pipeline..."
echo ""

# Build genome arguments
GENOME_ARGS=""
for genome in "${CLEANED_GENOMES[@]}"; do
    GENOME_ARGS="$GENOME_ARGS --genomes $genome"
done

# Run with --all flag
python3 unified_pipeline.py \
    $GENOME_ARGS \
    --pangenome-reference-fasta "$REFERENCE_FASTA" \
    --all \
    --quast-threads "$THREADS" \
    --busco-threads "$THREADS" \
    --pangenome-threads "$THREADS" \
    --pggb-threads "$THREADS" \
    --pggb-segment "$PGGB_SEGMENT" \
    --pggb-pid "$PGGB_IDENTITY" \
    --pggb-passes "$PGGB_PASSES" \
    --busco-lineage "$BUSCO_LINEAGE" \
    --expected-size "$EXPECTED_SIZE" \
    --expected-gc "$EXPECTED_GC"

echo ""
echo "========================================="
echo "PIPELINE COMPLETE!"
echo "========================================="
echo ""
echo "End time: $(date)"
echo ""
echo "Check outputs in:"
echo "  - results/reports/"
echo "  - results/pggb/"
echo "  - results/minigraph_cactus/"
echo "  - results/quast_pre/, results/quast_post/"
echo "  - results/busco_pre/, results/busco_post/"
echo ""
