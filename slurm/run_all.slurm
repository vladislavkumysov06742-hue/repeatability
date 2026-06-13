#!/bin/bash
#SBATCH --job-name=mtDNA_repeatability
#SBATCH --output=logs/run_all_%j.out
#SBATCH --error=logs/run_all_%j.err
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@example.com

# ==================================================
#   Full pipeline: clean → stats → plots
# ==================================================

set -e  # exit on error

# Load modules (adjust to your cluster)
# module load python/3.9
# module load gcc

# Change to the repository root
cd /path/to/repeatability || exit 1

# Create necessary directories
mkdir -p logs data/2_derived/01KP data/3_results

# Setup Python environment (optional)
# python -m venv venv
# source venv/bin/activate
# pip install -r requirements.txt

# Set variables
RAW_FASTA="data/1_raw/Homo_sapiens.mtDNA.fasta"
RAW_REPEAT_DIR="data/2_derived/01KP"
CLEAN_DIR="data/3_results/cleaned"
STATS_DIR="data/3_results/stats"
PLOT_DIR="data/3_results/plots"

echo "========================================"
echo "Starting full repeatability pipeline"
echo "========================================"

# Step 1: Generate raw repeat data (if not already done)
# This step is time‑consuming; run only once.
# Uncomment if you need to run the repeat finder on all positions.
# echo "Step 1: Generating raw repeats for all positions..."
# python -m repeat_analysis.cli generate \
#     --fasta "$RAW_FASTA" \
#     --outdir "$RAW_REPEAT_DIR" \
#     --max-workers 32

# Step 2: Clean the raw data (filter self‑repeats, nested repeats, add GC/H-bonds)
echo "Step 2: Cleaning raw repeat files..."
python -m repeat_analysis.cli clean \
    --input_dir "$RAW_REPEAT_DIR" \
    --output_dir "$CLEAN_DIR" \
    --reference "$RAW_FASTA" \
    --workers 32

# Step 3: Compute full statistics and generate all scientific plots
echo "Step 3: Running statistical analysis and plots..."
python -m repeat_analysis.cli stats \
    --input_dir "$CLEAN_DIR" \
    --output_dir "$STATS_DIR" \
    --reference "$RAW_FASTA" \
    --workers 32 \
    --window_size 200

# Step 4: (Optional) Generate simple ΔR bar plots for the 6 model SNPs
echo "Step 4: Generating delta‑R plots for selected SNPs..."
python -m repeat_analysis.cli plot \
    --input_dir "$CLEAN_DIR" \
    --output_dir "$PLOT_DIR" \
    --reference "$RAW_FASTA"

echo "========================================"
echo "All steps completed successfully"
echo "Results in: $STATS_DIR and $PLOT_DIR"
echo "========================================"
