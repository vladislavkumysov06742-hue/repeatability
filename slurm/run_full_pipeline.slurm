#!/bin/bash
#SBATCH --job-name=full_mtDNA_repeat
#SBATCH --output=logs/full_pipeline_%j.out
#SBATCH --error=logs/full_pipeline_%j.err
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@example.com

# ----------------------------------------------------------
# Parse arguments: --fasta (required), --pika_fasta (optional)
# ----------------------------------------------------------
FASTA_PATH=""
PIKA_FASTA=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta)
            FASTA_PATH="$2"
            shift 2
            ;;
        --pika_fasta)
            PIKA_FASTA="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: sbatch run_full_pipeline.slurm --fasta <path> [--pika_fasta <path>]"
            exit 1
            ;;
    esac
done

if [ -z "$FASTA_PATH" ]; then
    echo "ERROR: --fasta is required"
    exit 1
fi
if [ ! -f "$FASTA_PATH" ]; then
    echo "ERROR: FASTA file not found: $FASTA_PATH"
    exit 1
fi
if [ -n "$PIKA_FASTA" ] && [ ! -f "$PIKA_FASTA" ]; then
    echo "ERROR: Pika FASTA not found: $PIKA_FASTA"
    exit 1
fi

# ----------------------------------------------------------
# Environment setup (adjust to your cluster)
# ----------------------------------------------------------
set -e  # exit on error

# Load modules (if needed)
# module load python/3.9
# module load gcc

# Change to repository root
cd /path/to/repeatability || exit 1

# Create directories
mkdir -p logs data/2_derived/01KP data/3_results/{cleaned,stats,plots,extra_figures}

# Activate virtual environment if used
# source venv/bin/activate
# pip install -r requirements.txt

# ----------------------------------------------------------
# Define input/output directories for different steps
# ----------------------------------------------------------
RAW_DIR="data/2_derived/01KP"
CLEAN_DIR="data/3_results/cleaned"
STATS_DIR="data/3_results/stats"
EXTRA_DIR="data/3_results/extra_figures"

echo "=========================================="
echo "Full mtDNA repeatability pipeline"
echo "FASTA: $FASTA_PATH"
if [ -n "$PIKA_FASTA" ]; then
    echo "Pika FASTA: $PIKA_FASTA"
    echo "Cross‑species correlation will be performed"
fi
echo "=========================================="

# ----------------------------------------------------------
# PART 1: Run repeat_analysis CLI (generate, clean, stats)
# ----------------------------------------------------------
echo "--- Part 1: repeat_analysis pipeline ---"

# 1.1 Generate raw repeat files (all positions, all nucleotides)
echo "Step 1.1: Generating raw repeats..."
python -m repeat_analysis.cli generate \
    --fasta "$FASTA_PATH" \
    --outdir "$RAW_DIR" \
    --workers "$SLURM_CPUS_PER_TASK"

# 1.2 Clean raw files (remove self‑repeats, add GC/H-bonds)
echo "Step 1.2: Cleaning repeats..."
python -m repeat_analysis.cli clean \
    --input_dir "$RAW_DIR" \
    --output_dir "$CLEAN_DIR" \
    --reference "$FASTA_PATH" \
    --workers "$SLURM_CPUS_PER_TASK"

# 1.3 Run full statistical analysis (all figures + reports)
echo "Step 1.3: Running statistics and plots..."
python -m repeat_analysis.cli stats \
    --input_dir "$CLEAN_DIR" \
    --output_dir "$STATS_DIR" \
    --reference "$FASTA_PATH" \
    --workers "$SLURM_CPUS_PER_TASK" \
    --window_size 200

# 1.4 Generate quick ΔR plots for the six model SNPs
echo "Step 1.4: Generating SNP delta‑R plots..."
python -m repeat_analysis.cli plot \
    --input_dir "$CLEAN_DIR" \
    --output_dir "$STATS_DIR/plots" \
    --reference "$FASTA_PATH"

# ----------------------------------------------------------
# PART 2: Run the converted notebook scripts
# ----------------------------------------------------------
echo "--- Part 2: Running additional analyses from notebooks ---"

# 2.1 01C: Mutation analysis (requires raw files)
echo "Step 2.1: 01C – Mutation repeatability metrics..."
python scripts/01C.KP.TryRepeatabilityMetrics.py \
    --fasta "$FASTA_PATH" \
    --input_dir "$RAW_DIR" \
    --output "$EXTRA_DIR/01C_results" \
    --max_degradation 20

# 2.2 01D: Statistics and plots from cleaned files
echo "Step 2.2: 01D – Statistics and plots (histogram, heatmap, trend, GC vs length)..."
python scripts/01D.KP.StatisticsAndPlots.py \
    --input_dir "$CLEAN_DIR" \
    --fasta "$FASTA_PATH" \
    --output_dir "$EXTRA_DIR/01D_figures" \
    --window_size 200 \
    --dpi 150

# 2.3 01E: Robust repeatability metric for specific positions
echo "Step 2.3: 01E – Robust repeatability metric for key positions..."
python scripts/01E.KP.Preliminary_Analysis.py \
    --input_dir "$RAW_DIR" \
    --positions 8251 8472 8473 12705 14798 16223 \
    --output "$EXTRA_DIR/01E_repeatability_metrics.csv"

# 2.4 analysis.py (full ΔR + GC context + heatmaps + text report)
echo "Step 2.4: analysis.py – Comprehensive analysis..."
python scripts/analysis.py \
    --input_dir "$CLEAN_DIR" \
    --fasta "$FASTA_PATH" \
    --outdir "$EXTRA_DIR/full_analysis"

# 2.5 prepare_key_positions.py (six figures for key positions)
echo "Step 2.5: prepare_key_positions.py – Key positions visualisation..."
python scripts/prepare_key_positions.py \
    --input_dir "$RAW_DIR" \
    --fasta "$FASTA_PATH" \
    --output_dir "$EXTRA_DIR/key_positions_figures"

# 2.6 scatter_plot.py (basic ΔR scatter and histogram)
echo "Step 2.6: scatter_plot.py – Scatter plot for all positions..."
python scripts/scatter_plot.py \
    --input_dir "$CLEAN_DIR" \
    --fasta "$FASTA_PATH" \
    --output "$EXTRA_DIR/scatter_histogram.png" \
    --csv "$EXTRA_DIR/delta_R_values.csv"

# 2.7 correlation_analysis.py (cross‑species, if pika FASTA provided)
if [ -n "$PIKA_FASTA" ]; then
    echo "Step 2.7: correlation_analysis.py – Comparing human and pika..."
    # First generate repeats for pika (clean then compare)
    PIKA_RAW="data/pika/01KP"
    PIKA_CLEAN="data/pika/cleaned"
    python -m repeat_analysis.cli generate --fasta "$PIKA_FASTA" --outdir "$PIKA_RAW" --workers "$SLURM_CPUS_PER_TASK"
    python -m repeat_analysis.cli clean --input_dir "$PIKA_RAW" --output_dir "$PIKA_CLEAN" --reference "$PIKA_FASTA" --workers "$SLURM_CPUS_PER_TASK"
    python scripts/correlation_analysis.py \
        --human_dir "$CLEAN_DIR" \
        --pika_dir "$PIKA_CLEAN" \
        --output_dir "$EXTRA_DIR/correlation_results" \
        --window 100
else
    echo "Step 2.7: Skipping correlation analysis (no --pika_fasta)"
fi

# ----------------------------------------------------------
# Final summary
# ----------------------------------------------------------
echo "=========================================="
echo "ALL ANALYSES COMPLETED SUCCESSFULLY"
echo "Main results: $STATS_DIR"
echo "Additional figures: $EXTRA_DIR"
echo "=========================================="
