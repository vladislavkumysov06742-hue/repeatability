#!/bin/bash
#SBATCH --job-name=notebooks_then_repeat
#SBATCH --output=logs/notebooks_then_repeat_%j.out
#SBATCH --error=logs/notebooks_then_repeat_%j.err
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=48:00:00

# ----------------------------------------------------------
# Parse arguments: --fasta (required), --pika_fasta (optional)
# ----------------------------------------------------------
FASTA_PATH=""
PIKA_FASTA=""
while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta) FASTA_PATH="$2"; shift 2 ;;
        --pika_fasta) PIKA_FASTA="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done
if [ -z "$FASTA_PATH" ]; then echo "ERROR: --fasta required"; exit 1; fi
if [ ! -f "$FASTA_PATH" ]; then echo "ERROR: $FASTA_PATH not found"; exit 1; fi
if [ -n "$PIKA_FASTA" ] && [ ! -f "$PIKA_FASTA" ]; then echo "ERROR: $PIKA_FASTA not found"; exit 1; fi

set -e
cd /path/to/repeatability || exit 1

# Create directories
mkdir -p logs data/2_derived/01KP data/3_results/{cleaned,stats,notebook_outputs,correlation}

RAW_DIR="data/2_derived/01KP"
CLEAN_DIR="data/3_results/cleaned"
STATS_DIR="data/3_results/stats"
NOTEBOOK_DIR="data/3_results/notebook_outputs"
CORR_DIR="data/3_results/correlation"

# ----------------------------------------------------------
# 1. Generate raw repeat files using 01B (notebook version)
# ----------------------------------------------------------
echo "=========================================="
echo "Step 1: 01B – Generate raw repeats for all positions"
echo "=========================================="
python scripts/01B_run_whole_mtDNA.py \
    --fasta "$FASTA_PATH" \
    --outdir "$RAW_DIR" \
    --workers "$SLURM_CPUS_PER_TASK"

# ----------------------------------------------------------
# 2. Clean raw files using repeat_analysis clean (needed by 01D, analysis, scatter)
# ----------------------------------------------------------
echo "=========================================="
echo "Step 2: Clean raw files (repeat_analysis clean) – but not stats yet"
echo "=========================================="
python -m repeat_analysis.cli clean \
    --input_dir "$RAW_DIR" \
    --output_dir "$CLEAN_DIR" \
    --reference "$FASTA_PATH" \
    --workers "$SLURM_CPUS_PER_TASK"

# ----------------------------------------------------------
# 3. Run all notebook‑derived scripts (they now use the cleaned files)
# ----------------------------------------------------------
echo "=========================================="
echo "Step 3: 01C – Mutation analysis (uses raw files)"
echo "=========================================="
python scripts/01C.KP.TryRepeatabilityMetrics.py \
    --fasta "$FASTA_PATH" \
    --input_dir "$RAW_DIR" \
    --output "$NOTEBOOK_DIR/01C_results"

echo "=========================================="
echo "Step 4: 01D – Statistics & plots (uses cleaned files)"
echo "=========================================="
python scripts/01D.KP.StatisticsAndPlots.py \
    --input_dir "$CLEAN_DIR" \
    --fasta "$FASTA_PATH" \
    --output_dir "$NOTEBOOK_DIR/01D_figures" \
    --window_size 200

echo "=========================================="
echo "Step 5: 01E – Robust repeatability for key positions (uses raw files)"
echo "=========================================="
python scripts/01E.KP.Preliminary_Analysis.py \
    --input_dir "$RAW_DIR" \
    --positions 8251 8472 8473 12705 14798 16223 \
    --output "$NOTEBOOK_DIR/01E_metrics.csv"

echo "=========================================="
echo "Step 6: analysis.py – Full ΔR, GC context, heatmaps (uses cleaned files)"
echo "=========================================="
python scripts/analysis.py \
    --input_dir "$CLEAN_DIR" \
    --fasta "$FASTA_PATH" \
    --outdir "$NOTEBOOK_DIR/full_analysis"

echo "=========================================="
echo "Step 7: prepare_key_positions.py – Six publication figures (uses raw files)"
echo "=========================================="
python scripts/prepare_key_positions.py \
    --input_dir "$RAW_DIR" \
    --fasta "$FASTA_PATH" \
    --output_dir "$NOTEBOOK_DIR/key_positions_figures"

echo "=========================================="
echo "Step 8: scatter_plot.py – Simple scatter & histogram (uses cleaned files)"
echo "=========================================="
python scripts/scatter_plot.py \
    --input_dir "$CLEAN_DIR" \
    --fasta "$FASTA_PATH" \
    --output "$NOTEBOOK_DIR/scatter_histogram.png" \
    --csv "$NOTEBOOK_DIR/delta_R_values.csv"

echo "=========================================="
echo "Step 9: correlation_analysis.py – Cross‑species comparison (if pika FASTA provided)"
echo "=========================================="
if [ -n "$PIKA_FASTA" ]; then
    PIKA_RAW="data/pika/01KP"
    PIKA_CLEAN="data/pika/cleaned"
    echo "Generating pika repeats..."
    python scripts/01B_run_whole_mtDNA.py --fasta "$PIKA_FASTA" --outdir "$PIKA_RAW" --workers "$SLURM_CPUS_PER_TASK"
    python -m repeat_analysis.cli clean --input_dir "$PIKA_RAW" --output_dir "$PIKA_CLEAN" --reference "$PIKA_FASTA" --workers "$SLURM_CPUS_PER_TASK"
    python scripts/correlation_analysis.py \
        --human_dir "$CLEAN_DIR" \
        --pika_dir "$PIKA_CLEAN" \
        --output_dir "$CORR_DIR" \
        --window 100
else
    echo "Skipping correlation (no --pika_fasta)"
fi

# ----------------------------------------------------------
# 4. Finally, run the full repeat_analysis stats (advanced part)
# ----------------------------------------------------------
echo "=========================================="
echo "Step 10: repeat_analysis stats – Advanced statistics and figures"
echo "=========================================="
python -m repeat_analysis.cli stats \
    --input_dir "$CLEAN_DIR" \
    --output_dir "$STATS_DIR" \
    --reference "$FASTA_PATH" \
    --workers "$SLURM_CPUS_PER_TASK" \
    --window_size 200

echo "=========================================="
echo "Step 11: repeat_analysis plot – Quick ΔR bar plots for 6 SNPs"
echo "=========================================="
python -m repeat_analysis.cli plot \
    --input_dir "$CLEAN_DIR" \
    --output_dir "$STATS_DIR/plots" \
    --reference "$FASTA_PATH"

echo "=========================================="
echo "ALL STEPS COMPLETED"
echo "Notebook outputs: $NOTEBOOK_DIR"
echo "Advanced stats: $STATS_DIR"
if [ -n "$PIKA_FASTA" ]; then echo "Correlation: $CORR_DIR"; fi
echo "=========================================="
