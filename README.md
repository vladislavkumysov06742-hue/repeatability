# Repeatability analysis of human mtDNA

We define **repeatability** of a nucleotide position as the maximum effective length of a direct repeat that overlaps that position (allowing a small fraction of mismatches). This property can help understand how single‑nucleotide variants (SNVs) affect the formation of deletions and other structural rearrangements.

This repository contains a full pipeline to:
- Compute all approximate repeats covering every position and every alternative nucleotide.
- Clean the raw data (remove self‑repeats, nested repeats, add GC content and H‑bonds).
- Perform statistical tests (MAD Z‑score) to identify positions where alternative alleles significantly change repeatability.
- Generate a wide range of publication‑ready figures (Manhattan plots, heatmaps, scatter plots, boxplots, sliding windows, etc.).
- Compare two genomes (e.g., human vs. pika) using the same metrics.

## Repository structure

```
repeatability/
├── data/                     # Input data and raw/derived results
├── figures/                  # Generated figures (optional)
├── logs/                     # SLURM job logs
├── notebooks/                # Standalone analysis scripts (01A – 01E, analysis.py, etc.)
├── repeat_analysis/          # The main Python package (cli, stats, plot, clean, etc.)
├── scripts/                  # Legacy/utility scripts
├── slurm/                    # SLURM submission script
├── .gitignore
├── LICENSE
├── README.md
└── requirements.txt
```

## Requirements

- Python ≥ 3.8
- Packages listed in `requirements.txt` (pandas, numpy, scipy, matplotlib, seaborn, statsmodels, numba, tqdm, biopython)
- Optional but recommended: `statsmodels` for FDR correction, `numba` for speed.

Install with:

```bash
pip install -r requirements.txt
```

## Data

Place the reference mitochondrial genome (FASTA) in `data/1_raw/Homo_sapiens.mtDNA.fasta`.  
Example: `NC_012920.1`

## Usage

All functionality is accessible via the command‑line interface `repeat_analysis.cli`.  
Run `python -m repeat_analysis.cli --help` to see available commands.

### Main commands

| Command   | Description |
|-----------|-------------|
| `clean`   | Filter raw output: remove self‑repeats, nested repeats, add GC/H‑bonds. |
| `stats`   | Compute all statistical tests and generate ~20 figures. |
| `plot`    | Quick ΔR bar plot and forest plot for six model SNPs. |
| `all`     | Run clean + stats in one go (but not the raw repeat generation). |
| `compare` | Compare two genomes (e.g., human vs. pika) – produces comparison plots. |

### Typical workflow (from scratch)

```bash
# 1. Generate raw repeat files (using the notebook script, not the CLI)
python notebooks/01B.KP.RunTheFunctionOnWholeMtDna.py \
    --fasta data/1_raw/Homo_sapiens.mtDNA.fasta \
    --outdir data/2_derived/01KP \
    --workers 32

# 2. Clean
python -m repeat_analysis.cli clean \
    --input_dir data/2_derived/01KP \
    --output_dir data/3_results/cleaned \
    --reference data/1_raw/Homo_sapiens.mtDNA.fasta \
    --workers 32

# 3. Full statistical analysis and figures
python -m repeat_analysis.cli stats \
    --input_dir data/3_results/cleaned \
    --output_dir data/3_results/stats \
    --reference data/1_raw/Homo_sapiens.mtDNA.fasta \
    --workers 32 \
    --window_size 200

# 4. (Optional) Simple plots for the six model SNPs
python -m repeat_analysis.cli plot \
    --input_dir data/3_results/cleaned \
    --output_dir data/3_results/plots \
    --reference data/1_raw/Homo_sapiens.mtDNA.fasta
```

### Running on a SLURM cluster

A ready‑to‑use submission script is provided in `slurm/run_all.slurm`.  
Adjust the paths, partition, email, and module loads, then run:

```bash
sbatch slurm/run_all.slurm --fasta /path/to/mtDNA.fasta
```

Optionally, provide a second FASTA for cross‑species comparison:

```bash
sbatch slurm/run_all.slurm --fasta human.fasta --pika_fasta pika.fasta
```

The script will:
- generate raw repeat files (step 1),
- clean them,
- run all notebook‑derived scripts (`01C`, `01D`, `01E`, `analysis.py`, `prepare_key_positions.py`, `scatter_plot.py`),
- and finally run the advanced `repeat_analysis.cli stats` and `repeat_analysis.cli plot`.

### Comparing two genomes

```bash
python -m repeat_analysis.cli compare \
    --input_dir1 data/human/cleaned \
    --reference1 data/human/ref.fasta \
    --input_dir2 data/pika/cleaned \
    --reference2 data/pika/ref.fasta \
    --name1 Human --name2 Pika \
    --output_dir comparison_results
```

## Analysis scripts

The `notebooks/` directory contains several standalone Python scripts that can be used independently or together:

| Script | Description |
|--------|-------------|
| `01A.KP.PrepareTheFunction.py` | Core function to compute repeats for a single position/nucleotide. |
| `01B.KP.RunTheFunctionOnWholeMtDna.py` | Parallel wrapper to run `01A` for the whole mtDNA genome. |
| `01C.KP.TryRepeatabilityMetrics.py` | Compare Ref vs Alt alleles, generate mutation tables and delta‑repeatability statistics. |
| `01D.KP.StatisticsAndPlots.py` | Generate histograms, heatmaps, boxplots, and GC‑content analysis from cleaned files. |
| `01E.KP.Preliminary_Analysis.py` | Extract longest perfect repeat length for specific positions. |
| `analysis.py` | Full ΔR analysis, GC context, and publication‑ready figures (scatter, heatmaps, mutation‑type breakdown). |
| `correlation_analysis.py` | Cross‑species correlation (sliding window and scatter). |
| `prepare_key_positions.py` | Focused analysis on a predefined set of key mtDNA positions (8251, 8472, 8473, 12705, 16223). |
| `scatter_plot.py` | Quick scatter (ref length vs. |ΔR|) and histogram of ΔR. |

All these scripts can be run directly, for example:

```bash
python notebooks/analysis.py --input_dir data/3_results/cleaned --fasta data/1_raw/Homo_sapiens.mtDNA.fasta --outdir my_results
```

## Output

All results are saved in the specified output directories. For the `stats` command, typical outputs include:

- `manhattan_metrics.png` – Manhattan plots for 7 repeatability metrics.
- `coverage_plot.png` – weighted coverage along the genome.
- `heatmap_max_perfect.png` – heatmap of Δ max perfect repeat.
- `scatter_delta_perfect.png` – scatter plot with regression.
- `boxplot_99percentile_*.png` – 99th percentile of |Δ| for transitions/transversions.
- `windowed_*` – sliding window bar plots (if `--window_size > 0`).
- `stats_run.log` – detailed log of all tests.

The notebook‑derived scripts produce their own output directories as specified by their `--output` or `--output_dir` arguments.

## Citation

If you use this pipeline, please cite the original concept and this implementation.  
(Add your own citation when published.)
