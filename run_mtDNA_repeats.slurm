#!/bin/bash
#SBATCH --job-name=mtDNA_repeats
#SBATCH --output=slurm_logs/mtDNA_repeats_%j.out
#SBATCH --error=slurm_logs/mtDNA_repeats_%j.err
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT
#SBATCH --mail-user=quarecute@gmail.com
#SBATCH --chdir=/home/gdzhobadze/projects/MitoRepeatability/Repeatability
#SBATCH --export=ALL

# =============================================================================
# SLURM Submission Script for mtDNA Repeat Analysis
# 
# Script: run_mtDNA_repeats.slurm
# Description: Run repeat analysis for all positions in human mtDNA
#              with all alternative alleles using 32 CPU cores and 100GB RAM
# 
# Usage: sbatch run_mtDNA_repeats.slurm
# 
# Author: Bioinformatician
# Date: $(date)
# =============================================================================

echo "================================================================================"
echo "mtDNA Repeat Analysis - SLURM Job"
echo "================================================================================"
echo "Job ID:            $SLURM_JOB_ID"
echo "Job Name:          $SLURM_JOB_NAME"
echo "Node:              $SLURM_JOB_NODELIST"
echo "Partition:         $SLURM_JOB_PARTITION"
echo "CPU Cores:         $SLURM_CPUS_PER_TASK"
echo "Memory:            $SLURM_MEM_PER_NODE"
echo "Time Limit:        $SLURM_TIMELIMIT seconds"
echo "Working Directory: $SLURM_SUBMIT_DIR"
echo "Submit Time:       $(date)"
echo "================================================================================"

# =============================================================================
# 1. MODULE LOADING
# =============================================================================
echo "1. Loading required modules..."

# Очистка модулей
module purge

module load python/3.9.0          # Python 3.9
module load gcc/11.2.0            # Компилятор GCC
module load openmpi/4.1.1         # MPI

echo "   Loaded modules:"
module list

# =============================================================================
# 2. ENVIRONMENT SETUP
# =============================================================================
echo ""
echo "2. Setting up environment..."

source /home/gdzhobadze/projects/MitoRepeatability/Repeatability/venv/bin/activate

# Создание директорий для логов
mkdir -p slurm_logs
mkdir -p logs
mkdir -p data/2_derived/01KP

# Проверка Python окружения
echo "   Python version: $(python3 --version)"
echo "   Python path: $(which python3)"

# =============================================================================
# 3. DEPENDENCY CHECK
# =============================================================================
echo ""
echo "3. Checking dependencies..."

# Проверка наличия необходимых Python пакетов
python3 -c "
import sys
required = {'numpy', 'pandas', 'Bio', 'tqdm'}
missing = []
for pkg in required:
    try:
        if pkg == 'Bio':
            __import__('Bio')
        else:
            __import__(pkg)
    except ImportError:
        missing.append(pkg)
if missing:
    print(f'Missing packages: {missing}')
    sys.exit(1)
else:
    print('All required packages are installed')
"

if [ $? -ne 0 ]; then
    echo "   Installing missing packages..."
    pip3 install --user numpy pandas biopython tqdm
fi

# =============================================================================
# 4. DATA VALIDATION
# =============================================================================
echo ""
echo "4. Validating input data..."

# Проверка FASTA файла
FASTA_FILE="data/1_raw/Homo_sapients.mtDNA.fasta"

if [ ! -f "$FASTA_FILE" ]; then
    echo "   ERROR: FASTA file not found: $FASTA_FILE"
    echo "   Download from: https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=fasta"
    exit 1
fi

# Проверка длины последовательности
SEQ_LENGTH=$(grep -v '^>' "$FASTA_FILE" | tr -d '\n' | wc -c)
echo "   FASTA file: $FASTA_FILE"
echo "   Sequence length: $SEQ_LENGTH bp"

# Проверка скрипта
SCRIPT_FILE="run_all_variants.py"
if [ ! -f "$SCRIPT_FILE" ]; then
    echo "   ERROR: Script not found: $SCRIPT_FILE"
    exit 1
fi
echo "   Script file: $SCRIPT_FILE"

# =============================================================================
# 5. RESOURCE MONITORING SETUP
# =============================================================================
echo ""
echo "5. Setting up resource monitoring..."

# Создаем файл для мониторинга ресурсов
MONITOR_FILE="logs/resource_usage_${SLURM_JOB_ID}.csv"
echo "timestamp,cpu_percent,mem_percent,mem_used_gb,mem_available_gb" > "$MONITOR_FILE"

# Функция для мониторинга ресурсов
monitor_resources() {
    while true; do
        TIMESTAMP=$(date +"%Y-%m-%d %H:%M:%S")
        CPU_PERCENT=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)
        MEM_TOTAL=$(free -g | awk '/^Mem:/ {print $2}')
        MEM_USED=$(free -g | awk '/^Mem:/ {print $3}')
        MEM_AVAIL=$(free -g | awk '/^Mem:/ {print $7}')
        MEM_PERCENT=$(echo "scale=2; $MEM_USED * 100 / $MEM_TOTAL" | bc)
        
        echo "$TIMESTAMP,$CPU_PERCENT,$MEM_PERCENT,$MEM_USED,$MEM_AVAIL" >> "$MONITOR_FILE"
        sleep 300  # Записываем каждые 5 минут
    done
}

# Запускаем мониторинг в фоновом режиме
monitor_resources &
MONITOR_PID=$!

# =============================================================================
# 6. EXECUTION OF MAIN ANALYSIS
# =============================================================================
echo ""
echo "6. Starting main analysis..."
echo "   Start time: $(date)"
echo "================================================================================"

# Устанавливаем переменные окружения для оптимальной работы
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Запускаем основной скрипт с указанием всех параметров
python3 "$SCRIPT_FILE" \
    --fasta "$FASTA_FILE" \
    --outdir "data/2_derived/01KP" \
    --max-workers 32 \
    --skip-existing \
    --log-file "logs/processing_summary_$(date +%Y%m%d_%H%M%S).txt"

ANALYSIS_EXIT_CODE=$?

echo ""
echo "================================================================================"
echo "7. Analysis completed"
echo "================================================================================"
echo "   Exit code: $ANALYSIS_EXIT_CODE"
echo "   End time: $(date)"

# Останавливаем мониторинг ресурсов
kill $MONITOR_PID 2>/dev/null
wait $MONITOR_PID 2>/dev/null

# =============================================================================
# 8. POST-PROCESSING AND CLEANUP
# =============================================================================
echo ""
echo "8. Post-processing results..."

# Подсчет результатов
if [ -d "data/2_derived/01KP" ]; then
    RESULT_COUNT=$(find "data/2_derived/01KP" -name "01KP.*.txt" -type f | wc -l)
    TOTAL_SIZE=$(find "data/2_derived/01KP" -name "01KP.*.txt" -type f -exec stat -c%s {} + | awk '{sum+=$1} END {print sum}')
    
    echo "   Result files: $RESULT_COUNT"
    echo "   Total size: $(echo "scale=2; $TOTAL_SIZE / (1024*1024*1024)" | bc) GB"
    
    # Создание сводного отчета
    cat > "logs/job_summary_${SLURM_JOB_ID}.txt" << EOF
========================================
mtDNA Repeat Analysis - Job Summary
========================================
Job Information:
  Job ID: $SLURM_JOB_ID
  Job Name: $SLURM_JOB_NAME
  Partition: $SLURM_JOB_PARTITION
  Node: $SLURM_JOB_NODELIST
  CPU Cores: $SLURM_CPUS_PER_TASK
  Memory: $SLURM_MEM_PER_NODE
  Time Limit: $SLURM_TIMELIMIT

Timing:
  Start Time: $(grep "Start time" slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.out | tail -1)
  End Time: $(date)
  Duration: $(($(date +%s) - $(grep "Start time" slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.out | tail -1 | date +%s -f -))) seconds

Results:
  Total Files: $RESULT_COUNT
  Total Size: $(echo "scale=2; $TOTAL_SIZE / (1024*1024*1024)" | bc) GB
  Output Directory: data/2_derived/01KP

Status:
  Exit Code: $ANALYSIS_EXIT_CODE
  Completion: $(if [ $ANALYSIS_EXIT_CODE -eq 0 ]; then echo "SUCCESS"; else echo "FAILED"; fi)
========================================
EOF
else
    echo "   WARNING: Results directory not found"
fi

# Проверка на ошибки
if [ -f "slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.err" ]; then
    ERROR_COUNT=$(grep -c "ERROR\|Error\|error" "slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.err")
    if [ $ERROR_COUNT -gt 0 ]; then
        echo "   WARNING: Found $ERROR_COUNT error(s) in error log"
        echo "   Check: slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.err"
    fi
fi

# =============================================================================
# 9. FINALIZATION
# =============================================================================
echo ""
echo "9. Finalizing job..."

# Создание файла завершения
COMPLETION_FILE="logs/job_${SLURM_JOB_ID}_completed.txt"
echo "Job $SLURM_JOB_ID completed at $(date)" > "$COMPLETION_FILE"
echo "Exit code: $ANALYSIS_EXIT_CODE" >> "$COMPLETION_FILE"

echo ""
echo "================================================================================"
echo "Job completed successfully!"
echo "================================================================================"
echo ""
echo "Output files location:"
echo "  - Results: data/2_derived/01KP/"
echo "  - Logs: slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.out"
echo "  - Errors: slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.err"
echo "  - Resource usage: $MONITOR_FILE"
echo "  - Job summary: logs/job_summary_${SLURM_JOB_ID}.txt"
echo ""
echo "To check results:"
echo "  ls -lh data/2_derived/01KP/ | head -20"
echo "  tail -f slurm_logs/mtDNA_repeats_${SLURM_JOB_ID}.out"
echo ""
echo "================================================================================"

exit $ANALYSIS_EXIT_CODE
