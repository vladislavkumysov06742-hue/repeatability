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
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your.email@example.com  # ЗАМЕНИТЕ на ваш email

# ============================================================================
# SLURM Script for Running mtDNA Repeat Analysis
# 
# This script runs the repeat analysis for all positions in human mtDNA
# with all alternative alleles using 32 CPU cores and 100GB RAM.
# 
# Usage: sbatch run_mtDNA_repeats_slurm.sh
# ============================================================================

echo "================================================================================"
echo "Starting mtDNA Repeat Analysis"
echo "================================================================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_JOB_NODELIST"
echo "Partition: $SLURM_JOB_PARTITION"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Memory: $SLURM_MEM_PER_NODE"
echo "Time Limit: $SLURM_TIMELIMIT"
echo "================================================================================"

# Загружаем необходимые модули (зависит от кластера)
# Если на кластере используется модульная система, раскомментируйте нужные строки:

# module load python/3.9.0
# module load gcc/11.2.0
# module load openmpi/4.1.1

# Проверяем, есть ли нужные модули Python
# Если нет модульной системы, используем conda или установленные пакеты

# Настройка Python окружения
PYTHON_PATH="/usr/bin/python3"  # Измените на путь к вашему Python 3
# или используйте conda:
# source activate bioinformatics_env

# Создаем директорию для логов
mkdir -p slurm_logs

# Переходим в директорию со скриптами
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd "$SCRIPT_DIR" || exit 1

echo "Current directory: $(pwd)"
echo "Script directory: $SCRIPT_DIR"

# Проверяем существование скрипта
MAIN_SCRIPT="scripts/run_all_variants.py"
if [ ! -f "$MAIN_SCRIPT" ]; then
    echo "ERROR: Main script not found: $MAIN_SCRIPT"
    echo "Looking for script in:"
    find . -name "run_all_variants.py" -type f 2>/dev/null
    exit 1
fi

echo "Main script found: $MAIN_SCRIPT"

# Проверяем существование FASTA файла
FASTA_FILE="data/1_raw/Homo_sapients.mtDNA.fasta"
if [ ! -f "$FASTA_FILE" ]; then
    echo "ERROR: FASTA file not found: $FASTA_FILE"
    echo "Please download the mtDNA FASTA from:"
    echo "https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1?report=fasta"
    echo "and save it to $FASTA_FILE"
    exit 1
fi

echo "FASTA file found: $FASTA_FILE"
echo "Sequence length: $(grep -v '^>' "$FASTA_FILE" | tr -d '\n' | wc -c) bp"

# Создаем выходные директории
mkdir -p data/2_derived/01KP
mkdir -p logs

# Генерируем уникальный префикс для логов на основе времени
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="logs/run_all_variants_${TIMESTAMP}.log"
SUMMARY_FILE="logs/processing_summary_${TIMESTAMP}.txt"

echo "Log file: $LOG_FILE"
echo "Summary file: $SUMMARY_FILE"

# Устанавливаем переменные окружения для оптимальной работы Python
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=1  # Отключаем OpenMP параллелизм внутри процессов
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Функция для вывода информации о системе
print_system_info() {
    echo "================================================================================"
    echo "System Information"
    echo "================================================================================"
    echo "Hostname: $(hostname)"
    echo "OS: $(uname -s)"
    echo "Kernel: $(uname -r)"
    echo "CPU threads available: $(nproc)"
    echo "Total memory: $(free -h | awk '/^Mem:/ {print $2}')"
    echo "Python version: $($PYTHON_PATH --version 2>&1)"
    echo "Python path: $(which $PYTHON_PATH)"
    echo "================================================================================"
}

# Функция для установки Python зависимостей
setup_python_environment() {
    echo "Setting up Python environment..."
    
    # Проверяем наличие virtualenv/conda или устанавливаем зависимости глобально
    if [ -f "requirements.txt" ]; then
        echo "Found requirements.txt"
        
        # Создаем виртуальное окружение если его нет
        if [ ! -d "venv" ]; then
            echo "Creating virtual environment..."
            $PYTHON_PATH -m venv venv
        fi
        
        # Активируем виртуальное окружение
        if [ -f "venv/bin/activate" ]; then
            source venv/bin/activate
            echo "Virtual environment activated"
            
            # Обновляем pip
            pip install --upgrade pip
            
            # Устанавливаем зависимости
            echo "Installing Python dependencies..."
            pip install -r requirements.txt
            
            # Дополнительные пакеты для производительности
            pip install numpy pandas biopython tqdm
        else
            echo "Warning: Virtual environment not found, using system Python"
        fi
    else
        echo "Warning: requirements.txt not found, trying to install basic packages..."
        pip install numpy pandas biopython tqdm
    fi
    
    echo "Python environment setup complete"
}

# Функция для запуска анализа с мониторингом
run_analysis() {
    echo "================================================================================"
    echo "Starting Repeat Analysis"
    echo "================================================================================"
    
    local START_TIME=$(date +%s)
    
    # Запускаем основной скрипт
    # Используем все доступные CPU (32) для параллельной обработки
    $PYTHON_PATH "$MAIN_SCRIPT" \
        --fasta "$FASTA_FILE" \
        --outdir data/2_derived/01KP \
        --max-workers 32 \
        --skip-existing \
        --log-file "$SUMMARY_FILE" \
        2>&1 | tee "$LOG_FILE"
    
    local EXIT_CODE=${PIPESTATUS[0]}
    local END_TIME=$(date +%s)
    local DURATION=$((END_TIME - START_TIME))
    
    echo "Analysis completed with exit code: $EXIT_CODE"
    echo "Total runtime: $DURATION seconds ($(($DURATION / 3600))h $((($DURATION % 3600) / 60))m $(($DURATION % 60))s)"
    
    return $EXIT_CODE
}

# Функция для пост-обработки результатов
post_processing() {
    echo "================================================================================"
    echo "Post-processing Results"
    echo "================================================================================"
    
    local RESULTS_DIR="data/2_derived/01KP"
    
    if [ -d "$RESULTS_DIR" ]; then
        echo "Counting result files..."
        local FILE_COUNT=$(find "$RESULTS_DIR" -name "01KP.*.txt" -type f | wc -l)
        echo "Total result files: $FILE_COUNT"
        
        echo "Calculating total size..."
        local TOTAL_SIZE=$(find "$RESULTS_DIR" -name "01KP.*.txt" -type f -exec stat -c%s {} + | awk '{sum+=$1} END {print sum}')
        echo "Total data size: $(numfmt --to=iec-i --suffix=B $TOTAL_SIZE)"
        
        # Создаем сводный файл
        echo "Creating summary report..."
        cat > "logs/analysis_summary_${TIMESTAMP}.txt" << EOF
mtDNA Repeat Analysis Summary
========================================
Timestamp: $(date)
Job ID: $SLURM_JOB_ID
Host: $(hostname)
Runtime: ${DURATION}s
Total files processed: $FILE_COUNT
Total data size: $(numfmt --to=iec-i --suffix=B $TOTAL_SIZE)
Output directory: $RESULTS_DIR
Log file: $LOG_FILE
EOF
        
        # Проверяем наличие файла с ошибками
        if [ -f "$LOG_FILE" ] && grep -q "ERROR" "$LOG_FILE"; then
            local ERROR_COUNT=$(grep -c "ERROR" "$LOG_FILE")
            echo "WARNING: Found $ERROR_COUNT errors in log file"
        fi
        
    else
        echo "WARNING: Results directory not found: $RESULTS_DIR"
    fi
}

# Функция для мониторинга ресурсов
monitor_resources() {
    echo "================================================================================"
    echo "Resource Monitoring"
    echo "================================================================================"
    
    # Создаем файл для мониторинга
    MONITOR_FILE="logs/resource_monitor_${TIMESTAMP}.csv"
    echo "timestamp,cpu_percent,mem_percent,mem_used_gb,mem_available_gb" > "$MONITOR_FILE"
    
    # Запускаем мониторинг в фоновом режиме
    (
        while true; do
            TIMESTAMP_MON=$(date +"%Y-%m-%d %H:%M:%S")
            CPU_PERCENT=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}' | cut -d'%' -f1)
            MEM_INFO=$(free -m | awk '/^Mem:/ {print $3, $7}')
            MEM_USED=$(echo "$MEM_INFO" | awk '{print $1}')
            MEM_AVAIL=$(echo "$MEM_INFO" | awk '{print $2}')
            MEM_PERCENT=$(echo "scale=2; $MEM_USED / ($MEM_USED + $MEM_AVAIL) * 100" | bc)
            
            echo "$TIMESTAMP_MON,$CPU_PERCENT,$MEM_PERCENT,$(echo "scale=2; $MEM_USED / 1024" | bc),$(echo "scale=2; $MEM_AVAIL / 1024" | bc)" >> "$MONITOR_FILE"
            
            sleep 60  # Обновляем каждую минуту
        done
    ) &
    
    MONITOR_PID=$!
    
    # Останавливаем мониторинг при завершении скрипта
    trap "kill $MONITOR_PID 2>/dev/null" EXIT
}

# Основной процесс выполнения
main() {
    echo "================================================================================"
    echo "Job Execution Started at: $(date)"
    echo "================================================================================"
    
    # Выводим информацию о системе
    print_system_info
    
    # Настраиваем Python окружение
    setup_python_environment
    
    # Запускаем мониторинг ресурсов
    monitor_resources
    
    # Запускаем анализ
    run_analysis
    ANALYSIS_EXIT_CODE=$?
    
    # Пост-обработка
    post_processing
    
    # Очистка (если нужно)
    # ...
    
    echo "================================================================================"
    echo "Job Execution Finished at: $(date)"
    echo "Exit Code: $ANALYSIS_EXIT_CODE"
    echo "================================================================================"
    
    # Создаем завершающий файл
    echo "Job completed" > "logs/job_completed_${TIMESTAMP}.txt"
    
    return $ANALYSIS_EXIT_CODE
}

# Запускаем основную функцию
main
exit $?
