import sys
from pathlib import Path

# Попытка импортировать pandas — популярная библиотека для работы с таблицами
try:
    import pandas as pd
except Exception as e:
    # Если импорт не удался, печатаем сообщение и выходим с кодом 2
    # (код выхода используется здесь для диагностики при запуске из оболочки)
    print('Pandas import error:', e)
    sys.exit(2)

# Импортируем функцию `analyze_mtDNA_repeats` из локального пакета `src.analysis`.
# Эта функция выполняет основной анализ и возвращает пути к сгенерированным файлам.
try:
    from src.analysis import analyze_mtDNA_repeats
except Exception as e:
    # Если импорт не удался, сообщаем и выходим — вероятно, проблема с путями или __init__.py
    print('Import error from src.analysis:', e)
    sys.exit(2)

# Определяем базовую директорию, где лежит этот скрипт
base = Path(__file__).resolve().parents[0]

# Собираем путь к входному FASTA-файлу (относительно директории с этим скриптом)
fasta = base / '..' / 'data' / '1_raw' / 'Homo_sapients.mtDNA.fasta'
fasta = fasta.resolve()  # приводим к абсолютному пути для ясности
print('Using fasta:', fasta)

# Вызов основной функции анализа.
# Передаём путь к FASTA и параметры варианта (позиция и нуклеотиды).
try:
    out = analyze_mtDNA_repeats(
        str(fasta),
        pos=8473,
        ref_nuc='T',
        alt_nuc='C',
        output_path=str(base / '..' / 'data' / '2_derived')
    )
    # `out` ожидается как словарь/магический объект с путями к файлам
    print('analyze returned:', out)
except Exception as e:
    # Если анализ упал — печатаем ошибку и выходим с кодом 3
    print('analyze_mtDNA_repeats failed:', e)
    sys.exit(3)

# Попробуем сравнить с референсными R-результатами (если они есть в репозитории)
# `generated_all` и `generated_summary` — предполагаемые пути, возвращённые Python-функцией
generated_all = out.get('all_repeats')
generated_summary = out.get('summary_top5_major_arc')

# Путь к директории, где в репозитории хранятся R-генерированные файлы для сравнения
r_dir = (base / '..' / 'data' / '2_derived' / 'pos8473_TtoC').resolve()
r_all = r_dir / 'my_mtDNA_repeat_all_repeats.csv'
r_summary = r_dir / 'my_mtDNA_repeat_major_arc_summary_top5.csv'

print('Generated files:', generated_all, generated_summary)
print('Reference files exist (R):', r_all.exists(), r_summary.exists())

# Сравнение табличных файлов выполняем только если оба файла существуют
if r_all.exists() and Path(generated_all).exists():
    try:
        # Читаем сгенерированный Python-файлом CSV и R-референс
        ga = pd.read_csv(generated_all)
        ra = pd.read_csv(r_all)
        print('all_repeats shapes gen/ref:', ga.shape, ra.shape)
        # Простая слияние/индикатор: строки, которые присутствуют только в одной из таблиц
        merged = ga.merge(ra.drop_duplicates(), how='outer', indicator=True)
        diff_count = (merged['_merge'] != 'both').sum()
        print('Rows that differ (symmetric):', diff_count)
    except Exception as e:
        # Если чтение или сравнение упало — печатаем ошибку, но не выходим с ошибкой
        print('Error comparing all_repeats:', e)

# Аналогичная проверка для summary таблиц
if r_summary.exists() and Path(generated_summary).exists():
    try:
        gs = pd.read_csv(generated_summary)
        rs = pd.read_csv(r_summary)
        print('summary shapes gen/ref:', gs.shape, rs.shape)
        merged = gs.merge(rs.drop_duplicates(), how='outer', indicator=True)
        diff_count = (merged['_merge'] != 'both').sum()
        print('Summary rows that differ (symmetric):', diff_count)
    except Exception as e:
        print('Error comparing summary:', e)

print('Done')
