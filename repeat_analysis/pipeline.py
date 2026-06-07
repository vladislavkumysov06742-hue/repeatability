from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from .reference import read_reference
from .clean import clean_file
from .stats import run_stats

def run_all(input_dir, output_dir, ref_path, workers, window_size=0, snv_file=None):
    cleaned_dir = Path(output_dir) / 'cleaned'
    stats_dir = Path(output_dir) / 'stats'
    ref_seq = read_reference(ref_path)
    MT_LEN = len(ref_seq)
    print("=== Этап 1: Очистка ===")
    txt_files = list(Path(input_dir).glob("*.txt"))
    if not txt_files:
        print("Нет .txt файлов в", input_dir); return
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(clean_file, f, cleaned_dir, MT_LEN): f for f in txt_files}
        for future in as_completed(futures):
            f = futures[future]
            try: future.result()
            except Exception as e: print(f"   Ошибка {f.name}: {e}")
    print("Очистка завершена.")
    print("=== Этап 2: Статистика ===")
    run_stats(cleaned_dir, stats_dir, ref_path, workers, window_size, snv_file)
    print("Полный пайплайн завершён.")

