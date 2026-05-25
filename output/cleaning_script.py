import argparse
import subprocess
import platform
import tempfile
import heapq
import os
import shutil
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
#  КОНФИГУРАЦИЯ
# ============================================================
MT_LEN = 16569
MAJOR_ARC_START = 5747
MAJOR_ARC_END = 407

PHENOTYPE = {
    "8251":  1,
    "8472":  1,
    "8473":  1,
    "12705": -1,
    "14798": -1,
    "16223": 1,
}
SNP_NAMES = ["8251", "8472", "8473", "12705", "14798", "16223"]

REF_SUFFIX = {"8251": "G", "8472": "C", "8473": "T", "12705": "C", "14798": "T", "16223": "C"}
ALT_SUFFIX = {"8251": "A", "8472": "T", "8473": "C", "12705": "T", "14798": "C", "16223": "T"}

# Пастельные цвета
COLOR_POS_BAR = "#a8e6cf"
COLOR_NEG_BAR = "#ffb3b3"
COLOR_PHEN_POS = "#89CFF0"
COLOR_PHEN_NEG = "#FADADD"

if platform.system() == 'Windows':
    HAS_SYSTEM_SORT = False
    print("Windows: будет использоваться встроенная сортировка (надёжнее).")
else:
    HAS_SYSTEM_SORT = shutil.which('sort') is not None
    if not HAS_SYSTEM_SORT:
        print("Системный sort не найден. Будет использоваться встроенная сортировка (медленнее).")
        
# ============================================================
#  ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ PLOT
# ============================================================
def interval_in_major_arc(start: int, end: int) -> bool:
    if start <= end:
        if start >= MAJOR_ARC_START:
            return end <= MT_LEN
        elif end <= MAJOR_ARC_END:
            return start >= 1
        else:
            return False
    else:
        in_upper = (start >= MAJOR_ARC_START) and (start <= MT_LEN)
        in_lower = (end >= 1) and (end <= MAJOR_ARC_END)
        return in_upper and in_lower

def get_max_perfect_repeat_details(file_path: Path):
    max_len = 0
    best_seq = ""
    best_start = 0
    best_end = 0
    with file_path.open(encoding='utf-8') as f:
        f.readline()
        for line in f:
            line = line.strip()
            if not line: continue
            parts = line.split('\t')
            if len(parts) < 10: continue
            try:
                hamming = int(parts[9])
                rstart = int(parts[7])
                rend   = int(parts[8])
            except ValueError:
                continue
            if hamming == 0 and interval_in_major_arc(rstart, rend):
                repeat_seq = parts[6]
                seq_len = len(repeat_seq)
                if seq_len > max_len:
                    max_len = seq_len
                    best_seq = repeat_seq
                    best_start = rstart
                    best_end = rend
    return max_len, best_seq, best_start, best_end

def plot_delta_r(delta_values: dict, snp_names: list, output_dir: Path):
    snps = snp_names
    deltas = [delta_values[snp] for snp in snps]
    colors = [COLOR_POS_BAR if d > 0 else COLOR_NEG_BAR if d < 0 else "#dddddd" for d in deltas]

    plt.figure(figsize=(10, 8))
    bars = plt.bar(snps, deltas, color=colors, edgecolor='black', linewidth=0.8)
    plt.axhline(0, color='black', linewidth=0.8, linestyle='--')
    plt.xlabel('SNP')
    plt.ylabel('ΔR (п.н.) = Alt_len - Ref_len')
    plt.title('Изменение повторяемости самого длинного совершенного повтора (только major arc)', fontsize=12)

    for bar, delta in zip(bars, deltas):
        height = bar.get_height()
        if delta >= 0:
            va = 'bottom'
            y_offset = 0
        else:
            va = 'top'
            y_offset = 0
        plt.text(bar.get_x() + bar.get_width()/2, height + y_offset,
                 f'{delta}', ha='center', va=va, fontsize=9, fontweight='bold')
    plt.tight_layout()
    out_file = output_dir / 'delta_R_bars.png'
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f"   Столбчатая диаграмма: {out_file}")

def plot_forest(delta_values: dict, snp_names: list, phenotype: dict, output_dir: Path):
    snps = snp_names
    deltas = [delta_values[snp] for snp in snps]
    phenos = [phenotype[snp] for snp in snps]
    y_pos = np.arange(len(snps))

    plt.figure(figsize=(6, 4))
    colors = [COLOR_PHEN_POS if p == 1 else COLOR_PHEN_NEG for p in phenos]
    plt.scatter(deltas, y_pos, c=colors, s=100, edgecolors='black', linewidth=0.8, zorder=3)
    plt.axvline(0, color='gray', linestyle='--', linewidth=0.8, alpha=0.7)
    plt.yticks(y_pos, snps)
    plt.xlabel('ΔR (п.н.)')
    plt.title('Forest plot ΔR для каждого SNP')
    x_min = min(deltas) - 2 if min(deltas) < 0 else -2
    x_max = max(deltas) + 2 if max(deltas) > 0 else 2
    plt.xlim(x_min, x_max)
    for i, (snp, delta) in enumerate(zip(snps, deltas)):
        plt.annotate(f'{delta}', (delta, i), xytext=(5, 0), textcoords='offset points', fontsize=8)
    plt.tight_layout()
    out_file = output_dir / 'forest_plot_deltaR.png'
    plt.savefig(out_file, dpi=150)
    plt.close()
    print(f"   Лесной график: {out_file}")

# ============================================================
#  ОЧИСТКА: УДАЛЕНИЕ ВЛОЖЕННЫХ ПОВТОРОВ
# ============================================================

def norm_interval(start: int, end: int) -> tuple[int, int]:
    """Преобразует crossing-интервал в линейный (start, end+MT_LEN) для сравнения."""
    if start <= end:
        return start, end
    return start, end + MT_LEN

def check_contained(m_s1, m_e1, r_s1, r_e1, m_s2, m_e2, r_s2, r_e2):
    """Возвращает True, если интервал 2 полностью вложен в интервал 1."""
    return (m_s1 <= m_s2 and m_e2 <= m_e1 and
            r_s1 <= r_s2 and r_e2 <= r_e1)

def process_group(group_lines):
    """
    Принимает список строк одной группы (одинаковые cluster_id).
    Формат строки: cluster_id\tm_start\tr_start\tm_end\tr_end\toriginal_line
    Возвращает список original_line, прошедших фильтр.
    """
    parsed = []
    for line in group_lines:
        parts = line.split('\t', 5)
        m_s = int(parts[1])
        r_s = int(parts[2])
        m_e = int(parts[3])
        r_e = int(parts[4])
        original = parts[5]
        # Определяем hamming из оригинальной строки
        orig_parts = original.split('\t')
        try:
            hamming = int(orig_parts[9])
        except (IndexError, ValueError):
            hamming = -1
        parsed.append((m_s, m_e, r_s, r_e, hamming, original))

    kept_perfect = []   # совершенные повторы (hamming=0)
    kept_imperfect = [] # несовершенные (hamming>0)

    for m_s, m_e, r_s, r_e, hamming, orig in parsed:
        if hamming == 0:
            # Совершенный повтор: проверяем, не поглощён ли другими совершенными
            absorbed = False
            for km_s, km_e, kr_s, kr_e, _, _ in kept_perfect:
                if check_contained(km_s, km_e, kr_s, kr_e, m_s, m_e, r_s, r_e):
                    absorbed = True
                    break
            if not absorbed:
                # Удаляем из kept_perfect все, которые поглощаются текущим
                new_perfect = []
                for k in kept_perfect:
                    km_s, km_e, kr_s, kr_e, _, _ = k
                    if not check_contained(m_s, m_e, r_s, r_e, km_s, km_e, kr_s, kr_e):
                        new_perfect.append(k)
                kept_perfect = new_perfect
                kept_perfect.append((m_s, m_e, r_s, r_e, hamming, orig))
        else:
            # Несовершенный повтор: проверяем, не поглощён ли он сохранёнными совершенными
            # или другими несовершенными (из kept_imperfect)
            absorbed = False
            # проверяем по совершенным
            for km_s, km_e, kr_s, kr_e, _, _ in kept_perfect:
                if check_contained(km_s, km_e, kr_s, kr_e, m_s, m_e, r_s, r_e):
                    absorbed = True
                    break
            if not absorbed:
                for km_s, km_e, kr_s, kr_e, _, _ in kept_imperfect:
                    if check_contained(km_s, km_e, kr_s, kr_e, m_s, m_e, r_s, r_e):
                        absorbed = True
                        break
            if not absorbed:
                # Удаляем из kept_imperfect все, которые поглощаются текущим
                new_imperfect = []
                for k in kept_imperfect:
                    km_s, km_e, kr_s, kr_e, _, _ = k
                    if not check_contained(m_s, m_e, r_s, r_e, km_s, km_e, kr_s, kr_e):
                        new_imperfect.append(k)
                kept_imperfect = new_imperfect
                kept_imperfect.append((m_s, m_e, r_s, r_e, hamming, orig))

    # Собираем результат: все совершенные + оставшиеся несовершенные
    result = [orig for (_, _, _, _, _, orig) in kept_perfect]
    result.extend(orig for (_, _, _, _, _, orig) in kept_imperfect)
    return result

def cluster_by_repeat_overlap(records):
    """
    Принимает список (m_s, m_e, r_s, r_e, original_line).
    Сортирует по r_s, объединяет в кластеры по ПЕРЕКРЫТИЮ интервалов повторов (без зазора).
    Возвращает список cluster_id.
    """
    indexed = list(enumerate(records))
    indexed.sort(key=lambda x: x[1][2])  # r_s

    cluster_ids = [0] * len(records)
    if not indexed:
        return cluster_ids

    current_cluster = [indexed[0][0]]
    current_max_r_e = indexed[0][1][3]
    next_cluster_id = 0
    for orig_idx, (_, _, r_s, r_e, _) in indexed[1:]:
        if r_s <= current_max_r_e:      # перекрытие
            current_cluster.append(orig_idx)
            if r_e > current_max_r_e:
                current_max_r_e = r_e
        else:
            for i in current_cluster:
                cluster_ids[i] = next_cluster_id
            next_cluster_id += 1
            current_cluster = [orig_idx]
            current_max_r_e = r_e
    for i in current_cluster:
        cluster_ids[i] = next_cluster_id
    return cluster_ids

def _sort_key(line):
    """Ключ для сортировки: cluster_id, m_start, r_start, -m_end, -r_end."""
    parts = line.split('\t', 5)
    return (int(parts[0]), int(parts[1]), int(parts[2]), -int(parts[3]), -int(parts[4]))


def external_sort(input_path, output_path):
    """
    Сортирует временный файл.
    Ключи: cluster_id asc, m_start asc, r_start asc, m_end desc, r_end desc.
    Использует системный sort, если доступен, иначе fallback на Python.
    """
    if HAS_SYSTEM_SORT:
        cmd = [
            'sort', '-t', '\t',
            '-k1,1n',              # cluster_id asc
            '-k2,2n',              # m_start asc
            '-k3,3n',              # r_start asc
            '-k4,4rn',             # m_end desc
            '-k5,5rn',             # r_end desc
            input_path
        ]
        with open(output_path, 'w', encoding='utf-8') as out:
            subprocess.run(cmd, stdout=out, check=True)
        return

    # Fallback: чанковая сортировка и слияние
    chunk_size = 500_000
    chunk_files = []

    with open(input_path, 'r', encoding='utf-8') as fin:
        chunk = []
        for line in fin:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                chunk.sort(key=_sort_key)
                tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8', suffix='.chunk')
                tmp.writelines(chunk)
                tmp.close()
                chunk_files.append(tmp.name)
                chunk = []
        if chunk:
            chunk.sort(key=_sort_key)
            tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8', suffix='.chunk')
            tmp.writelines(chunk)
            tmp.close()
            chunk_files.append(tmp.name)

    class LineWrapper:
        __slots__ = ('line', 'key')
        def __init__(self, line):
            self.line = line
            self.key = _sort_key(line)
        def __lt__(self, other):
            return self.key < other.key

    file_handles = []
    iterators = []
    try:
        for fname in chunk_files:
            f = open(fname, 'r', encoding='utf-8')
            file_handles.append(f)
            iterators.append((LineWrapper(line) for line in f))
        with open(output_path, 'w', encoding='utf-8') as fout:
            for wrapped in heapq.merge(*iterators):
                fout.write(wrapped.line)
    finally:
        for f in file_handles:
            f.close()
        for fname in chunk_files:
            try:
                os.unlink(fname)
            except OSError:
                pass

def clean_file(file_path: Path, output_dir: Path) -> None:
    print(f"   Обработка {file_path.name} ...")
    content = file_path.read_text(encoding='utf-8').splitlines()
    if len(content) < 2:
        return
    header = content[0]
    lines = content[1:]

    # Первичный фильтр
    valid = []
    for line in lines:
        line = line.strip()
        if not line: continue
        parts = line.split('\t')
        if len(parts) < 10: continue
        try:
            motif_length = int(parts[3])
            motif_start = int(parts[4])
            motif_end = int(parts[5])
            repeat_start = int(parts[7])
            repeat_end = int(parts[8])
        except ValueError:
            continue
        if motif_length < 5: continue
        if motif_start == repeat_start and motif_end == repeat_end: continue
        valid.append(line)

    if not valid:
        output_dir.mkdir(parents=True, exist_ok=True)
        out_path = output_dir / file_path.name
        out_path.write_text(header + '\n', encoding='utf-8')
        return

    # Нормализация координат
    records = []  # (m_s, m_e, r_s, r_e, original_line)
    for line in valid:
        parts = line.split('\t')
        m_s, m_e = norm_interval(int(parts[4]), int(parts[5]))
        r_s, r_e = norm_interval(int(parts[7]), int(parts[8]))
        records.append((m_s, m_e, r_s, r_e, line))

    # Кластеризация по перекрытию повторов (gap=0)
    cluster_ids = cluster_by_repeat_overlap(records)

    # Временный файл: cluster_id, m_start, r_start, m_end, r_end, original_line
    tmp_input = tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8', suffix='.presort')
    for cid, (m_s, m_e, r_s, r_e, line) in zip(cluster_ids, records):
        tmp_input.write(f"{cid}\t{m_s}\t{r_s}\t{m_e}\t{r_e}\t{line}\n")
    tmp_input.close()

    # Сортировка временного файла
    tmp_sorted = tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8', suffix='.sorted')
    tmp_sorted.close()
    external_sort(tmp_input.name, tmp_sorted.name)

    kept_overall = []
    with open(tmp_sorted.name, 'r', encoding='utf-8') as f:
        current_cid = None
        group = []
        for line in f:
            line = line.rstrip('\n')
            parts = line.split('\t', 5)
            cid = int(parts[0])
            if cid != current_cid:
                if group:
                    kept_overall.extend(process_group(group))
                current_cid = cid
                group = [line]
            else:
                group.append(line)
        if group:
            kept_overall.extend(process_group(group))

    # Финальная сортировка по эффективной длине (motif.length - hamming)
    def effective_length(line):
        parts = line.split('\t')
        try:
            return int(parts[3]) - int(parts[9])
        except (IndexError, ValueError):
            return 0
    kept_overall.sort(key=effective_length, reverse=True)

    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / file_path.name
    with out_path.open('w', encoding='utf-8') as f:
        f.write(header + '\n')
        for line in kept_overall:
            f.write(line + '\n')

    os.unlink(tmp_input.name)
    os.unlink(tmp_sorted.name)
    print(f"   Готово: {out_path} (оставлено {len(kept_overall)} повторов)")

# ============================================================
#  РЕЖИМ PLOT
# ============================================================
def run_plot(input_dir: Path, output_dir: Path):
    print(f"Режим plot: чтение файлов из {input_dir}")
    delta_values = {}
    details = {}
    missing = []

    for snp in SNP_NAMES:
        ref_path = input_dir / f"01KP.{snp}.{REF_SUFFIX[snp]}.txt"
        alt_path = input_dir / f"01KP.{snp}.{ALT_SUFFIX[snp]}.txt"
        if not ref_path.exists():
            missing.append(str(ref_path))
            continue
        if not alt_path.exists():
            missing.append(str(alt_path))
            continue
        ref_len, ref_seq, ref_start, ref_end = get_max_perfect_repeat_details(ref_path)
        alt_len, alt_seq, alt_start, alt_end = get_max_perfect_repeat_details(alt_path)
        delta = alt_len - ref_len
        delta_values[snp] = delta
        details[snp] = {
            'ref': (ref_len, ref_seq, ref_start, ref_end),
            'alt': (alt_len, alt_seq, alt_start, alt_end)
        }
        print(f"   {snp}: Ref_len={ref_len}, Alt_len={alt_len}, ΔR={delta}")

    if missing:
        print("Предупреждение: не найдены файлы:", missing)

    output_dir.mkdir(parents=True, exist_ok=True)

    log_path = output_dir / "used_repeats_log.txt"
    with open(log_path, 'w', encoding='utf-8') as log:
        log.write("SNP\tAllele\tLength\tRepeat_Sequence\tStart\tEnd\n")
        for snp in SNP_NAMES:
            if snp not in details:
                continue
            r_len, r_seq, r_st, r_en = details[snp]['ref']
            a_len, a_seq, a_st, a_en = details[snp]['alt']
            log.write(f"{snp}\tRef\t{r_len}\t{r_seq}\t{r_st}\t{r_en}\n")
            log.write(f"{snp}\tAlt\t{a_len}\t{a_seq}\t{a_st}\t{a_en}\n")
    print(f"   Лог сохранён: {log_path}")

    plot_delta_r(delta_values, SNP_NAMES, output_dir)
    plot_forest(delta_values, SNP_NAMES, PHENOTYPE, output_dir)

# ============================================================
#  ОСНОВНОЙ БЛОК
# ============================================================
def main():
    parser = argparse.ArgumentParser(description='Обработка TSV-файлов с повторами')
    parser.add_argument('mode', choices=['clean', 'plot'], help='Режим: clean или plot')
    parser.add_argument('--input_dir', type=str, help='Папка с исходными .txt файлами')
    parser.add_argument('--output_dir', type=str, help='Для clean: папка для очищенных; для plot: папка для графиков')
    parser.add_argument('--workers', type=int, default=8, help='Количество потоков/процессов (для clean)')
    args = parser.parse_args()

    if args.mode == 'clean':
        if args.input_dir is None:
            input_dir = Path(r"D:\repeatability\output\01KP")
        else:
            input_dir = Path(args.input_dir)
        if args.output_dir is None:
            output_dir = input_dir / "cleaned"
        else:
            output_dir = Path(args.output_dir)

        print(f"Режим clean: {input_dir}/*.txt -> {output_dir}")
        txt_files = list(input_dir.glob("*.txt"))
        if not txt_files:
            print("Нет .txt файлов.")
            return

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(clean_file, f, output_dir): f for f in txt_files}
            for future in as_completed(futures):
                f = futures[future]
                try:
                    future.result()
                except Exception as e:
                    print(f"   Ошибка при обработке {f.name}: {e}")
        print("Все файлы очищены.")

    elif args.mode == 'plot':
        if args.input_dir is None:
            input_dir = Path(r"D:\repeatability\mut")
        else:
            input_dir = Path(args.input_dir)
        if args.output_dir is None:
            output_dir = input_dir / "plots"
        else:
            output_dir = Path(args.output_dir)

        print(f"Режим plot: файлы = {input_dir}, графики = {output_dir}")
        run_plot(input_dir, output_dir)
        print("Графики построены.")

if __name__ == "__main__":
    main()