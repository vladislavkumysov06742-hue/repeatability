#!/usr/bin/env python3
"""
Полный пайплайн анализа повторяемости мтДНК.
Режимы: clean, plot, stats.
"""

import argparse, subprocess, tempfile, heapq, os, shutil, platform
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde, nbinom, poisson, norm

# ============================================================
#  КОНФИГУРАЦИЯ
# ============================================================
MT_LEN = 16569
MAJOR_ARC_START = 5747
MAJOR_ARC_END = 407

PHENOTYPE = {
    "8251":  1, "8472":  1, "8473":  1,
    "12705": -1, "14798": -1, "16223": 1,
}
SNP_NAMES = ["8251", "8472", "8473", "12705", "14798", "16223"]

REF_SUFFIX = {"8251": "G", "8472": "C", "8473": "T", "12705": "C", "14798": "T", "16223": "C"}
ALT_SUFFIX = {"8251": "A", "8472": "T", "8473": "C", "12705": "T", "14798": "C", "16223": "T"}

COLOR_POS_BAR = "#a8e6cf"
COLOR_NEG_BAR = "#ffb3b3"
COLOR_PHEN_POS = "#89CFF0"
COLOR_PHEN_NEG = "#FADADD"

# Глобальная проверка sort
if platform.system() == 'Windows':
    HAS_SYSTEM_SORT = False
else:
    HAS_SYSTEM_SORT = shutil.which('sort') is not None

# ============================================================
#  ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ============================================================
def norm_interval(start, end):
    return (start, end) if start <= end else (start, end + MT_LEN)

def check_contained(m_s1, m_e1, r_s1, r_e1, m_s2, m_e2, r_s2, r_e2):
    return (m_s1 <= m_s2 and m_e2 <= m_e1 and r_s1 <= r_s2 and r_e2 <= r_e1)

def process_group(group_lines):
    parsed = []
    for line in group_lines:
        parts = line.split('\t', 5)
        m_s, r_s = int(parts[1]), int(parts[2])
        m_e, r_e = int(parts[3]), int(parts[4])
        original = parts[5]
        orig_parts = original.split('\t')
        try:
            hamming = int(orig_parts[9])
        except:
            hamming = -1
        parsed.append((m_s, m_e, r_s, r_e, hamming, original))
    kept_perfect, kept_imperfect = [], []
    for m_s, m_e, r_s, r_e, hamming, orig in parsed:
        if hamming == 0:
            absorbed = any(check_contained(k[0],k[1],k[2],k[3], m_s,m_e,r_s,r_e) for k in kept_perfect)
            if not absorbed:
                kept_perfect = [k for k in kept_perfect if not check_contained(m_s,m_e,r_s,r_e, k[0],k[1],k[2],k[3])]
                kept_perfect.append((m_s,m_e,r_s,r_e, hamming, orig))
        else:
            absorbed = any(check_contained(k[0],k[1],k[2],k[3], m_s,m_e,r_s,r_e) for k in kept_perfect + kept_imperfect)
            if not absorbed:
                kept_imperfect = [k for k in kept_imperfect if not check_contained(m_s,m_e,r_s,r_e, k[0],k[1],k[2],k[3])]
                kept_imperfect.append((m_s,m_e,r_s,r_e, hamming, orig))
    return [orig for (_,_,_,_,_,orig) in kept_perfect] + [orig for (_,_,_,_,_,orig) in kept_imperfect]

def cluster_by_repeat_overlap(records):
    indexed = list(enumerate(records))
    indexed.sort(key=lambda x: x[1][2])
    if not indexed:
        return []
    cluster_ids = [0] * len(records)
    current_cluster = [indexed[0][0]]
    current_max_r_e = indexed[0][1][3]
    next_id = 0
    for orig_idx, (_, _, r_s, r_e, _) in indexed[1:]:
        if r_s <= current_max_r_e:
            current_cluster.append(orig_idx)
            if r_e > current_max_r_e:
                current_max_r_e = r_e
        else:
            for i in current_cluster:
                cluster_ids[i] = next_id
            next_id += 1
            current_cluster = [orig_idx]
            current_max_r_e = r_e
    for i in current_cluster:
        cluster_ids[i] = next_id
    return cluster_ids

def _sort_key(line):
    parts = line.split('\t', 5)
    return (int(parts[0]), int(parts[1]), int(parts[2]), -int(parts[3]), -int(parts[4]))

def external_sort(input_path, output_path):
    if HAS_SYSTEM_SORT:
        cmd = ['sort', '-t', '\t', '-k1,1n', '-k2,2n', '-k3,3n', '-k4,4rn', '-k5,5rn', input_path]
        subprocess.run(cmd, stdout=open(output_path, 'w'), check=True)
        return
    chunk_size = 500_000
    chunk_files = []
    with open(input_path, 'r') as fin:
        chunk = []
        for line in fin:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                chunk.sort(key=_sort_key)
                tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.chunk')
                tmp.writelines(chunk); tmp.close()
                chunk_files.append(tmp.name)
                chunk = []
        if chunk:
            chunk.sort(key=_sort_key)
            tmp = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.chunk')
            tmp.writelines(chunk); tmp.close()
            chunk_files.append(tmp.name)
    class LineWrapper:
        __slots__ = ('line', 'key')
        def __init__(self, line): self.line = line; self.key = _sort_key(line)
        def __lt__(self, o): return self.key < o.key
    file_handles = [open(f, 'r') for f in chunk_files]
    iterators = [(LineWrapper(l) for l in f) for f in file_handles]
    with open(output_path, 'w') as fout:
        for w in heapq.merge(*iterators):
            fout.write(w.line)
    for f in file_handles: f.close()
    for f in chunk_files: os.unlink(f)

def gc_content(seq):
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) if seq else 0.0

def h_bonds(seq):
    bonds = {'A':2, 'T':2, 'G':3, 'C':3}
    return sum(bonds.get(b, 0) for b in seq.upper())

# ============================================================
#  CLEAN
# ============================================================
def clean_file(file_path, output_dir):
    print(f"   Обработка {file_path.name} ...")
    lines = file_path.read_text().splitlines()
    if len(lines) < 2: return
    header, lines = lines[0], lines[1:]
    valid = []
    for l in lines:
        if not l.strip(): continue
        p = l.split('\t')
        if len(p) < 10: continue
        try:
            ml, ms, me, rs, re = int(p[3]), int(p[4]), int(p[5]), int(p[7]), int(p[8])
        except: continue
        if ml < 5: continue
        if ms == rs and me == re: continue
        valid.append(l)
    if not valid:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / file_path.name).write_text(header + '\tmotif_GC\trepeat_GC\tmotif_H_bonds\trepeat_H_bonds\n')
        return
    records = [(norm_interval(int(p[4]),int(p[5]))[0], norm_interval(int(p[4]),int(p[5]))[1],
                norm_interval(int(p[7]),int(p[8]))[0], norm_interval(int(p[7]),int(p[8]))[1], l)
               for l in valid for p in [l.split('\t')]]
    cids = cluster_by_repeat_overlap(records)
    tmp_in = tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.presort')
    for cid, (ms, me, rs, re, line) in zip(cids, records):
        tmp_in.write(f"{cid}\t{ms}\t{rs}\t{me}\t{re}\t{line}\n")
    tmp_in.close()
    tmp_sorted = tempfile.NamedTemporaryFile(delete=False, suffix='.sorted')
    tmp_sorted.close()
    external_sort(tmp_in.name, tmp_sorted.name)
    kept = []
    with open(tmp_sorted.name) as f:
        cur_cid, group = None, []
        for line in f:
            line = line.rstrip()
            parts = line.split('\t', 5)
            cid = int(parts[0])
            if cid != cur_cid:
                if group: kept.extend(process_group(group))
                cur_cid = cid; group = [line]
            else: group.append(line)
        if group: kept.extend(process_group(group))
    kept.sort(key=lambda l: int(l.split('\t')[3]) - int(l.split('\t')[9]), reverse=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / file_path.name
    with out.open('w') as f:
        f.write(header + '\tmotif_GC\trepeat_GC\tmotif_H_bonds\trepeat_H_bonds\n')
        for l in kept:
            p = l.split('\t')
            mgc = round(gc_content(p[2]), 4); rgc = round(gc_content(p[6]), 4)
            mhb = h_bonds(p[2]); rhb = h_bonds(p[6])
            f.write('\t'.join(p + [str(mgc), str(rgc), str(mhb), str(rhb)]) + '\n')
    os.unlink(tmp_in.name); os.unlink(tmp_sorted.name)
    print(f"   Готово: {out} (оставлено {len(kept)} повторов)")

# ============================================================
#  PLOT (без изменений)
# ============================================================
def interval_in_major_arc(start, end):
    if start <= end:
        if start >= MAJOR_ARC_START: return end <= MT_LEN
        elif end <= MAJOR_ARC_END: return start >= 1
        return False
    in_upper = start >= MAJOR_ARC_START and start <= MT_LEN
    in_lower = end >= 1 and end <= MAJOR_ARC_END
    return in_upper and in_lower

def get_max_perfect_repeat_details(file_path):
    max_len = 0; best_seq = ""; best_start = 0; best_end = 0
    with file_path.open() as f:
        f.readline()
        for line in f:
            if not line.strip(): continue
            p = line.strip().split('\t')
            if len(p) < 10: continue
            try: hamming = int(p[9]); rstart = int(p[7]); rend = int(p[8])
            except: continue
            if hamming == 0 and interval_in_major_arc(rstart, rend):
                seq_len = len(p[6])
                if seq_len > max_len:
                    max_len = seq_len; best_seq = p[6]; best_start = rstart; best_end = rend
    return max_len, best_seq, best_start, best_end

def plot_delta_r(delta_values, snp_names, output_dir):
    deltas = [delta_values[s] for s in snp_names]
    colors = [COLOR_POS_BAR if d>0 else COLOR_NEG_BAR if d<0 else "#dddddd" for d in deltas]
    fig, ax = plt.subplots(figsize=(10,8))
    bars = ax.bar(snp_names, deltas, color=colors, edgecolor='black')
    ax.axhline(0, color='black', linestyle='--')
    ax.set_xlabel('SNP'); ax.set_ylabel('ΔR (п.н.)')
    ax.set_title('Изменение самого длинного совершенного повтора (major arc)')
    for bar, d in zip(bars, deltas):
        va = 'bottom' if d>=0 else 'top'
        ax.text(bar.get_x()+bar.get_width()/2, bar.get_height(), f'{d}', ha='center', va=va, fontweight='bold')
    fig.tight_layout()
    fig.savefig(output_dir / 'delta_R_bars.png', dpi=150); plt.close()

def plot_forest(delta_values, snp_names, phenotype, output_dir):
    deltas = [delta_values[s] for s in snp_names]
    phenos = [phenotype[s] for s in snp_names]
    fig, ax = plt.subplots(figsize=(6,4))
    colors = [COLOR_PHEN_POS if p==1 else COLOR_PHEN_NEG for p in phenos]
    ax.scatter(deltas, range(len(deltas)), c=colors, s=100, edgecolors='black', zorder=3)
    ax.axvline(0, color='gray', linestyle='--')
    ax.set_yticks(range(len(deltas))); ax.set_yticklabels(snp_names)
    ax.set_xlabel('ΔR (п.н.)'); ax.set_title('Forest plot ΔR')
    for i, d in enumerate(deltas):
        ax.annotate(f'{d}', (d,i), xytext=(5,0), textcoords='offset points', fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / 'forest_plot_deltaR.png', dpi=150); plt.close()

def run_plot(input_dir, output_dir):
    print(f"Режим plot: {input_dir}")
    delta_values, details, missing = {}, {}, []
    for snp in SNP_NAMES:
        ref_path = input_dir / f"01KP.{snp}.{REF_SUFFIX[snp]}.txt"
        alt_path = input_dir / f"01KP.{snp}.{ALT_SUFFIX[snp]}.txt"
        if not ref_path.exists(): missing.append(str(ref_path))
        elif not alt_path.exists(): missing.append(str(alt_path))
        else:
            rlen, rseq, rs, re = get_max_perfect_repeat_details(ref_path)
            alen, aseq, as_, ae = get_max_perfect_repeat_details(alt_path)
            delta = alen - rlen
            delta_values[snp] = delta
            details[snp] = {'ref':(rlen, rseq, rs, re), 'alt':(alen, aseq, as_, ae)}
            print(f"   {snp}: Ref_len={rlen}, Alt_len={alen}, ΔR={delta}")
    if missing: print("Предупреждение: не найдены файлы:", missing)
    if not delta_values: print("Нет данных для графиков."); return
    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "used_repeats_log.txt", 'w') as log:
        log.write("SNP\tAllele\tLength\tRepeat_Sequence\tStart\tEnd\n")
        for s in SNP_NAMES:
            if s not in details: continue
            r, a = details[s]['ref'], details[s]['alt']
            log.write(f"{s}\tRef\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n")
            log.write(f"{s}\tAlt\t{a[0]}\t{a[1]}\t{a[2]}\t{a[3]}\n")
    plot_delta_r(delta_values, SNP_NAMES, output_dir)
    plot_forest(delta_values, SNP_NAMES, PHENOTYPE, output_dir)

# ============================================================
#  STATS
# ============================================================
def read_reference(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'): continue
            seq.append(line.strip().upper())
    return ''.join(seq)

def load_cleaned_files(input_dir, reference_seq, log_print):
    per_pos = defaultdict(dict)
    all_repeats_raw, all_repeats_weighted = [], []
    ref_nuc = {}
    metrics = defaultdict(dict)
    n_files = 0
    for fpath in input_dir.glob("01KP.*.*.txt"):
        parts = fpath.stem.split('.')
        if len(parts) != 3: continue
        try: pos = int(parts[1])
        except: continue
        nuc = parts[2].upper()
        if nuc not in 'ACGT': continue
        total = perfect = max_eff = sum_eff = 0
        lines = []
        with fpath.open() as f:
            f.readline()
            for line in f:
                line = line.strip()
                if not line: continue
                p = line.split('\t')
                if len(p) < 11: continue
                lines.append(line)
                total += 1
                try: hamming = int(p[9]); motif_len = int(p[3])
                except: continue
                eff = motif_len - hamming
                if hamming == 0: perfect += 1
                if eff > max_eff: max_eff = eff
                sum_eff += eff
                try: rs, re = int(p[7]), int(p[8])
                except: continue
                if rs <= re:
                    all_repeats_raw.append((rs, re))
                    all_repeats_weighted.append((rs, re, eff))
                else:
                    all_repeats_raw.append((rs, MT_LEN)); all_repeats_raw.append((1, re))
                    all_repeats_weighted.append((rs, MT_LEN, eff)); all_repeats_weighted.append((1, re, eff))
        per_pos[pos][nuc] = lines
        metrics[pos][nuc] = {'total': total, 'perfect': perfect, 'max_eff': max_eff, 'sum_eff': sum_eff}
        if 1 <= pos <= len(reference_seq):
            ref_nuc[pos] = reference_seq[pos-1]
        n_files += 1
    log_print(f"Загружено файлов: {n_files}")
    log_print(f"Всего повторов (raw): {len(all_repeats_raw)}")
    return per_pos, ref_nuc, all_repeats_raw, all_repeats_weighted, metrics

def parametric_pvalue(observed, background, metric):
    """
    Автоматический выбор распределения и расчёт двустороннего p-value.
    """
    if len(background) < 5:
        return 1.0
    bg = np.array(background, dtype=float)
    mu = np.mean(bg)
    sigma = np.std(bg)
    if metric in ('total', 'perfect'):
        # Счётные данные: выбор между Poisson и Negative Binomial
        var = np.var(bg)
        if var > mu + 1:  # overdispersion -> Negative Binomial
            # Параметры NB: r, p
            # r = mu^2 / (var - mu), p = mu / var
            if var == mu:
                var = mu + 0.1  # избежать деления на ноль
            p_param = mu / var
            r_param = mu * p_param / (1 - p_param) if p_param < 1 else 100
            # Двусторонний тест
            if observed <= mu:
                left = nbinom.cdf(observed, r_param, p_param)
                right_obs = int(2*mu - observed)
                if right_obs >= 0:
                    right = nbinom.sf(right_obs - 1, r_param, p_param)
                else:
                    right = 1.0
                pval = left + right
            else:
                right = nbinom.sf(observed - 1, r_param, p_param)
                left_obs = int(2*mu - observed)
                if left_obs >= 0:
                    left = nbinom.cdf(left_obs, r_param, p_param)
                else:
                    left = 0.0
                pval = left + right
            return min(max(pval, 0), 1.0)
        else:
            # Poisson
            if mu < 1e-6:
                mu = 0.1
            if observed <= mu:
                left = poisson.cdf(observed, mu)
                right_obs = int(2*mu - observed)
                right = poisson.sf(right_obs - 1, mu) if right_obs >= 0 else 1.0
                pval = left + right
            else:
                right = poisson.sf(observed - 1, mu)
                left_obs = int(2*mu - observed)
                left = poisson.cdf(left_obs, mu) if left_obs >= 0 else 0.0
                pval = left + right
            return min(max(pval, 0), 1.0)
    else:
        # Непрерывные метрики (sum_eff, max_eff) – нормальное приближение
        if sigma == 0:
            return 1.0 if observed == mu else 0.0
        z = (observed - mu) / sigma
        pval = 2 * norm.sf(abs(z))
        return pval

def add_significance_lines(ax, pvals, log_print=None):
    n = len(pvals)
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='Nominal α=0.05')
    sorted_p = np.sort(pvals)
    sig = None
    for i, p in enumerate(sorted_p):
        if p <= (i+1)/n*0.05:
            sig = p
    if sig is not None:
        ax.axhline(-np.log10(sig), color='red', linestyle='--', alpha=0.5, label='FDR threshold')
        if log_print: log_print(f"   FDR порог: p≤{sig:.2e}")
    else:
        if log_print: log_print("   FDR: значимых точек нет")
    ax.legend()

def compute_coverage(all_repeats, weighted=False):
    cov = np.zeros(MT_LEN+1, dtype=float)
    for item in all_repeats:
        if weighted: s, e, w = item
        else: s, e = item; w = 1
        if e > MT_LEN: e = MT_LEN
        cov[s:e+1] += w
    return cov

def run_stats(input_dir, output_dir, ref_path, workers):
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "stats_run.log"
    logf = open(log_path, 'w', encoding='utf-8')
    def log_print(msg):
        print(msg)
        print(msg, file=logf)

    log_print("Загрузка референсной последовательности...")
    ref_seq = read_reference(ref_path)
    if len(ref_seq) != MT_LEN:
        log_print(f"Внимание: длина {len(ref_seq)} вместо {MT_LEN}.")
    log_print("Чтение очищенных файлов...")
    per_pos, ref_nuc, all_repeats_raw, all_repeats_weighted, metrics = load_cleaned_files(input_dir, ref_seq, log_print)
    log_print(f"Позиций с данными: {len(per_pos)}")

    # ----- фон -----
    log_print("Подготовка фонов...")
    ref_total, ref_perfect, ref_max_eff, ref_sum_eff = [], [], [], []
    for pos, nucs in per_pos.items():
        ref = ref_nuc.get(pos)
        if ref and ref in nucs:
            m = metrics[pos][ref]
            ref_total.append(m['total']); ref_perfect.append(m['perfect'])
            ref_max_eff.append(m['max_eff']); ref_sum_eff.append(m['sum_eff'])
    log_print(f" Референсных файлов в фоне: {len(ref_total)}")

    # ----- Манхеттены с фильтрацией по MAD и параметрическим p-value -----
    log_print("Построение манхеттенов...")
    metrics_names = ['total','perfect','max_eff','sum_eff']
    backgrounds = [ref_total, ref_perfect, ref_max_eff, ref_sum_eff]
    titles = ['Total repeats', 'Perfect repeats', 'Max effective length', 'Sum effective length']
    fig, axes = plt.subplots(4, 1, figsize=(14, 20))
    for idx, (metric, bg, title) in enumerate(zip(metrics_names, backgrounds, titles)):
        ax = axes[idx]
        # Собираем абсолютные Δ для порога MAD
        all_deltas = []
        for pos in sorted(per_pos.keys()):
            ref = ref_nuc.get(pos)
            if not ref or ref not in per_pos[pos]: continue
            ref_val = metrics[pos][ref][metric]
            for nuc in 'ACGT':
                if nuc == ref or nuc not in per_pos[pos]: continue
                all_deltas.append(metrics[pos][nuc][metric] - ref_val)
        if len(all_deltas) > 1:
            mad = np.median(np.abs(np.array(all_deltas) - np.median(all_deltas)))
            threshold = 1.5 * mad
        else:
            threshold = 0
        log_print(f" {title}: MAD={mad:.2f}, порог |Δ| > {threshold:.2f}")

        positions, pvals = [], []
        tested_total, tested_passed = 0, 0
        for pos in sorted(per_pos.keys()):
            ref = ref_nuc.get(pos)
            if not ref or ref not in per_pos[pos]: continue
            ref_val = metrics[pos][ref][metric]
            for nuc in 'ACGT':
                if nuc == ref or nuc not in per_pos[pos]: continue
                tested_total += 1
                alt_val = metrics[pos][nuc][metric]
                delta = alt_val - ref_val
                if abs(delta) <= threshold: continue
                tested_passed += 1
                # Параметрическое p-value
                p = parametric_pvalue(alt_val, bg, metric)
                positions.append(pos)
                pvals.append(p)
        log_print(f"   Всего тестов: {tested_total}, после фильтрации: {tested_passed}")
        if not pvals:
            ax.text(0.5,0.5,'No points', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{title}')
            continue
        logp = -np.log10(np.maximum(pvals, 1e-300))
        ax.scatter(positions, logp, c='black', s=5)
        ax.set_title(f'Manhattan plot – {title}')
        add_significance_lines(ax, pvals, log_print=log_print)
    plt.tight_layout()
    plt.savefig(output_dir / 'manhattan_metrics.png', dpi=150); plt.close()
    log_print(" Манхеттены метрик сохранены.")

    # ----- Покрытия (эмпирический p-value оставим, но для coverage его достаточно) -----
    log_print("Расчёт покрытий...")
    cov_raw = compute_coverage(all_repeats_raw, weighted=False)
    cov_weighted = compute_coverage(all_repeats_weighted, weighted=True)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 10))
    for ax, cov, title in [(ax1, cov_raw, 'Raw coverage'), (ax2, cov_weighted, 'Weighted coverage')]:
        bg = cov[1:]
        pvals_cov = []
        for pos in range(1, MT_LEN+1):
            # Для coverage используем эмпирический двусторонний p-value
            if len(bg) == 0:
                pvals_cov.append(1.0)
            else:
                med = np.median(bg)
                dev_obs = abs(cov[pos] - med)
                dev_bg = np.abs(bg - med)
                p = (np.sum(dev_bg >= dev_obs) + 1) / (len(bg) + 1)
                pvals_cov.append(min(p, 1.0))
        logp_cov = -np.log10(np.maximum(pvals_cov, 1e-300))
        ax.scatter(range(1, MT_LEN+1), logp_cov, c='steelblue', s=2)
        ax.set_title(f'Manhattan – {title}')
        add_significance_lines(ax, pvals_cov, log_print=log_print)
    plt.tight_layout()
    plt.savefig(output_dir / 'manhattan_coverage.png', dpi=150); plt.close()

    # Сглаженный coverage
    fig, ax = plt.subplots(figsize=(14,4))
    ax.plot(range(1,MT_LEN+1), cov_weighted[1:], linewidth=0.3, color='darkblue', alpha=0.5, label='Raw')
    smooth = np.convolve(cov_weighted[1:], np.ones(500)/500, mode='same')
    ax.plot(range(1,MT_LEN+1), smooth, linewidth=1.5, color='red', label='Smoothed')
    ax.legend(); ax.set_title('Weighted coverage along mtDNA')
    plt.tight_layout(); plt.savefig(output_dir / 'coverage_plot.png', dpi=150); plt.close()

    # ----- Гистограммы -----
    log_print("Гистограммы...")
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,5))
    vals = cov_weighted[1:]
    ax1.hist(vals, bins=50, density=True, color='lightblue', edgecolor='black', alpha=0.7)
    if len(vals)>1:
        kde = gaussian_kde(vals)
        xr = np.linspace(vals.min(), vals.max(), 500)
        ax1.plot(xr, kde(xr), 'r-')
    ax1.set_title('Weighted coverage')
    ax1.text(0.95,0.95, f'Mean={np.mean(vals):.1f}\nMedian={np.median(vals):.1f}\nSkew={stats.skew(vals):.2f}\nKurt={stats.kurtosis(vals):.2f}',
             transform=ax1.transAxes, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.8))
    all_cnt = [len(v) for pos_data in per_pos.values() for v in pos_data.values()]
    ax2.hist(all_cnt, bins=50, density=True, color='lightgreen', edgecolor='black', alpha=0.7)
    if len(all_cnt)>1:
        kde = gaussian_kde(all_cnt)
        xr = np.linspace(min(all_cnt), max(all_cnt), 500)
        ax2.plot(xr, kde(xr), 'r-')
    ax2.set_title('Repeats per file')
    ax2.text(0.95,0.95, f'Mean={np.mean(all_cnt):.1f}\nMedian={np.median(all_cnt):.1f}\nSkew={stats.skew(all_cnt):.2f}\nKurt={stats.kurtosis(all_cnt):.2f}',
             transform=ax2.transAxes, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.8))
    plt.tight_layout(); plt.savefig(output_dir / 'histograms.png', dpi=150); plt.close()

    # ----- Scatter -----
    log_print("Scatter plots...")
    def max_eff_len(lines, perfect_only=False):
        mx = 0
        for l in lines:
            p = l.split('\t')
            try: h = int(p[9]); ml = int(p[3])
            except: continue
            if perfect_only and h != 0: continue
            eff = ml - h
            if eff > mx: mx = eff
        return mx
    scatter_data = []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref or ref not in per_pos[pos]: continue
        ref_lines = per_pos[pos][ref]
        rl_perf = max_eff_len(ref_lines, True)
        rl_imperf = max_eff_len(ref_lines, False)
        for nuc in 'ACGT':
            if nuc == ref or nuc not in per_pos[pos]: continue
            alt_lines = per_pos[pos][nuc]
            al_perf = max_eff_len(alt_lines, True)
            al_imperf = max_eff_len(alt_lines, False)
            scatter_data.append((rl_perf, al_perf - rl_perf, rl_imperf, al_imperf - rl_imperf, al_perf, al_imperf))
    log_print(f" Точек для scatter: {len(scatter_data)}")

    def scatter_with_regression(ax, x, y, color, alpha=0.3):
        if len(x)==0: return
        x_jit = x + np.random.normal(0, 0.5, size=len(x))
        y_jit = y + np.random.normal(0, 0.5, size=len(y))
        ax.scatter(x_jit, y_jit, color=color, alpha=alpha, s=8)
        if len(x)>1:
            slope, intercept, r, p, _ = stats.linregress(x, y)
            line_x = np.array([x.min(), x.max()])
            ax.plot(line_x, slope*line_x+intercept, 'k--')
            ax.text(0.02,0.98, f'R²={r**2:.3f}\np={p:.2e}', transform=ax.transAxes, va='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=7)

    # delta perfect
    x_p = np.array([d[0] for d in scatter_data])
    y_p = np.array([d[1] for d in scatter_data])
    fig, ax = plt.subplots(figsize=(8,6))
    scatter_with_regression(ax, x_p, y_p, 'blue')
    ax.set_xlabel('Ref max eff. perfect'); ax.set_ylabel('Δ perfect')
    ax.axhline(0, color='gray', linestyle='--')
    ax.set_title('Perfect: Δ vs ref')
    plt.tight_layout(); plt.savefig(output_dir / 'scatter_delta_perfect.png', dpi=150); plt.close()

    # delta imperfect
    x_i = np.array([d[2] for d in scatter_data])
    y_i = np.array([d[3] for d in scatter_data])
    fig, ax = plt.subplots(figsize=(8,6))
    scatter_with_regression(ax, x_i, y_i, 'red')
    ax.set_xlabel('Ref max eff. imperfect'); ax.set_ylabel('Δ imperfect')
    ax.axhline(0, color='gray', linestyle='--')
    ax.set_title('Imperfect: Δ vs ref')
    plt.tight_layout(); plt.savefig(output_dir / 'scatter_delta_imperfect.png', dpi=150); plt.close()

    # alt vs ref perfect
    x_ref = np.array([d[0] for d in scatter_data])
    y_alt = np.array([d[4] for d in scatter_data])
    fig, ax = plt.subplots(figsize=(8,6))
    scatter_with_regression(ax, x_ref, y_alt, 'green')
    ax.plot([0, max(x_ref)], [0, max(x_ref)], 'k--', lw=0.5)
    ax.set_xlabel('Ref max eff. perfect'); ax.set_ylabel('Alt max eff. perfect')
    ax.set_title('Perfect: alt vs ref')
    plt.tight_layout(); plt.savefig(output_dir / 'scatter_alt_vs_ref_perfect.png', dpi=150); plt.close()

    # ----- Боксплоты 99-го перцентиля -----
    log_print("Боксплоты 99-го перцентиля...")
    ti_perf, tv_perf, ti_imperf, tv_imperf = [], [], [], []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref or ref not in per_pos[pos]: continue
        ref_lines = per_pos[pos][ref]
        rp = max_eff_len(ref_lines, True)
        ri = max_eff_len(ref_lines, False)
        for nuc in 'ACGT':
            if nuc == ref or nuc not in per_pos[pos]: continue
            alt_lines = per_pos[pos][nuc]
            ap = max_eff_len(alt_lines, True)
            ai = max_eff_len(alt_lines, False)
            is_ti = {ref, nuc} in ({'A','G'}, {'C','T'})
            if is_ti:
                ti_perf.append(abs(ap-rp)); ti_imperf.append(abs(ai-ri))
            else:
                tv_perf.append(abs(ap-rp)); tv_imperf.append(abs(ai-ri))
    def plot_perc(ax, data, labels, title):
        vals = [np.percentile(d, 99) if d else np.nan for d in data]
        ax.bar(labels, vals, color=['blue','red'])
        for i,v in enumerate(vals):
            if not np.isnan(v): ax.text(i, v, f'{v:.1f}', ha='center', va='bottom')
        ax.set_title(title)
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(10,4))
    plot_perc(ax1, [ti_perf, tv_perf], ['Transition','Transversion'], '99th %ile |Δ| perfect')
    plot_perc(ax2, [ti_imperf, tv_imperf], ['Transition','Transversion'], '99th %ile |Δ| imperfect')
    plt.tight_layout(); plt.savefig(output_dir / 'boxplot_99percentile_delta.png', dpi=150); plt.close()

    # GC/H differences
    gc_ti, gc_tv, hb_ti, hb_tv = [], [], [], []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref or ref not in per_pos[pos]: continue
        ref_lines = per_pos[pos][ref]
        rgc = np.mean([float(l.split('\t')[11]) for l in ref_lines if len(l.split('\t'))>=15])
        rhb = np.mean([int(l.split('\t')[13]) for l in ref_lines if len(l.split('\t'))>=15])
        for nuc in 'ACGT':
            if nuc == ref or nuc not in per_pos[pos]: continue
            alt_lines = per_pos[pos][nuc]
            agc = np.mean([float(l.split('\t')[11]) for l in alt_lines if len(l.split('\t'))>=15])
            ahb = np.mean([int(l.split('\t')[13]) for l in alt_lines if len(l.split('\t'))>=15])
            if np.isnan(agc) or np.isnan(rgc): continue
            is_ti = {ref, nuc} in ({'A','G'}, {'C','T'})
            if is_ti:
                gc_ti.append(abs(agc-rgc)); hb_ti.append(abs(ahb-rhb))
            else:
                gc_tv.append(abs(agc-rgc)); hb_tv.append(abs(ahb-rhb))
    fig, (ax1,ax2) = plt.subplots(1,2, figsize=(10,4))
    plot_perc(ax1, [gc_ti, gc_tv], ['Transition','Transversion'], '99th %ile |ΔGC|')
    plot_perc(ax2, [hb_ti, hb_tv], ['Transition','Transversion'], '99th %ile |ΔH-bonds|')
    plt.tight_layout(); plt.savefig(output_dir / 'boxplot_99percentile_gc_hb.png', dpi=150); plt.close()
    log_print(" Боксплоты сохранены.")

    # ----- Heatmaps -----
    log_print("Heatmap...")
    transitions_max = defaultdict(list)
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref or ref not in per_pos[pos]: continue
        ref_len_perf = max_eff_len(per_pos[pos][ref], True)
        for nuc in 'ACGT':
            if nuc == ref or nuc not in per_pos[pos]: continue
            alt_len_perf = max_eff_len(per_pos[pos][nuc], True)
            transitions_max[(ref, nuc)].append(alt_len_perf - ref_len_perf)
    nuc_order = ['A','C','G','T']
    hm = np.zeros((4,4)); pm = np.ones((4,4))
    for i,ref in enumerate(nuc_order):
        for j,alt in enumerate(nuc_order):
            if (ref,alt) in transitions_max and len(transitions_max[(ref,alt)])>=2:
                d = transitions_max[(ref,alt)]
                _, pval = stats.ttest_1samp(d, 0)
                pm[i,j] = pval
                hm[i,j] = np.mean(d)
    pf = pm.flatten(); sig = np.zeros_like(pf, dtype=bool)
    sidx = np.argsort(pf)
    for r,idx in enumerate(sidx):
        if pf[idx] <= (r+1)/len(pf)*0.05: sig[idx] = True
    sig = sig.reshape(4,4)
    fig, ax = plt.subplots(figsize=(6,5))
    im = ax.imshow(hm, cmap='RdBu_r', aspect='auto')
    ax.set_xticks(range(4)); ax.set_xticklabels(nuc_order)
    ax.set_yticks(range(4)); ax.set_yticklabels(nuc_order)
    plt.colorbar(im, label='Mean Δ max perfect')
    for i in range(4):
        for j in range(4):
            if (nuc_order[i], nuc_order[j]) in transitions_max:
                txt = f'{hm[i,j]:.2f}'
                if sig[i,j]: txt += '*'
                ax.text(j, i, txt, ha='center', va='center', fontsize=9,
                        color='white' if abs(hm[i,j])>np.max(np.abs(hm))/2 else 'black')
    ax.set_title('Heatmap max perfect Δ (* FDR<0.05)')
    plt.tight_layout(); plt.savefig(output_dir / 'heatmap_max_perfect.png', dpi=150); plt.close()

    # all repeats
    transitions_all = defaultdict(list)
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref or ref not in per_pos[pos]: continue
        ref_effs = [int(l.split('\t')[3])-int(l.split('\t')[9]) for l in per_pos[pos][ref] if len(l.split('\t'))>=10]
        avg_ref = np.mean(ref_effs) if ref_effs else 0
        for nuc in 'ACGT':
            if nuc == ref or nuc not in per_pos[pos]: continue
            alt_effs = [int(l.split('\t')[3])-int(l.split('\t')[9]) for l in per_pos[pos][nuc] if len(l.split('\t'))>=10]
            avg_alt = np.mean(alt_effs) if alt_effs else 0
            transitions_all[(ref, nuc)].append(avg_alt - avg_ref)
    hm2 = np.zeros((4,4)); pm2 = np.ones((4,4))
    for i,ref in enumerate(nuc_order):
        for j,alt in enumerate(nuc_order):
            if (ref,alt) in transitions_all and len(transitions_all[(ref,alt)])>=2:
                d = transitions_all[(ref,alt)]
                _, pval = stats.ttest_1samp(d, 0)
                pm2[i,j] = pval
                hm2[i,j] = np.mean(d)
    pf2 = pm2.flatten(); sig2 = np.zeros_like(pf2, dtype=bool)
    sidx2 = np.argsort(pf2)
    for r,idx in enumerate(sidx2):
        if pf2[idx] <= (r+1)/len(pf2)*0.05: sig2[idx] = True
    sig2 = sig2.reshape(4,4)
    fig, ax = plt.subplots(figsize=(6,5))
    im = ax.imshow(hm2, cmap='RdBu_r', aspect='auto')
    ax.set_xticks(range(4)); ax.set_xticklabels(nuc_order)
    ax.set_yticks(range(4)); ax.set_yticklabels(nuc_order)
    plt.colorbar(im, label='Mean Δ all repeats')
    for i in range(4):
        for j in range(4):
            if (nuc_order[i], nuc_order[j]) in transitions_all:
                txt = f'{hm2[i,j]:.2f}'
                if sig2[i,j]: txt += '*'
                ax.text(j, i, txt, ha='center', va='center', fontsize=9,
                        color='white' if abs(hm2[i,j])>np.max(np.abs(hm2))/2 else 'black')
    ax.set_title('Heatmap all repeats Δ (* FDR<0.05)')
    plt.tight_layout(); plt.savefig(output_dir / 'heatmap_all_repeats.png', dpi=150); plt.close()
    log_print(" Heatmap сохранены.")

    # ----- GC/H-bonds средние -----
    log_print("Инфографика GC/H-bonds...")
    rgc_m, rgc_r, agc_m, agc_r = [], [], [], []
    rhb_m, rhb_r, ahb_m, ahb_r = [], [], [], []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref: continue
        for nuc, lines in per_pos[pos].items():
            is_ref = (nuc == ref)
            for l in lines:
                p = l.split('\t')
                if len(p) < 15: continue
                try:
                    mgc=float(p[11]); rgc=float(p[12]); mhb=int(p[13]); rhb=int(p[14])
                except: continue
                if is_ref:
                    rgc_m.append(mgc); rgc_r.append(rgc); rhb_m.append(mhb); rhb_r.append(rhb)
                else:
                    agc_m.append(mgc); agc_r.append(rgc); ahb_m.append(mhb); ahb_r.append(rhb)
    def safe_mean_std(lst): return (np.mean(lst), np.std(lst)) if lst else (np.nan, np.nan)
    fig, axes = plt.subplots(2,2, figsize=(10,8))
    m1,s1 = safe_mean_std(rgc_m); m2,s2 = safe_mean_std(agc_m)
    axes[0,0].bar(['Ref','Alt'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[0,0].set_title('Motif GC')
    m1,s1 = safe_mean_std(rgc_r); m2,s2 = safe_mean_std(agc_r)
    axes[0,1].bar(['Ref','Alt'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[0,1].set_title('Repeat GC')
    m1,s1 = safe_mean_std(rhb_m); m2,s2 = safe_mean_std(ahb_m)
    axes[1,0].bar(['Ref','Alt'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[1,0].set_title('Motif H-bonds')
    m1,s1 = safe_mean_std(rhb_r); m2,s2 = safe_mean_std(ahb_r)
    axes[1,1].bar(['Ref','Alt'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[1,1].set_title('Repeat H-bonds')
    plt.suptitle('Mean ± SD')
    plt.tight_layout(); plt.savefig(output_dir / 'gc_hbonds_infographic.png', dpi=150); plt.close()

    log_print(f"\nВсе результаты в {output_dir}")
    logf.close()

# ============================================================
#  MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(description='Анализ повторяемости мтДНК')
    parser.add_argument('mode', choices=['clean','plot','stats'])
    parser.add_argument('--input_dir', type=str)
    parser.add_argument('--output_dir', type=str)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--reference', type=str)
    args = parser.parse_args()

    if args.mode == 'clean':
        input_dir = Path(args.input_dir or "01KP")
        output_dir = Path(args.output_dir or input_dir / "cleaned")
        print(f"Режим clean: {input_dir} -> {output_dir}")
        files = list(input_dir.glob("*.txt"))
        if not files: print("Нет .txt файлов."); return
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(clean_file, f, output_dir): f for f in files}
            for future in as_completed(futures):
                f = futures[future]
                try: future.result()
                except Exception as e: print(f"   Ошибка {f.name}: {e}")
        print("Все файлы очищены.")

    elif args.mode == 'plot':
        input_dir = Path(args.input_dir or "mut")
        output_dir = Path(args.output_dir or input_dir / "plots")
        run_plot(input_dir, output_dir)

    elif args.mode == 'stats':
        if not args.input_dir or not args.reference:
            print("Укажите --input_dir и --reference"); return
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir or input_dir / "stats")
        run_stats(input_dir, output_dir, args.reference, args.workers)

if __name__ == "__main__":
    main()
