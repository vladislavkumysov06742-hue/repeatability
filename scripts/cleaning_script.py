#!/usr/bin/env python3
"""
Полный пайплайн анализа повторяемости мтДНК.
Режимы: clean, plot, stats, all.
"""

import argparse, subprocess, tempfile, heapq, os, shutil, platform, sys, csv
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from scipy.stats import gaussian_kde, norm, mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# Попытка импорта statsmodels (для FDR в окнах)
try:
    from statsmodels.stats.multitest import multipletests
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False
    print("statsmodels не установлен. Оконный анализ может не работать. Установите: pip install statsmodels")

# ============================================================
#  ГЛОБАЛЬНЫЕ КОНСТАНТЫ
# ============================================================
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

if platform.system() == 'Windows':
    HAS_SYSTEM_SORT = False
else:
    HAS_SYSTEM_SORT = shutil.which('sort') is not None

# ============================================================
#  ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
# ============================================================
def norm_interval(start, end, MT_LEN):
    if start <= end:
        return start, end
    return start, end + MT_LEN

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
#  ЗАГРУЗКА РЕФЕРЕНСА
# ============================================================
def read_reference(path):
    seq = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'): continue
            seq.append(line.strip().upper())
    return ''.join(seq)

# ============================================================
#  CLEAN
# ============================================================
def clean_file(file_path, output_dir, MT_LEN):
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
    records = []
    for l in valid:
        p = l.split('\t')
        m_s, m_e = norm_interval(int(p[4]), int(p[5]), MT_LEN)
        r_s, r_e = norm_interval(int(p[7]), int(p[8]), MT_LEN)
        records.append((m_s, m_e, r_s, r_e, l))
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
def interval_in_major_arc(start, end, MT_LEN):
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

def get_max_perfect_repeat_details(file_path, MT_LEN):
    max_len = 0; best_seq = ""; best_start = 0; best_end = 0
    with file_path.open() as f:
        f.readline()
        for line in f:
            if not line.strip(): continue
            p = line.strip().split('\t')
            if len(p) < 10: continue
            try: hamming = int(p[9]); rstart = int(p[7]); rend = int(p[8])
            except: continue
            if hamming == 0 and interval_in_major_arc(rstart, rend, MT_LEN):
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

def run_plot(input_dir, output_dir, MT_LEN):
    print(f"Режим plot: {input_dir}")
    delta_values, details, missing = {}, {}, []
    for snp in SNP_NAMES:
        ref_path = input_dir / f"01KP.{snp}.{REF_SUFFIX[snp]}.txt"
        alt_path = input_dir / f"01KP.{snp}.{ALT_SUFFIX[snp]}.txt"
        if not ref_path.exists(): missing.append(str(ref_path))
        elif not alt_path.exists(): missing.append(str(alt_path))
        else:
            rlen, rseq, rs, re = get_max_perfect_repeat_details(ref_path, MT_LEN)
            alen, aseq, as_, ae = get_max_perfect_repeat_details(alt_path, MT_LEN)
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
#  STATS (MAD Z-score + все графики + окна)
# ============================================================
def load_cleaned_files(input_dir, reference_seq, log_print, MT_LEN):
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
        count_ge7 = perfect_ge7 = count_imperf_ge10 = 0
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
                if hamming == 0:
                    perfect += 1
                    if eff >= 7: perfect_ge7 += 1
                else:
                    if eff >= 10: count_imperf_ge10 += 1
                if eff >= 7: count_ge7 += 1
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
        metrics[pos][nuc] = {
            'total': total, 'perfect': perfect, 'max_eff': max_eff, 'sum_eff': sum_eff,
            'count_ge7': count_ge7, 'perfect_ge7': perfect_ge7, 'imperfect_ge10': count_imperf_ge10
        }
        if 1 <= pos <= len(reference_seq):
            ref_nuc[pos] = reference_seq[pos-1]
        n_files += 1
    log_print(f"Загружено файлов: {n_files}")
    log_print(f"Всего повторов (raw): {len(all_repeats_raw)}")
    return per_pos, ref_nuc, all_repeats_raw, all_repeats_weighted, metrics

def mad_zscore_pvalue(observed, background):
    """Робастный MAD Z-score с двусторонним p-value."""
    if len(background) < 5:
        return 1.0, {'test': 'MAD-Z', 'params': 'N<5'}
    bg = np.array(background, dtype=float)
    med = np.median(bg)
    mad = np.median(np.abs(bg - med))
    if mad == 0:
        # если MAD=0, используем обычный z-score
        std = np.std(bg)
        if std == 0:
            return 1.0 if observed == med else 0.0, {'test': 'Z-score', 'params': f'std={std:.3f}'}
        z = (observed - med) / std
        pval = 2 * norm.sf(abs(z))
        return pval, {'test': 'Z-score', 'params': f'std={std:.3f}'}
    # Модифицированный Z-score (Iglewicz & Hoaglin)
    z = 0.6745 * (observed - med) / mad
    pval = 2 * norm.sf(abs(z))
    return pval, {'test': 'MAD-Z', 'params': f'med={med:.2f}, MAD={mad:.2f}'}

def add_significance_lines(ax, pvals, log_print=None, test_name=''):
    n = len(pvals)
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='Nominal α=0.05')
    if n > 0 and STATSMODELS_AVAILABLE:
        _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
        sig = np.min(np.array(pvals)[qvals < 0.05]) if np.any(qvals < 0.05) else None
        if sig is not None:
            ax.axhline(-np.log10(sig), color='red', linestyle='--', alpha=0.5, label='FDR threshold')
            if log_print: log_print(f"   FDR порог: p≤{sig:.2e}")
        else:
            if log_print: log_print("   FDR: значимых точек нет")
    ax.legend(title=test_name)

def windowed_analysis(per_pos, ref_nuc, metrics, metric_name, MT_LEN, window_size=200, step_size=100, n_perm=1000, log_print=None):
    positions_all = sorted(per_pos.keys())
    n_pos = len(positions_all)
    ref_vals = np.zeros(n_pos)
    for i, pos in enumerate(positions_all):
        ref = ref_nuc.get(pos)
        if ref and ref in per_pos[pos]:
            ref_vals[i] = metrics[pos][ref].get(metric_name, 0)
    windows = []
    for start in range(1, MT_LEN - window_size + 2, step_size):
        end = start + window_size - 1
        indices = [i for i, p in enumerate(positions_all) if start <= p <= end]
        if len(indices) == 0: continue
        alt_sum = 0.0
        for i in indices:
            pos = positions_all[i]
            ref = ref_nuc.get(pos)
            for nuc in 'ACGT':
                if nuc == ref or nuc not in per_pos[pos]: continue
                alt_sum += metrics[pos][nuc].get(metric_name, 0)
        windows.append((start, end, indices, alt_sum))
    log_print(f"   Оконный анализ ({metric_name}) для {len(windows)} окон...")
    pvals = []
    for win_idx, (start, end, indices, alt_sum) in enumerate(windows):
        perm_sums = []
        for _ in range(n_perm):
            rand_indices = np.random.choice(n_pos, size=len(indices), replace=False)
            perm_sum = np.sum(ref_vals[rand_indices])
            perm_sums.append(perm_sum)
        p = (np.sum(np.array(perm_sums) >= alt_sum) + 1) / (n_perm + 1)
        pvals.append(p)
        if (win_idx+1) % 20 == 0:
            log_print(f"   ... обработано {win_idx+1}/{len(windows)} окон")
    if STATSMODELS_AVAILABLE:
        _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
    else:
        qvals = np.ones_like(pvals)
    positions = [(s+e)//2 for s,e,_ in windows]
    logp = -np.log10(np.maximum(pvals, 1e-300))
    fig, ax = plt.subplots(figsize=(14,5))
    ax.scatter(positions, logp, c='black', s=5)
    sig_wins = [i for i, q in enumerate(qvals) if q < 0.05]
    if sig_wins:
        ax.scatter([positions[i] for i in sig_wins], [logp[i] for i in sig_wins], c='red', s=10)
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='Nominal α=0.05')
    ax.axhline(-np.log10(np.max(pvals[np.array(qvals) < 0.05] if sig_wins else [1])), color='red', linestyle='--', label='FDR threshold')
    ax.set_title(f'Windowed Manhattan – {metric_name} (window {window_size}bp, step {step_size}bp)')
    ax.set_xlabel('Position'); ax.set_ylabel('-log10(p)')
    ax.legend(title='Permutation test + FDR')
    plt.tight_layout()
    return fig, ax

import heapq

def compute_max_coverage(all_repeats, MT_LEN):
    max_len = np.zeros(MT_LEN + 1, dtype=np.int32)   # теперь индекс 0 – не используется
    if not all_repeats:
        return max_len
    events = []
    for start, end, eff in all_repeats:
        if end > MT_LEN:
            end = MT_LEN
        events.append((start, 1, eff))
        events.append((end + 1, -1, eff))
    events.sort(key=lambda x: (x[0], -x[1]))
    max_heap = []
    remove_heap = []
    event_idx = 0
    n_events = len(events)
    for pos in range(1, MT_LEN + 1):
        while event_idx < n_events and events[event_idx][0] == pos:
            _, typ, eff = events[event_idx]
            if typ == 1:
                heapq.heappush(max_heap, -eff)
            else:
                heapq.heappush(remove_heap, -eff)
            event_idx += 1
        while remove_heap and max_heap and remove_heap[0] == max_heap[0]:
            heapq.heappop(remove_heap)
            heapq.heappop(max_heap)
        if max_heap:
            max_len[pos] = -max_heap[0]   # записываем в позицию pos
    return max_len

def compute_coverage(all_repeats, MT_LEN, weighted=False):
    """Быстрое вычисление покрытия через разностный массив и cumsum.
    Возвращает массив длиной MT_LEN+1, где индекс 0 игнорируется,
    а индексы 1..MT_LEN соответствуют позициям."""
    cov = np.zeros(MT_LEN + 2, dtype=np.float64)
    if not all_repeats:
        return cov[:MT_LEN+1]

    starts = np.array([item[0] for item in all_repeats], dtype=np.int32)
    ends   = np.array([item[1] for item in all_repeats], dtype=np.int32)
    ends = np.clip(ends, 0, MT_LEN)

    if weighted:
        weights = np.array([item[2] for item in all_repeats], dtype=np.float64)
    else:
        weights = np.ones(len(starts), dtype=np.float64)

    np.add.at(cov, starts, weights)

    valid = ends + 1 <= MT_LEN
    if valid.any():
        np.add.at(cov, ends[valid] + 1, -weights[valid])

    cov = np.cumsum(cov)[:MT_LEN+1]
    return cov


def run_stats(input_dir, output_dir, ref_path, workers, window_size=0, snv_file=None):
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "stats_run.log"
    logf = open(log_path, 'w', encoding='utf-8')
    def log_print(msg):
        print(msg)
        print(msg, file=logf)

    log_print("Загрузка референсной последовательности...")
    ref_seq = read_reference(ref_path)
    MT_LEN = len(ref_seq)
    log_print(f"Длина референса: {MT_LEN}")

    log_print("Чтение очищенных файлов...")
    per_pos, ref_nuc, all_repeats_raw, all_repeats_weighted, metrics = load_cleaned_files(input_dir, ref_seq, log_print, MT_LEN)
    log_print(f"Позиций с данными: {len(per_pos)}")

    # ----- Фон -----
    ref_metrics_names = ['total', 'perfect', 'max_eff', 'sum_eff',
                         'count_ge7', 'perfect_ge7', 'imperfect_ge10']
    ref_backgrounds = {name: [] for name in ref_metrics_names}
    for pos, nucs in per_pos.items():
        ref = ref_nuc.get(pos)
        if ref and ref in nucs:
            m = metrics[pos][ref]
            for name in ref_metrics_names:
                ref_backgrounds[name].append(m[name])
    log_print("Фоны подготовлены.")

    # ----- Манхеттены (MAD-Z тест) -----
    titles = {
        'total': 'Total repeats', 'perfect': 'Perfect repeats',
        'count_ge7': 'Repeats ≥7 (all)', 'perfect_ge7': 'Perfect repeats ≥7',
        'imperfect_ge10': 'Imperfect repeats ≥10',
        'max_eff': 'Max effective length', 'sum_eff': 'Sum effective length'
    }
    plot_order = ['total', 'perfect', 'count_ge7', 'perfect_ge7', 'imperfect_ge10', 'max_eff', 'sum_eff']
    fig, axes = plt.subplots(len(plot_order), 1, figsize=(14, 5*len(plot_order)))
    for idx, metric in enumerate(plot_order):
        ax = axes[idx]
        bg = ref_backgrounds[metric]
        title = titles[metric]
        # Порог фильтрации (MAD)
        all_deltas = []
        for pos in sorted(per_pos.keys()):
            ref = ref_nuc.get(pos)
            if not ref or ref not in per_pos[pos]: continue
            ref_val = metrics[pos][ref][metric]
            for nuc in 'ACGT':
                if nuc == ref or nuc not in per_pos[pos]: continue
                all_deltas.append(metrics[pos][nuc][metric] - ref_val)
        if len(all_deltas) > 1:
            med_delta = np.median(all_deltas)
            mad_delta = np.median(np.abs(np.array(all_deltas) - med_delta))
            threshold = 1.5 * mad_delta if metric != 'max_eff' else 1.0 * mad_delta
        else:
            threshold = 0
        log_print(f"\n {title}: MAD={mad_delta:.2f}, порог |Δ| > {threshold:.2f}")
        positions, pvals = [], []
        tested_total, tested_passed = 0, 0
        last_report = 0
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
                pval, info = mad_zscore_pvalue(alt_val, bg)
                positions.append(pos)
                pvals.append(pval)
                if tested_passed - last_report >= 1000:
                    log_print(f"   ... обработано {tested_passed} тестов (pmin={min(pvals):.2e})")
                    last_report = tested_passed
        log_print(f"   Всего тестов: {tested_total}, после фильтрации: {tested_passed}")
        if pvals:
            info_list = sorted(zip(positions, pvals), key=lambda x: x[1])
            log_print("   Топ-5 p-values:")
            for pos, p in info_list[:5]:
                log_print(f"     pos={pos}, p={p:.2e}")
            log_print(f"   Минимальное p-value до FDR: {info_list[0][1]:.2e}")
            # Тест для легенды
            _, info_test = mad_zscore_pvalue(bg[0], bg)  # чтобы получить тип теста
            test_name = info_test.get('test', 'MAD-Z')
        else:
            log_print("   Нет точек после фильтрации.")
        if not pvals:
            ax.text(0.5,0.5,'No points', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{title}')
            continue
        logp = -np.log10(np.maximum(pvals, 1e-300))
        ax.scatter(positions, logp, c='black', s=5)
        ax.set_title(f'Manhattan plot – {title}')
        add_significance_lines(ax, pvals, log_print=log_print, test_name=test_name)
    plt.tight_layout()
    plt.savefig(output_dir / 'manhattan_metrics.png', dpi=150); plt.close()
    log_print("\n Манхеттены метрик сохранены.")

    # ----- Покрытия (MAD-Z для coverage) -----
    cov_raw = compute_coverage(all_repeats_raw, MT_LEN, weighted=False)
    cov_weighted = compute_coverage(all_repeats_weighted, MT_LEN, weighted=True)
    for cov, name in [(cov_raw, 'Raw coverage'), (cov_weighted, 'Weighted coverage')]:
        bg_cov = cov[1:]
        pvals_cov = []
        for pos in range(1, MT_LEN+1):
            pval, _ = mad_zscore_pvalue(cov[pos], bg_cov)
            pvals_cov.append(pval)
        logp_cov = -np.log10(np.maximum(pvals_cov, 1e-300))
        fig, ax = plt.subplots(figsize=(14,5))
        ax.scatter(range(1, MT_LEN+1), logp_cov, c='steelblue', s=2)
        ax.set_title(f'Manhattan – {name}')
        add_significance_lines(ax, pvals_cov, log_print=log_print, test_name='MAD-Z')
        plt.tight_layout()
        plt.savefig(output_dir / f'manhattan_{name.replace(" ","_")}.png', dpi=150); plt.close()

    # Сглаженный coverage
    fig, ax = plt.subplots(figsize=(14,4))
    ax.plot(range(1,MT_LEN+1), cov_weighted[1:], linewidth=0.3, color='darkblue', alpha=0.5, label='Raw')
    smooth = np.convolve(cov_weighted[1:], np.ones(500)/500, mode='same')
    ax.plot(range(1,MT_LEN+1), smooth, linewidth=1.5, color='red', label='Smoothed')
    ax.legend(); ax.set_title('Weighted coverage along mtDNA')
    plt.tight_layout(); plt.savefig(output_dir / 'coverage_plot.png', dpi=150); plt.close()

    # ----- Гистограммы (быстрые, без KDE) -----
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,5))
    vals = cov_weighted[1:]
    ax1.hist(vals, bins=50, density=True, color='lightblue', edgecolor='black', alpha=0.7)
    ax1.set_title('Weighted coverage')
    ax1.text(0.95,0.95, f'Mean={np.mean(vals):.1f}\nMedian={np.median(vals):.1f}\nSkew={stats.skew(vals):.2f}\nKurt={stats.kurtosis(vals):.2f}',
             transform=ax1.transAxes, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.8))

    all_cnt = [metrics[pos][nuc]['total'] for pos in per_pos for nuc in per_pos[pos]]
    ax2.hist(all_cnt, bins=50, density=True, color='lightgreen', edgecolor='black', alpha=0.7)
    ax2.set_title('Repeats per file')
    ax2.text(0.95,0.95, f'Mean={np.mean(all_cnt):.1f}\nMedian={np.median(all_cnt):.1f}\nSkew={stats.skew(all_cnt):.2f}\nKurt={stats.kurtosis(all_cnt):.2f}',
             transform=ax2.transAxes, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.8))
    plt.tight_layout(); plt.savefig(output_dir / 'histograms.png', dpi=150); plt.close()

    
    # ----- Scatter plots (jitter + регрессия) -----
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
    ax.axhline(0, color='gray', linestyle='--'); ax.set_title('Perfect: Δ vs ref')
    plt.tight_layout(); plt.savefig(output_dir / 'scatter_delta_perfect.png', dpi=150); plt.close()
    # delta imperfect
    x_i = np.array([d[2] for d in scatter_data])
    y_i = np.array([d[3] for d in scatter_data])
    fig, ax = plt.subplots(figsize=(8,6))
    scatter_with_regression(ax, x_i, y_i, 'red')
    ax.set_xlabel('Ref max eff. imperfect'); ax.set_ylabel('Δ imperfect')
    ax.axhline(0, color='gray', linestyle='--'); ax.set_title('Imperfect: Δ vs ref')
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

    # ----- Боксплоты 99-го перцентиля (Δ, GC, H-bonds, longest) -----
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
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    plot_perc(ax1, [ti_perf, tv_perf], ['Transition','Transversion'], '99th %ile |Δ| perfect')
    plot_perc(ax2, [ti_imperf, tv_imperf], ['Transition','Transversion'], '99th %ile |Δ| imperfect')
    plt.tight_layout(); plt.savefig(output_dir / 'boxplot_99percentile_delta.png', dpi=150); plt.close()

    # GC/H differences (all repeats)
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
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    plot_perc(ax1, [gc_ti, gc_tv], ['Transition','Transversion'], '99th %ile |ΔGC|')
    plot_perc(ax2, [hb_ti, hb_tv], ['Transition','Transversion'], '99th %ile |ΔH-bonds|')
    plt.tight_layout(); plt.savefig(output_dir / 'boxplot_99percentile_gc_hb.png', dpi=150); plt.close()

    # GC/H для самых длинных повторов (longest)
    gc_long_ti, gc_long_tv, hb_long_ti, hb_long_tv = [], [], [], []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref or ref not in per_pos[pos]: continue
        ref_lines = per_pos[pos][ref]
        max_eff_ref = max_eff_len(ref_lines, False)
        for nuc in 'ACGT':
            if nuc == ref or nuc not in per_pos[pos]: continue
            alt_lines = per_pos[pos][nuc]
            max_eff_alt = max_eff_len(alt_lines, False)
            alt_gc_long = []; alt_hb_long = []
            for l in alt_lines:
                p = l.split('\t')
                try: eff = int(p[3]) - int(p[9])
                except: continue
                if eff == max_eff_alt and len(p) >= 15:
                    alt_gc_long.append(float(p[11])); alt_hb_long.append(int(p[13]))
            ref_gc_long = []; ref_hb_long = []
            for l in ref_lines:
                p = l.split('\t')
                try: eff = int(p[3]) - int(p[9])
                except: continue
                if eff == max_eff_ref and len(p) >= 15:
                    ref_gc_long.append(float(p[11])); ref_hb_long.append(int(p[13]))
            if alt_gc_long and ref_gc_long:
                diff_gc = abs(np.mean(alt_gc_long) - np.mean(ref_gc_long))
                diff_hb = abs(np.mean(alt_hb_long) - np.mean(ref_hb_long))
                is_ti = {ref, nuc} in ({'A','G'}, {'C','T'})
                if is_ti:
                    gc_long_ti.append(diff_gc); hb_long_ti.append(diff_hb)
                else:
                    gc_long_tv.append(diff_gc); hb_long_tv.append(diff_hb)
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    plot_perc(ax1, [gc_long_ti, gc_long_tv], ['Transition','Transversion'], '99th %ile |ΔGC| (longest)')
    plot_perc(ax2, [hb_long_ti, hb_long_tv], ['Transition','Transversion'], '99th %ile |ΔH-bonds| (longest)')
    plt.tight_layout(); plt.savefig(output_dir / 'boxplot_99percentile_gc_hb_longest.png', dpi=150); plt.close()

    # ----- Heatmap (max perfect Δ и all repeats Δ) -----
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

    # ----- GC/H-bonds средние -----
    rgc_m, rgc_r, agc_m, agc_r = [], [], [], []
    rhb_m, rhb_r, ahb_m, ahb_r = [], [], [], []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref: continue
        for nuc, lines in per_pos[pos].items():
            is_ref = (nuc == ref)
            for l in lines:
                p = l.split('\t')
                if len(p) < 14: continue          # было < 15  --  строка из 14 полей
                try:
                    mgc = float(p[10])            # motif_GC
                    rgc = float(p[11])            # repeat_GC
                    mhb = int(p[12])              # motif_H_bonds
                    rhb = int(p[13])              # repeat_H_bonds
                except:
                    continue
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

    # ----- Зависимость от major arc (расстояние) -----
    distances = []
    met_vals = []
    for pos in per_pos.keys():
        if not ref_nuc.get(pos): continue
        dist_start = min(abs(pos - MAJOR_ARC_START), MT_LEN - abs(pos - MAJOR_ARC_START))
        dist_end = min(abs(pos - MAJOR_ARC_END), MT_LEN - abs(pos - MAJOR_ARC_END))
        dist = min(dist_start, dist_end)
        ref_n = ref_nuc[pos]
        if ref_n in per_pos[pos]:
            # ИСПРАВЛЕНИЕ: заменён несуществующий 'long_count' на 'count_ge7'
            val = metrics[pos][ref_n]['count_ge7']
            distances.append(dist)
            met_vals.append(val)
    bins = np.arange(0, max(distances)+100, 100)
    bin_centers = (bins[:-1] + bins[1:])/2
    binned_mean = [np.mean(np.array(met_vals)[(np.array(distances)>=bins[i]) & (np.array(distances)<bins[i+1])]) for i in range(len(bins)-1)]
    fig, ax = plt.subplots(figsize=(10,4))
    ax.plot(bin_centers, binned_mean, 'o-')
    ax.set_xlabel('Distance to nearest major arc boundary (bp)')
    ax.set_ylabel('Mean count_ge7 (ref)')
    ax.set_title('Repeatability vs distance from major arc boundaries')
    plt.tight_layout(); plt.savefig(output_dir / 'major_arc_distance.png', dpi=150); plt.close()

    # ----- Оконный анализ -----
    if window_size > 0 and STATSMODELS_AVAILABLE:
        log_print("\n=== Оконный анализ ===")
        for metric in ['total', 'perfect', 'count_ge7', 'perfect_ge7', 'imperfect_ge10']:
            fig, _ = windowed_analysis(per_pos, ref_nuc, metrics, metric, MT_LEN, window_size=window_size, log_print=log_print)
            fig.savefig(output_dir / f'manhattan_window_{metric}.png', dpi=150)
            plt.close(fig)

    # ----- Forest plot для CSV с мутациями -----
    if snv_file:
        log_print(f"\n=== Forest plot для мутаций из {snv_file} ===")
        snv_data = []
        with open(snv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pos = int(row['position'])
                ref_allele = row['ref_allele'].upper()
                alt_allele = row['alt_allele'].upper()
                if pos in per_pos and ref_allele in per_pos[pos] and alt_allele in per_pos[pos]:
                    ref_lines = per_pos[pos][ref_allele]
                    alt_lines = per_pos[pos][alt_allele]
                    rp = max_eff_len(ref_lines, True)
                    ap = max_eff_len(alt_lines, True)
                    snv_data.append((f"{pos}{ref_allele}>{alt_allele}", ap - rp))
        if snv_data:
            labels, deltas = zip(*snv_data)
            sorted_idx = np.argsort(deltas)
            labels = [labels[i] for i in sorted_idx]
            deltas = [deltas[i] for i in sorted_idx]
            fig, ax = plt.subplots(figsize=(8, max(6, len(labels)*0.3)))
            ax.barh(range(len(labels)), deltas, color=[COLOR_POS_BAR if d>0 else COLOR_NEG_BAR for d in deltas])
            ax.set_yticks(range(len(labels)))
            ax.set_yticklabels(labels, fontsize=7)
            ax.axvline(0, color='black', linestyle='--')
            ax.set_xlabel('Δ max perfect length')
            ax.set_title('Forest plot for SNV (Δ max perfect)')
            plt.tight_layout()
            plt.savefig(output_dir / 'forest_snv.png', dpi=150); plt.close()
            log_print(f" Forest plot сохранён ({len(labels)} мутаций).")

    log_print(f"\nВсе результаты в {output_dir}")
    logf.close()

# ============================================================
#  ALL
# ============================================================
def run_all(input_dir, output_dir, ref_path, workers, window_size=0, snv_file=None):
    cleaned_dir = output_dir / 'cleaned'
    stats_dir = output_dir / 'stats'
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

# ============================================================
#  MAIN
# ============================================================
def main():
    parser = argparse.ArgumentParser(description='Анализ повторяемости мтДНК')
    parser.add_argument('mode', choices=['clean','plot','stats','all'])
    parser.add_argument('--input_dir', type=str)
    parser.add_argument('--output_dir', type=str)
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--reference', type=str, required=True)
    parser.add_argument('--window_size', type=int, default=0)
    parser.add_argument('--snv_file', type=str)
    args = parser.parse_args()

    if args.mode == 'clean':
        input_dir = Path(args.input_dir or "01KP")
        output_dir = Path(args.output_dir or input_dir / "cleaned")
        ref_seq = read_reference(args.reference)
        MT_LEN = len(ref_seq)
        print(f"Режим clean: {input_dir} -> {output_dir}")
        txt_files = list(input_dir.glob("*.txt"))
        if not txt_files: print("Нет .txt файлов."); return
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(clean_file, f, output_dir, MT_LEN): f for f in txt_files}
            for future in as_completed(futures):
                f = futures[future]
                try: future.result()
                except Exception as e: print(f"   Ошибка {f.name}: {e}")
        print("Все файлы очищены.")

    elif args.mode == 'plot':
        input_dir = Path(args.input_dir or "mut")
        output_dir = Path(args.output_dir or input_dir / "plots")
        ref_seq = read_reference(args.reference)
        MT_LEN = len(ref_seq)
        run_plot(input_dir, output_dir, MT_LEN)

    elif args.mode == 'stats':
        if not args.input_dir: print("Укажите --input_dir"); return
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir or input_dir / "stats")
        run_stats(input_dir, output_dir, args.reference, args.workers,
                  window_size=args.window_size, snv_file=args.snv_file)

    elif args.mode == 'all':
        if not args.input_dir: print("Укажите --input_dir"); return
        input_dir = Path(args.input_dir)
        output_dir = Path(args.output_dir or input_dir / "results")
        run_all(input_dir, output_dir, args.reference, args.workers,
                window_size=args.window_size, snv_file=args.snv_file)

if __name__ == "__main__":
    main()
