import heapq
import csv
import gzip
import pickle
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import gaussian_kde, norm, mannwhitneyu
from collections import defaultdict
from pathlib import Path

from .reference import read_reference
from .constants import MAJOR_ARC_START, MAJOR_ARC_END, COLOR_POS_BAR, COLOR_NEG_BAR

warnings.filterwarnings('ignore')

try:
    from statsmodels.stats.multitest import multipletests
    STATSMODELS_AVAILABLE = True
except ImportError:
    STATSMODELS_AVAILABLE = False
    print("statsmodels не установлен. Оконный анализ может не работать. Установите: pip install statsmodels")

# ------------------------------------------------------------
# Загрузка данных (с вычислением покрытия на лету)
# ------------------------------------------------------------
def load_cleaned_files(input_dir, reference_seq, log_print, MT_LEN):
    per_pos = defaultdict(dict)
    ref_nuc = {}
    metrics = defaultdict(dict)
    n_files = 0
    
    cov_raw = np.zeros(MT_LEN + 2, dtype=np.int32)
    cov_weighted = np.zeros(MT_LEN + 2, dtype=np.float64)
    
    for fpath in Path(input_dir).glob("01KP.*.*.txt"):
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
                    cov_raw[rs] += 1
                    cov_raw[re+1] -= 1
                    cov_weighted[rs] += eff
                    cov_weighted[re+1] -= eff
                else:
                    cov_raw[rs] += 1
                    cov_raw[MT_LEN+1] -= 1
                    cov_raw[1] += 1
                    cov_raw[re+1] -= 1
                    cov_weighted[rs] += eff
                    cov_weighted[MT_LEN+1] -= eff
                    cov_weighted[1] += eff
                    cov_weighted[re+1] -= eff
        per_pos[pos][nuc] = lines
        metrics[pos][nuc] = {
            'total': total, 'perfect': perfect, 'max_eff': max_eff, 'sum_eff': sum_eff,
            'count_ge7': count_ge7, 'perfect_ge7': perfect_ge7, 'imperfect_ge10': count_imperf_ge10
        }
        if 1 <= pos <= len(reference_seq):
            ref_nuc[pos] = reference_seq[pos-1]
        n_files += 1
    
    cov_raw = np.cumsum(cov_raw)[:MT_LEN+1]
    cov_weighted = np.cumsum(cov_weighted)[:MT_LEN+1]
    
    log_print(f"Загружено файлов: {n_files}")
    log_print(f"Покрытие вычислено (без сохранения списков повторов).")
    return per_pos, ref_nuc, cov_raw, cov_weighted, metrics

def load_cleaned_files_cached(input_dir, reference_seq, log_print, MT_LEN, cache_path):
    cache_path = Path(cache_path)
    if cache_path.exists():
        log_print(f"Загрузка кэша из {cache_path} ...")
        with gzip.open(cache_path, 'rb') as f:
            data = pickle.load(f)
        log_print("Кэш загружен.")
        return data
    log_print("Кэш не найден, полная загрузка...")
    per_pos, ref_nuc, cov_raw, cov_weighted, metrics = load_cleaned_files(
        input_dir, reference_seq, log_print, MT_LEN)
    log_print(f"Сохранение кэша в {cache_path} ...")
    with gzip.open(cache_path, 'wb') as f:
        pickle.dump((per_pos, ref_nuc, cov_raw, cov_weighted, metrics),
                    f, protocol=pickle.HIGHEST_PROTOCOL)
    log_print("Кэш сохранён.")
    return per_pos, ref_nuc, cov_raw, cov_weighted, metrics

# ------------------------------------------------------------
# Статистические функции
# ------------------------------------------------------------
def mad_zscore_pvalue(observed, background):
    if len(background) < 5:
        return 1.0, {'test': 'MAD-Z', 'params': 'N<5'}
    bg = np.array(background, dtype=float)
    med = np.median(bg)
    mad = np.median(np.abs(bg - med))
    if mad == 0:
        std = np.std(bg)
        if std == 0:
            return 1.0 if observed == med else 0.0, {'test': 'Z-score', 'params': f'std={std:.3f}'}
        z = (observed - med) / std
        pval = 2 * norm.sf(abs(z))
        return pval, {'test': 'Z-score', 'params': f'std={std:.3f}'}
    z = 0.6745 * (observed - med) / mad
    pval = 2 * norm.sf(abs(z))
    return pval, {'test': 'MAD-Z', 'params': f'med={med:.2f}, MAD={mad:.2f}'}

def add_significance_lines(ax, pvals, log_print=None, test_name=''):
    ax.axhline(-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='Номинальный α=0,05')
    if len(pvals) > 0 and STATSMODELS_AVAILABLE:
        _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
        sig = np.min(np.array(pvals)[qvals < 0.05]) if np.any(qvals < 0.05) else None
        if sig is not None:
            ax.axhline(-np.log10(sig), color='red', linestyle='--', alpha=0.5, label='Порог FDR')
            if log_print: log_print(f"   FDR порог: p≤{sig:.2e}")
        else:
            if log_print: log_print("   FDR: значимых точек нет")
    ax.legend(title=test_name)

def windowed_analysis(per_pos, ref_nuc, metrics, metric_name, MT_LEN, window_size=200, step_size=100,
                      n_perm=1000, log_print=None):
    metric_names_ru = {
        'total': 'общее число повторов',
        'perfect': 'число совершенных повторов',
        'count_ge7': 'число повторов с эффективной длиной ≥7',
        'perfect_ge7': 'число совершенных повторов с длиной ≥7',
        'imperfect_ge10': 'число несовершенных повторов с эффективной длиной ≥10'
    }
    ru_name = metric_names_ru.get(metric_name, metric_name)
    
    positions_all = sorted(per_pos.keys())
    n_pos = len(positions_all)
    ref_vals = np.zeros(n_pos)
    for i, pos in enumerate(positions_all):
        ref = ref_nuc.get(pos)
        if ref and ref in per_pos[pos]:
            ref_vals[i] = metrics[pos][ref].get(metric_name, 0)
    windows = []
    alt_sums = []
    for start in range(1, MT_LEN - window_size + 2, step_size):
        end = start + window_size - 1
        indices = [i for i, p in enumerate(positions_all) if start <= p <= end]
        if len(indices) == 0:
            continue
        alt_sum = 0.0
        for i in indices:
            pos = positions_all[i]
            ref = ref_nuc.get(pos)
            for nuc in 'ACGT':
                if nuc == ref or nuc not in per_pos[pos]:
                    continue
                alt_sum += metrics[pos][nuc].get(metric_name, 0)
        alt_sums.append(alt_sum)
        windows.append((start, end, indices))
    log_print(f"   Оконный анализ для {ru_name}: {len(windows)} окон...")
    pvals = []
    for win_idx, (start, end, indices) in enumerate(windows):
        alt_sum = alt_sums[win_idx]
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
    # barplot
    x_labels = [f"{s}-{e}" for (s,e,_) in windows]
    step_label = max(1, len(x_labels)//20)
    xticks = range(0, len(x_labels), step_label)
    xticklabels = [x_labels[i] if i < len(x_labels) else '' for i in xticks]
    fig, ax = plt.subplots(figsize=(max(10, len(x_labels)*0.3), 6))
    ax.bar(range(len(alt_sums)), alt_sums, width=0.8, color='steelblue')
    for i, (q, val) in enumerate(zip(qvals, alt_sums)):
        if q < 0.05:
            ax.text(i, val + 0.02*max(alt_sums), '*', ha='center', fontsize=12, fontweight='bold')
    ax.set_xlabel('Границы окна (начало-конец) п.н.')
    ax.set_ylabel(ru_name)
    ax.set_title(f'Оконный анализ: {ru_name} (окно {window_size} п.н.)')
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=45, ha='right', fontsize=8)
    plt.tight_layout()
    return fig, ax

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

# ------------------------------------------------------------
# Анализ одного генома (все графики)
# ------------------------------------------------------------
def analyze_single_genome(input_dir, ref_path, output_dir, workers, window_size, snv_file, genome_name, log_print):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_print(f"\n=== Анализ генома: {genome_name} ===")
    
    ref_seq = read_reference(ref_path)
    MT_LEN = len(ref_seq)
    log_print(f"Длина референса: {MT_LEN} п.н.")
    
    cache_path = output_dir / 'cleaned_data.pkl.gz'
    per_pos, ref_nuc, cov_raw, cov_weighted, metrics = load_cleaned_files_cached(
        input_dir, ref_seq, log_print, MT_LEN, cache_path)
    
    log_print(f"Позиций с данными: {len(per_pos)}")
    
    # ----- Фоны -----
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
    
    # ----- Манхэттены метрик -----
    titles_ru = {
        'total': 'общее число повторов',
        'perfect': 'число совершенных повторов',
        'count_ge7': 'число повторов с эффективной длиной ≥7',
        'perfect_ge7': 'число совершенных повторов длиной ≥7',
        'imperfect_ge10': 'число несовершенных повторов с эффективной длиной ≥10',
        'max_eff': 'максимальная эффективная длина повтора',
        'sum_eff': 'сумма эффективных длин всех повторов'
    }
    plot_order = ['total', 'perfect', 'count_ge7', 'perfect_ge7', 'imperfect_ge10', 'max_eff', 'sum_eff']
    fig, axes = plt.subplots(len(plot_order), 1, figsize=(14, 5*len(plot_order)))
    for idx, metric in enumerate(plot_order):
        ax = axes[idx]
        bg = ref_backgrounds[metric]
        title_ru = titles_ru[metric]
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
        log_print(f"\n {title_ru}: MAD={mad_delta:.2f}, порог |Δ| > {threshold:.2f}")
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
                    log_print(f"   ... обработано {tested_passed} тестов (min p={min(pvals):.2e})")
                    last_report = tested_passed
        log_print(f"   Всего тестов: {tested_total}, после фильтрации: {tested_passed}")
        if pvals:
            info_list = sorted(zip(positions, pvals), key=lambda x: x[1])
            log_print("   Топ-5 p-values:")
            for pos, p in info_list[:5]:
                log_print(f"     pos={pos}, p={p:.2e}")
            _, info_test = mad_zscore_pvalue(bg[0], bg)
            test_name = info_test.get('test', 'MAD-Z')
        else:
            log_print("   Нет точек после фильтрации.")
        if not pvals:
            ax.text(0.5,0.5,'Нет точек', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f'{title_ru}')
            continue
        logp = -np.log10(np.maximum(pvals, 1e-300))
        ax.scatter(positions, logp, c='black', s=5)
        ax.set_title(f'Манхэттен: {title_ru}')
        add_significance_lines(ax, pvals, log_print=log_print, test_name=test_name)
    plt.tight_layout()
    plt.savefig(output_dir / 'manhattan_metrics.png', dpi=150)
    plt.close()
    log_print("Манхэттены метрик сохранены.")
    
    # ----- Манхэттены покрытия -----
    for cov, name in [(cov_raw, 'Невзвешенное покрытие'), (cov_weighted, 'Взвешенное покрытие')]:
        vals = cov[1:]
        med = np.median(vals)
        mad_val = np.median(np.abs(vals - med))
        if mad_val == 0:
            std = np.std(vals)
            z = (vals - med) / std
        else:
            z = 0.6745 * (vals - med) / mad_val
        pvals = 2 * norm.sf(np.abs(z))
        fig, ax = plt.subplots(figsize=(14,5))
        ax.scatter(range(1, MT_LEN+1), -np.log10(np.maximum(pvals, 1e-300)), c='steelblue', s=2)
        ax.set_title(f'Манхэттен: {name}')
        add_significance_lines(ax, pvals, log_print=log_print, test_name='MAD-Z')
        plt.tight_layout()
        plt.savefig(output_dir / f'manhattan_{name.replace(" ","_")}.png', dpi=150)
        plt.close()
    
    # Сглаженное покрытие
    fig, ax = plt.subplots(figsize=(14,4))
    ax.plot(range(1, MT_LEN+1), cov_weighted[1:], linewidth=0.3, color='darkblue', alpha=0.5, label='Исходное')
    smooth = np.convolve(cov_weighted[1:], np.ones(500)/500, mode='same')
    ax.plot(range(1, MT_LEN+1), smooth, linewidth=1.5, color='red', label='Сглаженное (окно 500)')
    ax.legend()
    ax.set_title('Взвешенное покрытие вдоль мтДНК')
    ax.set_xlabel('Позиция в геноме (п.н.)')
    ax.set_ylabel('Сумма эффективных длин повторов')
    plt.tight_layout()
    plt.savefig(output_dir / 'coverage_plot.png', dpi=150)
    plt.close()
    
    # Гистограммы
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,5))
    vals = cov_weighted[1:]
    ax1.hist(vals, bins=50, density=True, color='lightblue', edgecolor='black', alpha=0.7)
    ax1.set_title('Взвешенное покрытие')
    ax1.text(0.95,0.95, f'Среднее={np.mean(vals):.1f}\nМедиана={np.median(vals):.1f}\nАсимметрия={stats.skew(vals):.2f}\nЭксцесс={stats.kurtosis(vals):.2f}',
             transform=ax1.transAxes, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.8))
    all_cnt = [metrics[pos][nuc]['total'] for pos in per_pos for nuc in per_pos[pos]]
    ax2.hist(all_cnt, bins=50, density=True, color='lightgreen', edgecolor='black', alpha=0.7)
    ax2.set_title('Число повторов в файле')
    ax2.text(0.95,0.95, f'Среднее={np.mean(all_cnt):.1f}\nМедиана={np.median(all_cnt):.1f}\nАсимметрия={stats.skew(all_cnt):.2f}\nЭксцесс={stats.kurtosis(all_cnt):.2f}',
             transform=ax2.transAxes, va='top', ha='right', bbox=dict(facecolor='white', alpha=0.8))
    plt.tight_layout()
    plt.savefig(output_dir / 'histograms.png', dpi=150)
    plt.close()
    
    # Гистограммы максимальных длин
    max_perf_lengths = []
    max_all_lengths = []
    for pos, nucs in per_pos.items():
        for nuc, lines in nucs.items():
            ml = max_eff_len(lines, perfect_only=True)
            if ml > 0: max_perf_lengths.append(ml)
            ml_all = max_eff_len(lines, perfect_only=False)
            if ml_all > 0: max_all_lengths.append(ml_all)
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12,5))
    ax1.hist(max_perf_lengths, bins=30, density=True, color='lightcoral', edgecolor='black', alpha=0.7)
    ax1.set_title('Максимальная длина совершенных повторов')
    ax1.set_xlabel('Длина (п.н.)')
    ax1.set_ylabel('Плотность')
    ax2.hist(max_all_lengths, bins=30, density=True, color='lightgreen', edgecolor='black', alpha=0.7)
    ax2.set_title('Максимальная длина всех повторов')
    ax2.set_xlabel('Длина (п.н.)')
    ax2.set_ylabel('Плотность')
    plt.tight_layout()
    plt.savefig(output_dir / 'histograms_max_dlin.png', dpi=150)
    plt.close()
    
    # Диаграммы рассеяния
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
    
    def scatter_with_regression(ax, x, y, color, xlabel, ylabel, title):
        if len(x)==0: return
        x_jit = x + np.random.normal(0, 0.5, size=len(x))
        y_jit = y + np.random.normal(0, 0.5, size=len(y))
        ax.scatter(x_jit, y_jit, color=color, alpha=0.3, s=8)
        if len(x)>1:
            slope, intercept, r, p, _ = stats.linregress(x, y)
            line_x = np.array([x.min(), x.max()])
            ax.plot(line_x, slope*line_x+intercept, 'k--')
            ax.text(0.02,0.98, f'R²={r**2:.3f}\np={p:.2e}', transform=ax.transAxes, va='top',
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=7)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.axhline(0, color='gray', linestyle='--')
    
    x_p = np.array([d[0] for d in scatter_data])
    y_p = np.array([d[1] for d in scatter_data])
    fig, ax = plt.subplots(figsize=(8,6))
    scatter_with_regression(ax, x_p, y_p, 'blue',
                            'Исходная длина совершенного повтора у референса (п.н.)',
                            'Изменение длины (альтернатива – референс) (п.н.)',
                            f'Совершенные повторы: зависимость изменения длины от исходной длины – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'scatter_delta_perfect.png', dpi=150)
    plt.close()
    
    x_i = np.array([d[2] for d in scatter_data])
    y_i = np.array([d[3] for d in scatter_data])
    fig, ax = plt.subplots(figsize=(8,6))
    scatter_with_regression(ax, x_i, y_i, 'red',
                            'Исходная длина любого повтора у референса (п.н.)',
                            'Изменение длины (альтернатива – референс) (п.н.)',
                            f'Все повторы: зависимость изменения длины от исходной длины – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'scatter_delta_imperfect.png', dpi=150)
    plt.close()
    
    x_ref = np.array([d[0] for d in scatter_data])
    y_alt = np.array([d[4] for d in scatter_data])
    fig, ax = plt.subplots(figsize=(8,6))
    scatter_with_regression(ax, x_ref, y_alt, 'green',
                            'Длина совершенного повтора у референса (п.н.)',
                            'Длина совершенного повтора у альтернативы (п.н.)',
                            f'Сравнение длин совершенных повторов – {genome_name}')
    ax.plot([0, max(x_ref)], [0, max(x_ref)], 'k--', lw=0.5, label='Линия равенства')
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / 'scatter_alt_vs_ref_perfect.png', dpi=150)
    plt.close()
    
    # Боксплоты 99-го процентиля
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
                ti_perf.append(abs(ap-rp))
                ti_imperf.append(abs(ai-ri))
            else:
                tv_perf.append(abs(ap-rp))
                tv_imperf.append(abs(ai-ri))
    
    def plot_percentile_99(ax, data, labels, title):
        vals = [np.percentile(d, 99) if d else np.nan for d in data]
        ax.bar(labels, vals, color=['#2c7bb6', '#d7191c'])
        for i,v in enumerate(vals):
            if not np.isnan(v): ax.text(i, v, f'{v:.1f}', ha='center', va='bottom')
        ax.set_title(title)
        ax.set_ylabel('99-й процентиль |Δ| длины (п.н.)')
    
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    plot_percentile_99(ax1, [ti_perf, tv_perf], ['Транзиции','Трансверсии'],
                       f'Совершенные повторы – 99-й процентиль |Δ| – {genome_name}')
    plot_percentile_99(ax2, [ti_imperf, tv_imperf], ['Транзиции','Трансверсии'],
                       f'Все повторы – 99-й процентиль |Δ| – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'boxplot_99percentile_delta.png', dpi=150)
    plt.close()
    
    # GC и H-связи (все повторы)
    gc_ti, gc_tv, hb_ti, hb_tv = [], [], [], []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref or ref not in per_pos[pos]: continue
        ref_lines = per_pos[pos][ref]
        rgc = np.mean([float(l.split('\t')[10]) for l in ref_lines if len(l.split('\t')) >= 14])
        rhb = np.mean([float(l.split('\t')[12]) for l in ref_lines if len(l.split('\t')) >= 14])
        for nuc in 'ACGT':
            if nuc == ref or nuc not in per_pos[pos]: continue
            alt_lines = per_pos[pos][nuc]
            agc = np.mean([float(l.split('\t')[10]) for l in alt_lines if len(l.split('\t')) >= 14])
            ahb = np.mean([float(l.split('\t')[12]) for l in alt_lines if len(l.split('\t')) >= 14])
            if np.isnan(agc) or np.isnan(rgc): continue
            is_ti = {ref, nuc} in ({'A','G'}, {'C','T'})
            if is_ti:
                gc_ti.append(abs(agc - rgc))
                hb_ti.append(abs(ahb - rhb))
            else:
                gc_tv.append(abs(agc - rgc))
                hb_tv.append(abs(ahb - rhb))
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    plot_percentile_99(ax1, [gc_ti, gc_tv], ['Транзиции','Трансверсии'],
                       f'Изменение GC-состава (|ΔGC|) – {genome_name}')
    plot_percentile_99(ax2, [hb_ti, hb_tv], ['Транзиции','Трансверсии'],
                       f'Изменение числа водородных связей (|ΔH|) – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'boxplot_99percentile_gc_hb.png', dpi=150)
    plt.close()
    
    # GC и H-связи для самых длинных повторов
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
                try:
                    eff = int(p[3]) - int(p[9])
                except: continue
                if eff == max_eff_alt and len(p) >= 14:
                    alt_gc_long.append(float(p[10]))
                    alt_hb_long.append(float(p[12]))
            ref_gc_long = []; ref_hb_long = []
            for l in ref_lines:
                p = l.split('\t')
                try:
                    eff = int(p[3]) - int(p[9])
                except: continue
                if eff == max_eff_ref and len(p) >= 14:
                    ref_gc_long.append(float(p[10]))
                    ref_hb_long.append(float(p[12]))
            if alt_gc_long and ref_gc_long:
                diff_gc = abs(np.mean(alt_gc_long) - np.mean(ref_gc_long))
                diff_hb = abs(np.mean(alt_hb_long) - np.mean(ref_hb_long))
                is_ti = {ref, nuc} in ({'A','G'}, {'C','T'})
                if is_ti:
                    gc_long_ti.append(diff_gc)
                    hb_long_ti.append(diff_hb)
                else:
                    gc_long_tv.append(diff_gc)
                    hb_long_tv.append(diff_hb)
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
    plot_percentile_99(ax1, [gc_long_ti, gc_long_tv], ['Транзиции','Трансверсии'],
                       f'Изменение GC-состава для самых длинных повторов – {genome_name}')
    plot_percentile_99(ax2, [hb_long_ti, hb_long_tv], ['Транзиции','Трансверсии'],
                       f'Изменение числа H-связей для самых длинных повторов – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'boxplot_99percentile_gc_hb_longest.png', dpi=150)
    plt.close()
    
    # Тепловые карты (изменение длины)
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
    ax.set_xlabel('Альтернативный нуклеотид')
    ax.set_ylabel('Референсный нуклеотид')
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Среднее изменение длины (п.н.)')
    for i in range(4):
        for j in range(4):
            if (nuc_order[i], nuc_order[j]) in transitions_max:
                txt = f'{hm[i,j]:.2f}'
                if sig[i,j]: txt += '*'
                ax.text(j, i, txt, ha='center', va='center', fontsize=9,
                        color='white' if abs(hm[i,j])>np.max(np.abs(hm))/2 else 'black')
    ax.set_title(f'Изменение максимальной длины совершенных повторов – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'heatmap_max_perfect.png', dpi=150)
    plt.close()
    
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
    ax.set_xlabel('Альтернативный нуклеотид')
    ax.set_ylabel('Референсный нуклеотид')
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Среднее изменение средней длины (п.н.)')
    for i in range(4):
        for j in range(4):
            if (nuc_order[i], nuc_order[j]) in transitions_all:
                txt = f'{hm2[i,j]:.2f}'
                if sig2[i,j]: txt += '*'
                ax.text(j, i, txt, ha='center', va='center', fontsize=9,
                        color='white' if abs(hm2[i,j])>np.max(np.abs(hm2))/2 else 'black')
    ax.set_title(f'Изменение средней эффективной длины всех повторов – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'heatmap_all_repeats.png', dpi=150)
    plt.close()
    
    # Инфографика GC/H
    rgc_m, rgc_r, agc_m, agc_r = [], [], [], []
    rhb_m, rhb_r, ahb_m, ahb_r = [], [], [], []
    for pos in sorted(per_pos.keys()):
        ref = ref_nuc.get(pos)
        if not ref: continue
        for nuc, lines in per_pos[pos].items():
            is_ref = (nuc == ref)
            for l in lines:
                p = l.split('\t')
                if len(p) < 14: continue
                try:
                    mgc = float(p[10])
                    rgc = float(p[11])
                    mhb = int(p[12])
                    rhb = int(p[13])
                except: continue
                if is_ref:
                    rgc_m.append(mgc); rgc_r.append(rgc); rhb_m.append(mhb); rhb_r.append(rhb)
                else:
                    agc_m.append(mgc); agc_r.append(rgc); ahb_m.append(mhb); ahb_r.append(rhb)
    def safe_mean_std(lst): return (np.mean(lst), np.std(lst)) if lst else (np.nan, np.nan)
    fig, axes = plt.subplots(2,2, figsize=(10,8))
    m1,s1 = safe_mean_std(rgc_m); m2,s2 = safe_mean_std(agc_m)
    axes[0,0].bar(['Референс','Альтернатива'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[0,0].set_title('GC-состав мотива')
    axes[0,0].set_ylabel('Доля GC')
    m1,s1 = safe_mean_std(rgc_r); m2,s2 = safe_mean_std(agc_r)
    axes[0,1].bar(['Референс','Альтернатива'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[0,1].set_title('GC-состав повтора')
    axes[0,1].set_ylabel('Доля GC')
    m1,s1 = safe_mean_std(rhb_m); m2,s2 = safe_mean_std(ahb_m)
    axes[1,0].bar(['Референс','Альтернатива'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[1,0].set_title('Число H-связей в мотиве')
    axes[1,0].set_ylabel('Число связей')
    m1,s1 = safe_mean_std(rhb_r); m2,s2 = safe_mean_std(ahb_r)
    axes[1,1].bar(['Референс','Альтернатива'],[m1,m2], yerr=[s1,s2], color=['gray','orange'])
    axes[1,1].set_title('Число H-связей в повторе')
    axes[1,1].set_ylabel('Число связей')
    plt.suptitle(f'Сравнение GC и H-связей – {genome_name}')
    plt.tight_layout()
    plt.savefig(output_dir / 'gc_hbonds_infographic.png', dpi=150)
    plt.close()
    
    # Расстояние до большой дуги
    distances = []
    met_vals = []
    for pos in per_pos.keys():
        if not ref_nuc.get(pos): continue
        dist_start = min(abs(pos - MAJOR_ARC_START), MT_LEN - abs(pos - MAJOR_ARC_START))
        dist_end = min(abs(pos - MAJOR_ARC_END), MT_LEN - abs(pos - MAJOR_ARC_END))
        dist = min(dist_start, dist_end)
        ref_n = ref_nuc[pos]
        if ref_n in per_pos[pos]:
            val = metrics[pos][ref_n]['count_ge7']
            distances.append(dist)
            met_vals.append(val)
    if distances:
        bins = np.arange(0, max(distances)+100, 100)
        bin_centers = (bins[:-1] + bins[1:])/2
        binned_mean = [np.mean(np.array(met_vals)[(np.array(distances)>=bins[i]) & (np.array(distances)<bins[i+1])]) for i in range(len(bins)-1)]
        fig, ax = plt.subplots(figsize=(10,4))
        ax.plot(bin_centers, binned_mean, 'o-')
        ax.set_xlabel('Расстояние до ближайшей границы большой дуги (п.н.)')
        ax.set_ylabel('Среднее число повторов с эффективной длиной ≥7 (референс)')
        ax.set_title(f'Зависимость повторяемости от расстояния до границ большой дуги – {genome_name}')
        plt.tight_layout()
        plt.savefig(output_dir / 'bolshaya_duga_rasstoyanie.png', dpi=150)
        plt.close()
    else:
        log_print("Нет данных для построения графика расстояния до большой дуги.")
    
    # Оконный анализ
    if window_size > 0 and STATSMODELS_AVAILABLE:
        log_print("\n=== Оконный анализ ===")
        for metric in ['total', 'perfect', 'count_ge7', 'perfect_ge7', 'imperfect_ge10']:
            fig, _ = windowed_analysis(per_pos, ref_nuc, metrics, metric, MT_LEN,
                                       window_size=window_size, log_print=log_print)
            fig.savefig(output_dir / f'windowed_barplot_{metric}.png', dpi=150)
            plt.close(fig)
    
    # Лесной график для SNV
    if snv_file:
        log_print(f"\n=== Лесной график для мутаций из {snv_file} ===")
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
                    snv_data.append((f"{pos}{ref_allele}→{alt_allele}", ap - rp))
        if snv_data:
            snv_data.sort(key=lambda x: x[1])
            for top_k in [30, 50]:
                half = top_k // 2
                neg = snv_data[:half]
                pos_sel = snv_data[-half:]
                selected = neg + pos_sel
                selected.sort(key=lambda x: x[1])
                labels, deltas = zip(*selected)
                fig, ax = plt.subplots(figsize=(8, max(4, len(labels)*0.3)))
                bars = ax.barh(range(len(labels)), deltas, 
                               color=[COLOR_POS_BAR if d>0 else COLOR_NEG_BAR for d in deltas], 
                               height=0.6)
                ax.set_yticks(range(len(labels)))
                ax.set_yticklabels(labels, fontsize=7)
                ax.axvline(0, color='black', linestyle='--')
                ax.set_xlabel('Изменение максимальной длины совершенного повтора (п.н.)')
                ax.set_title(f'Лесной график SNV (топ {half}) – {genome_name}')
                for i, (bar, d) in enumerate(zip(bars, deltas)):
                    xpos = d + (1 if d>0 else -1)
                    ha = 'left' if d>0 else 'right'
                    ax.text(xpos, i, f'{d}', ha=ha, va='center', fontsize=6)
                plt.tight_layout()
                plt.savefig(output_dir / f'lesnoj_grafik_snv_top{top_k}.png', dpi=150)
                plt.close()
            log_print(f"Лесные графики для топ-30 и топ-50 сохранены.")
        else:
            log_print("Нет данных для лесного графика.")
    
    log_print(f"Анализ генома {genome_name} завершён. Результаты в {output_dir}")
    return per_pos, ref_nuc, metrics, MT_LEN

# ------------------------------------------------------------
# Сравнение двух геномов (только сравнительные графики)
# ------------------------------------------------------------
def compare_two_genomes(input_dir1, input_dir2, ref_seq1, ref_seq2, name1, name2, output_dir, log_print, snv_file=None):
    MT_LEN1 = len(ref_seq1)
    MT_LEN2 = len(ref_seq2)
    log_print(f"Длина генома {name1}: {MT_LEN1} п.н., {name2}: {MT_LEN2} п.н.")
    
    # Загружаем данные из кэша (предполагается, что кэш уже создан при анализе)
    cache1 = output_dir / name1 / 'cleaned_data.pkl.gz'
    cache2 = output_dir / name2 / 'cleaned_data.pkl.gz'
    per_pos1, ref_nuc1, cov_raw1, cov_weighted1, metrics1 = load_cleaned_files_cached(
        input_dir1, ref_seq1, log_print, MT_LEN1, cache1)
    per_pos2, ref_nuc2, cov_raw2, cov_weighted2, metrics2 = load_cleaned_files_cached(
        input_dir2, ref_seq2, log_print, MT_LEN2, cache2)
    
    # ----- Плотность повторов (4 варианта) -----
    def repeat_density(per_pos, genome_len):
        total = sum(len(lines) for pos_dict in per_pos.values() for lines in pos_dict.values())
        return total / (genome_len / 1000.0)
    def repeat_density_perfect_all(per_pos, metrics, genome_len):
        total = 0
        for pos, nucs in per_pos.items():
            for nuc in nucs:
                total += metrics[pos][nuc].get('perfect', 0)
        return total / (genome_len / 1000.0)
    def repeat_density_perfect_ge7(per_pos, metrics, genome_len):
        total = 0
        for pos, nucs in per_pos.items():
            for nuc in nucs:
                total += metrics[pos][nuc].get('perfect_ge7', 0)
        return total / (genome_len / 1000.0)
    def repeat_density_imperfect_ge10(per_pos, metrics, genome_len):
        total = 0
        for pos, nucs in per_pos.items():
            for nuc in nucs:
                total += metrics[pos][nuc].get('imperfect_ge10', 0)
        return total / (genome_len / 1000.0)
    
    d_total1 = repeat_density(per_pos1, MT_LEN1)
    d_total2 = repeat_density(per_pos2, MT_LEN2)
    d_perf_all1 = repeat_density_perfect_all(per_pos1, metrics1, MT_LEN1)
    d_perf_all2 = repeat_density_perfect_all(per_pos2, metrics2, MT_LEN2)
    d_perf_ge7_1 = repeat_density_perfect_ge7(per_pos1, metrics1, MT_LEN1)
    d_perf_ge7_2 = repeat_density_perfect_ge7(per_pos2, metrics2, MT_LEN2)
    d_imperf_ge10_1 = repeat_density_imperfect_ge10(per_pos1, metrics1, MT_LEN1)
    d_imperf_ge10_2 = repeat_density_imperfect_ge10(per_pos2, metrics2, MT_LEN2)
    
    density_data = [
        (d_total1, d_total2, 'все повторы'),
        (d_perf_all1, d_perf_all2, 'совершенные повторы (любой длины)'),
        (d_perf_ge7_1, d_perf_ge7_2, 'совершенные повторы ≥7 п.н.'),
        (d_imperf_ge10_1, d_imperf_ge10_2, 'несовершенные повторы с эффективной длиной ≥10 п.н.')
    ]
    for d1, d2, label in density_data:
        fig, ax = plt.subplots(figsize=(6,5))
        ax.bar([name1, name2], [d1, d2], color=['steelblue', 'orange'])
        ax.set_ylabel('Число повторов на 1000 п.н.')
        ax.set_title(f'Плотность повторов: {label}')
        for i, v in enumerate([d1, d2]):
            ax.text(i, v + 0.05, f'{v:.2f}', ha='center', fontweight='bold')
        plt.tight_layout()
        safe_label = label.replace(' ', '_').replace('≥', 'ge').replace('(', '').replace(')', '').replace('п.н.', '')
        plt.savefig(output_dir / f'plotnost_{safe_label}.png', dpi=150)
        plt.close()
        log_print(f"Плотность {label}: {name1} = {d1:.2f}, {name2} = {d2:.2f}")
    
    # ----- Тепловые карты количества повторов (отдельно для каждого генома) -----
    for metric, metric_ru in [('perfect_ge7', 'совершенные повторы ≥7'), ('imperfect_ge10', 'несовершенные повторы с эффективной длиной ≥10')]:
        for (per_pos, ref_nuc, metrics, name) in [(per_pos1, ref_nuc1, metrics1, name1), (per_pos2, ref_nuc2, metrics2, name2)]:
            trans = defaultdict(list)
            for pos in per_pos:
                ref = ref_nuc.get(pos)
                if not ref or ref not in per_pos[pos]: continue
                ref_val = metrics[pos][ref].get(metric, 0)
                for nuc in 'ACGT':
                    if nuc == ref or nuc not in per_pos[pos]: continue
                    alt_val = metrics[pos][nuc].get(metric, 0)
                    trans[(ref, nuc)].append(alt_val - ref_val)
            nuc_order = ['A','C','G','T']
            hm = np.zeros((4,4)); pm = np.ones((4,4))
            for i,ref in enumerate(nuc_order):
                for j,alt in enumerate(nuc_order):
                    if (ref,alt) in trans and len(trans[(ref,alt)]) >= 2:
                        d = trans[(ref,alt)]
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
            ax.set_xlabel('Альтернативный нуклеотид')
            ax.set_ylabel('Референсный нуклеотид')
            cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label('Среднее изменение числа повторов')
            for i in range(4):
                for j in range(4):
                    if (nuc_order[i], nuc_order[j]) in trans:
                        txt = f'{hm[i,j]:.2f}'
                        if sig[i,j]: txt += '*'
                        ax.text(j, i, txt, ha='center', va='center', fontsize=9,
                                color='white' if abs(hm[i,j])>np.max(np.abs(hm))/2 else 'black')
            ax.set_title(f'Изменение числа {metric_ru} – {name}')
            plt.tight_layout()
            plt.savefig(output_dir / f'teplovaya_karta_kolichestva_{metric}_{name}.png', dpi=150)
            plt.close()
    
    # ----- Большая дуга: разность плотности (конец – начало) для каждого генома -----
    start_win = (MAJOR_ARC_START, MAJOR_ARC_START + 2000 - 1)   # (5747, 7746)
    end_win = (13978, 15977)
    def repeats_in_interval(per_pos, start, end, genome_len):
        count = 0
        if start <= end:
            for pos, nucs in per_pos.items():
                if start <= pos <= end:
                    for lines in nucs.values():
                        count += len(lines)
        else:
            for pos, nucs in per_pos.items():
                if pos >= start or pos <= end:
                    for lines in nucs.values():
                        count += len(lines)
        return count
    rep_start1 = repeats_in_interval(per_pos1, start_win[0], start_win[1], MT_LEN1)
    rep_end1   = repeats_in_interval(per_pos1, end_win[0], end_win[1], MT_LEN1)
    rep_start2 = repeats_in_interval(per_pos2, start_win[0], start_win[1], MT_LEN2)
    rep_end2   = repeats_in_interval(per_pos2, end_win[0], end_win[1], MT_LEN2)
    norm_start1 = rep_start1 / 2.0
    norm_end1   = rep_end1 / 2.0
    norm_start2 = rep_start2 / 2.0
    norm_end2   = rep_end2 / 2.0
    diff1 = norm_end1 - norm_start1
    diff2 = norm_end2 - norm_start2
    
    fig, ax = plt.subplots(figsize=(5,6))
    ax.bar([name1, name2], [diff1, diff2], color=['steelblue', 'orange'])
    ax.axhline(0, color='black', linestyle='--', linewidth=0.8)
    ax.set_ylabel('Разность плотности (конец – начало), на 1000 п.н.')
    ax.set_title('Большая дуга: разность плотности повторов между концом и началом')
    for i, v in enumerate([diff1, diff2]):
        ax.text(i, v + 0.05*(1 if v>=0 else -0.5), f'{v:.2f}', ha='center', fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / 'bolshaya_duga_raznost.png', dpi=150)
    plt.close()
    log_print(f"Разность плотности (конец-начало) для {name1}: {diff1:.2f}, для {name2}: {diff2:.2f}")
    
    # ----- Скользящие окна 10% длины генома -----
    window_ratio = 0.10
    step_ratio = 0.05
    def sliding_max_perfect(per_pos, ref_nuc, genome_len, window_ratio, step_ratio):
        window_size = int(genome_len * window_ratio)
        step = int(genome_len * step_ratio)
        starts = range(1, genome_len - window_size + 2, step)
        max_perfect = []
        rel_positions = []
        for start in starts:
            end = start + window_size - 1
            max_len = 0
            for pos, nucs in per_pos.items():
                if start <= pos <= end:
                    for nuc, lines in nucs.items():
                        if nuc == ref_nuc.get(pos):
                            ml = max_eff_len(lines, perfect_only=True)
                            if ml > max_len:
                                max_len = ml
            max_perfect.append(max_len)
            rel_positions.append((start + end) / 2 / genome_len)
        return rel_positions, max_perfect
    rel1, max1 = sliding_max_perfect(per_pos1, ref_nuc1, MT_LEN1, window_ratio, step_ratio)
    rel2, max2 = sliding_max_perfect(per_pos2, ref_nuc2, MT_LEN2, window_ratio, step_ratio)
    fig, ax = plt.subplots(figsize=(12,5))
    ax.plot(rel1, max1, 'o-', label=name1, markersize=3, linewidth=1)
    ax.plot(rel2, max2, 's-', label=name2, markersize=3, linewidth=1)
    ax.set_xlabel('Относительная позиция вдоль генома (0 – начало, 1 – конец)')
    ax.set_ylabel('Максимальная длина совершенного повтора в окне (п.н.)')
    ax.set_title('Скользящие окна 10% длины генома (шаг 5%)')
    ax.legend()
    plt.tight_layout()
    plt.savefig(output_dir / 'skolzyaschie_okna_10procentov.png', dpi=150)
    plt.close()
    
    # ----- Лесной график SNV для первого генома (если задан) -----
    if snv_file:
        log_print(f"\n=== Лесной график SNV для {name1} ===")
        snv_data = []
        with open(snv_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                pos = int(row['position'])
                ref_allele = row['ref_allele'].upper()
                alt_allele = row['alt_allele'].upper()
                if pos in per_pos1 and ref_allele in per_pos1[pos] and alt_allele in per_pos1[pos]:
                    ref_lines = per_pos1[pos][ref_allele]
                    alt_lines = per_pos1[pos][alt_allele]
                    rp = max_eff_len(ref_lines, True)
                    ap = max_eff_len(alt_lines, True)
                    snv_data.append((f"{pos}{ref_allele}→{alt_allele}", ap - rp))
        if snv_data:
            snv_data.sort(key=lambda x: x[1])
            for top_k in [30, 50]:
                half = top_k // 2
                neg = snv_data[:half]
                pos_sel = snv_data[-half:]
                selected = neg + pos_sel
                selected.sort(key=lambda x: x[1])
                labels, deltas = zip(*selected)
                fig, ax = plt.subplots(figsize=(8, max(4, len(labels)*0.3)))
                bars = ax.barh(range(len(labels)), deltas, 
                               color=[COLOR_POS_BAR if d>0 else COLOR_NEG_BAR for d in deltas], 
                               height=0.6)
                ax.set_yticks(range(len(labels)))
                ax.set_yticklabels(labels, fontsize=7)
                ax.axvline(0, color='black', linestyle='--')
                ax.set_xlabel('Изменение максимальной длины совершенного повтора (п.н.)')
                ax.set_title(f'Лесной график SNV (топ {half}) – {name1}')
                for i, (bar, d) in enumerate(zip(bars, deltas)):
                    xpos = d + (1 if d>0 else -1)
                    ha = 'left' if d>0 else 'right'
                    ax.text(xpos, i, f'{d}', ha=ha, va='center', fontsize=6)
                plt.tight_layout()
                plt.savefig(output_dir / f'lesnoj_grafik_snv_top{top_k}_{name1}.png', dpi=150)
                plt.close()
            log_print(f"Лесные графики для топ-30 и топ-50 сохранены.")
        else:
            log_print("Нет данных для лесного графика.")

# ------------------------------------------------------------
# Режим сравнения (точка входа)
# ------------------------------------------------------------
def run_comparison(input_dir1, ref_path1, input_dir2, ref_path2, output_dir, workers, window_size, snv_file, name1, name2, compare_only=False):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "sravnenie.log"
    with open(log_path, 'w', encoding='utf-8') as logf:
        def log_print(msg):
            print(msg)
            print(msg, file=logf)
        
        # Анализ первого генома (если не указан compare_only или нет кэша)
        out1 = output_dir / name1
        cache1 = out1 / 'cleaned_data.pkl.gz'
        if compare_only:
            log_print("Режим compare_only: пропускаем анализ первого генома, используем существующий кэш.")
        else:
            if cache1.exists() and (out1 / 'stats_run.log').exists():
                log_print(f"Кэш и результаты для {name1} уже существуют, пропускаем анализ.")
            else:
                analyze_single_genome(input_dir1, ref_path1, out1, workers, window_size, snv_file, name1, log_print)
        
        # Анализ второго генома (аналогично)
        out2 = output_dir / name2
        cache2 = out2 / 'cleaned_data.pkl.gz'
        if compare_only:
            log_print("Режим compare_only: пропускаем анализ второго генома, используем существующий кэш.")
        else:
            if cache2.exists() and (out2 / 'stats_run.log').exists():
                log_print(f"Кэш и результаты для {name2} уже существуют, пропускаем анализ.")
            else:
                analyze_single_genome(input_dir2, ref_path2, out2, workers, window_size, None, name2, log_print)
        
        # Сравнительные графики
        log_print("\n=== Построение сравнительных графиков ===")
        ref_seq1 = read_reference(ref_path1)
        ref_seq2 = read_reference(ref_path2)
        compare_two_genomes(
            input_dir1=Path(input_dir1),
            input_dir2=Path(input_dir2),
            ref_seq1=ref_seq1,
            ref_seq2=ref_seq2,
            name1=name1,
            name2=name2,
            output_dir=output_dir,
            log_print=log_print,
            snv_file=snv_file
        )
        log_print("\nСравнение завершено. Результаты:")
        log_print(f"  - Индивидуальная статистика: {out1}, {out2}")
        log_print(f"  - Сравнительные графики: {output_dir}")

# ------------------------------------------------------------
# Оригинальная run_stats (обратная совместимость)
# ------------------------------------------------------------
def run_stats(input_dir, output_dir, ref_path, workers, window_size=0, snv_file=None,
              input_dir2=None, ref_path2=None, name1='Genome1', name2='Genome2'):
    if input_dir2 is not None and ref_path2 is not None:
        # Если заданы два генома, запускаем режим сравнения (без compare_only)
        run_comparison(input_dir, ref_path, input_dir2, ref_path2, output_dir, workers, window_size, snv_file, name1, name2, compare_only=False)
        return
    # Иначе стандартный анализ одного генома
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "stats_run.log"
    with open(log_path, 'w', encoding='utf-8') as logf:
        def log_print(msg):
            print(msg)
            print(msg, file=logf)
        analyze_single_genome(input_dir, ref_path, output_dir, workers, window_size, snv_file, name1, log_print)
