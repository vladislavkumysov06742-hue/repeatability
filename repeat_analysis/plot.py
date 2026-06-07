import matplotlib.pyplot as plt
from pathlib import Path
from .constants import (
    MAJOR_ARC_START, MAJOR_ARC_END, SNP_NAMES, REF_SUFFIX, ALT_SUFFIX,
    PHENOTYPE, COLOR_POS_BAR, COLOR_NEG_BAR, COLOR_PHEN_POS, COLOR_PHEN_NEG
)

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
