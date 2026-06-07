import tempfile
import os
from pathlib import Path
from .utils import (
    norm_interval, cluster_by_repeat_overlap,
    external_sort, process_group, gc_content, h_bonds
)

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

