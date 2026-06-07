import heapq
import os
import subprocess
import tempfile
from .constants import HAS_SYSTEM_SORT

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
