import cProfile
import pstats
import io
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.analysis import get_repeatability_table

# Profile the pure-Python mode using the same synthetic sequence as bench
length = 20000
seq = ('ACGT' * (length // 4 + 1))[:length]
pos = length // 2

print(f"Profiling pure-Python get_repeatability_table on synthetic seq length={length} (this may take ~2 minutes)...")
pr = cProfile.Profile()
pr.enable()
get_repeatability_table(seq, pos, 'A', 'G', max_flank_left=20, max_flank_right=20, min_length=5, max_length=41, use_numpy=False, use_parallel=False)
pr.disable()

s = io.StringIO()
ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
ps.print_stats(20)

print(s.getvalue())
