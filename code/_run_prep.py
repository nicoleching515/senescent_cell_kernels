import sys
sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P
from joblib import Parallel, delayed
r = Parallel(n_jobs=6, prefer="processes", verbose=5)(
    delayed(P.prep)(s, True) for s in P.ALL_SECTIONS)
for x in r:
    print(x, flush=True)
print("PREP DONE", flush=True)
