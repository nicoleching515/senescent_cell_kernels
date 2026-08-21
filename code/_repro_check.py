"""Compare a re-run of the Phase 1 synthetic sweep against the stored results.

The README flagged, honestly, that numpy drifted 1.26.3 -> 2.4.6 partway through
the project and that seeded reproducibility had "not been re-verified end to
end".  This settles it: re-run code/sasp_sweep.py with RESULTS pointed at a fresh
directory, then run this.

Usage:  python3 _repro_check.py OLD_DIR NEW_DIR [OUT_CSV]
"""
import sys, os, glob
import numpy as np, pandas as pd

old_dir, new_dir = sys.argv[1], sys.argv[2]
out = sys.argv[3] if len(sys.argv) > 3 else None

rows, worst = [], 0.0
files = sorted(os.path.basename(f) for f in glob.glob(f"{old_dir}/*.csv"))
missing = [f for f in files if not os.path.exists(f"{new_dir}/{f}")]
extra = sorted(set(os.path.basename(f) for f in glob.glob(f"{new_dir}/*.csv")) - set(files))

for fn in files:
    if fn in missing:
        continue
    a = pd.read_csv(f"{old_dir}/{fn}")
    b = pd.read_csv(f"{new_dir}/{fn}")
    same_shape = a.shape == b.shape
    same_cols = list(a.columns) == list(b.columns)
    # Column ORDER is not part of the result -- these frames are assembled from
    # per-run dicts, so the order can differ between runs.  Compare by NAME and
    # report order separately.
    same_colset = set(a.columns) == set(b.columns)
    num = [c for c in a.columns if pd.api.types.is_numeric_dtype(a[c])]
    md, nexact, ncmp = 0.0, 0, 0
    if same_shape and same_colset:
        for c in num:
            x, y = a[c].to_numpy(float), b[c].to_numpy(float)
            m = np.isfinite(x) & np.isfinite(y)
            if (np.isfinite(x) != np.isfinite(y)).any():
                md = np.inf
            if m.sum():
                d = np.abs(x[m] - y[m])
                sc = np.maximum(np.abs(x[m]), 1e-12)
                md = max(md, float((d / sc).max()))
                ncmp += int(m.sum())
                nexact += int((x[m] == y[m]).sum())
    else:
        md = np.inf
    worst = max(worst, md if np.isfinite(md) else 1e9)
    rows.append(dict(file=fn, rows_old=len(a), rows_new=len(b),
                     same_shape=same_shape, same_columns=same_cols,
                     same_column_set=same_colset,
                     n_values=ncmp, n_bit_identical=nexact,
                     frac_bit_identical=(nexact / ncmp) if ncmp else np.nan,
                     max_rel_diff=md))

R = pd.DataFrame(rows)
print(f"files compared : {len(R)}")
print(f"missing in new : {missing if missing else 'none'}")
print(f"extra in new   : {extra if extra else 'none'}")
if len(R):
    print(f"shape/col match: {int(R.same_shape.sum())}/{len(R)} shape, "
          f"{int(R.same_column_set.sum())}/{len(R)} column sets, "
          f"{int(R.same_columns.sum())}/{len(R)} column ORDER")
    tot, ex = int(R.n_values.sum()), int(R.n_bit_identical.sum())
    print(f"values compared: {tot:,}")
    print(f"bit-identical  : {ex:,}  ({100.0 * ex / tot:.4f}%)" if tot else "")
    print(f"max rel diff   : {R.max_rel_diff.max():.3e}")
    bad = R[R.max_rel_diff > 1e-12].sort_values("max_rel_diff", ascending=False)
    if len(bad):
        print("\nfiles differing above 1e-12 relative:")
        print(bad[["file", "frac_bit_identical", "max_rel_diff"]].to_string(index=False))
    else:
        print("\nEVERY value in every file is bit-identical.")
if out:
    R.to_csv(out, index=False)
    print("\nwrote", out)
