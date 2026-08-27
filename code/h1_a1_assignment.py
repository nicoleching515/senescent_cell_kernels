#!/usr/bin/env python3
"""Phase 9 test A1 — transcript assignment rate on H1, the direct analogue of
`code/assignment_rate.py` (M1: 88.27 % on section 7259).

DECLARED DEVIATION: Phase 7 §12.3 says do NOT download `transcripts.parquet`, and the H1
acquisition did not.  A1 asks for the assignment rate, which cannot be computed from any
other deposited file (there is no `metrics_summary.csv` in this deposit — checked against
the GEO filelist).  Three sections were therefore fetched to
`data/raw_h1_transcripts/` (925 MB) purely for this test; nothing else uses them.

Method is `code/assignment_rate.py` verbatim, run per section.
Writes results/phase9_h1/a1_assignment_rate.csv.
"""
import glob, gzip, os, shutil, sys
import numpy as np, pandas as pd, pyarrow.parquet as pq

SRC = "/workspace/data/raw_h1_transcripts"
OUT = "/workspace/results/phase9_h1/a1_assignment_rate.csv"
rows = []
for gz in sorted(glob.glob(SRC + "/*_transcripts.parquet.gz")):
    sec = os.path.basename(gz).split("_")[1]
    pqf = gz[:-3]
    if not os.path.exists(pqf):
        print("decompressing", os.path.basename(gz), flush=True)
        with gzip.open(gz, "rb") as fi, open(pqf + ".part", "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 24)
        os.rename(pqf + ".part", pqf)
    f = pq.ParquetFile(pqf)
    cols = f.schema_arrow.names
    print(sec, "rows", f.metadata.num_rows, "cols", cols, flush=True)
    cell_col = next((c for c in ("cell_id", "cell_ID") if c in cols), None)
    qv_col = next((c for c in ("qv", "QV") if c in cols), None)
    ov_col = next((c for c in ("overlaps_nucleus",) if c in cols), None)
    assert cell_col, cols
    tot = asg = tot_q = asg_q = nuc_q = 0
    for b in f.iter_batches(batch_size=4_000_000,
                            columns=[c for c in (cell_col, qv_col, ov_col) if c]):
        cid = b.column(cell_col).to_numpy(zero_copy_only=False)
        if cid.dtype.kind in "iu":
            a = cid != -1
        else:
            s = np.asarray(cid, dtype=object).astype(str)
            a = ~np.isin(s, ["UNASSIGNED", "-1", "", "nan"])
        tot += a.size; asg += int(a.sum())
        if qv_col:
            q = b.column(qv_col).to_numpy(zero_copy_only=False) >= 20
            tot_q += int(q.sum()); asg_q += int((a & q).sum())
            if ov_col:
                ov = b.column(ov_col).to_numpy(zero_copy_only=False).astype(bool)
                nuc_q += int((a & q & ov).sum())
    r = dict(section=sec, transcripts_total=tot, assigned=asg,
             assigned_pct=round(100 * asg / tot, 2),
             unassigned_pct=round(100 * (tot - asg) / tot, 2),
             q20_total=tot_q, q20_pct_of_all=round(100 * tot_q / tot, 2) if tot_q else None,
             q20_assigned_pct=round(100 * asg_q / tot_q, 2) if tot_q else None,
             q20_assigned_in_nucleus_pct=round(100 * nuc_q / asg_q, 2) if nuc_q else None,
             verdict="PASS" if (100 * (tot_q - asg_q) / tot_q if tot_q else
                                100 * (tot - asg) / tot) <= 30 else "FAIL")
    rows.append(r); print(r, flush=True)
    os.remove(pqf)          # 1.5-2 GB each; the .gz stays for reproducibility
d = pd.DataFrame(rows); d.to_csv(OUT, index=False)
print("\nM1 reference (7259): 88.27 % assigned, 11.72 % unassigned, PASS at the 30 % threshold")
print(d.to_string()); print("wrote", OUT)
