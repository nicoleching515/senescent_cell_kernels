"""Section 8 Test 1: transcript assignment rate.

The plan flags >30% unassigned as evidence of bleed-through, which manufactures
spatial autocorrelation -- the artefact class this project accuses other methods
of failing to control. Computed on 7259, an ADMISSIBLE section.
"""
import pyarrow.parquet as pq, numpy as np, sys

P = "/workspace/data/raw/7259_liver_sbr_Male_26-U1/transcripts.parquet"
f = pq.ParquetFile(P)
print("columns:", f.schema_arrow.names)
print("rows:", f.metadata.num_rows)

cols = f.schema_arrow.names
cell_col = next((c for c in ("cell_id","cell_ID") if c in cols), None)
qv_col   = next((c for c in ("qv","QV") if c in cols), None)
ov_col   = next((c for c in ("overlaps_nucleus",) if c in cols), None)
assert cell_col, "no cell id column"

tot = asg = tot_q = asg_q = nuc_q = 0
for b in f.iter_batches(batch_size=2_000_000, columns=[c for c in (cell_col,qv_col,ov_col) if c]):
    cid = b.column(cell_col).to_numpy(zero_copy_only=False)
    # Xenium marks unassigned transcripts as "UNASSIGNED" (str) or -1 (int)
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

print(f"\n=== SECTION 8 TEST 1 — ASSIGNMENT RATE (7259, admissible) ===")
print(f"  transcripts total            {tot:,}")
print(f"  assigned to a cell           {asg:,}  ({asg/tot*100:.2f}%)")
print(f"  UNASSIGNED                   {tot-asg:,}  ({(tot-asg)/tot*100:.2f}%)")
if qv_col:
    print(f"  Q>=20 transcripts            {tot_q:,}  ({tot_q/tot*100:.2f}% of all)")
    print(f"  Q>=20 AND assigned           {asg_q:,}  ({asg_q/tot_q*100:.2f}% of Q>=20)")
    print(f"  Q>=20 unassigned             {tot_q-asg_q:,}  ({(tot_q-asg_q)/tot_q*100:.2f}%)")
    if ov_col:
        print(f"  Q>=20 assigned, in nucleus   {nuc_q:,}  ({nuc_q/asg_q*100:.2f}% of assigned)")
thr = (tot_q-asg_q)/tot_q*100 if qv_col else (tot-asg)/tot*100
print(f"\n  PLAN THRESHOLD: >30% unassigned => bleed-through concern")
print(f"  VERDICT: {thr:.2f}% unassigned -> {'FAIL' if thr>30 else 'PASS'}")
