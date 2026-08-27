#!/usr/bin/env python3
"""Phase 9 test A7 — the negative-control-probe kernel on H1.  MUST BE FLAT.

`code/run_a7_control_probes.py` transplanted to H1: the estimator is
`run_phase3_nulls.SectionFit` / `fit_cell` verbatim (100 um window, 40-point lambda grid,
MIN_RECEIVERS = 2000, nested base -> +N6 -> +N6+N5 designs, 400-replicate spatial block
bootstrap over 100 quantile blocks), with only the response matrix Y substituted.  The H1
cache is reached through `h1_sec`.

PER FEATURE FAMILY, SEPARATELY -- this is not cosmetic.  On M1 the pooled `all_controls`
response was NOT flat (-0.070 SD, p = 0.023) while the 40 negative control probes ALONE were
(-0.0225, p = 0.129), and `PREREG_PHASE8_genesets.md` §11 designates
`E_negative_control_probes` the PRE-REGISTERED PRIMARY technical null.  So the primary A7
result is the `neg_control_probe` row; the pooled row is reported beside it and must never be
called "negative control probes" (PREREG §10.1).

  neg_control_probe     control_probe_counts       40 Negative Control Probes  <- PRIMARY
  neg_control_codeword  control_codeword_counts    609 Negative Control Codewords
  genomic_control       genomic_control_counts     21 Genomic Controls
  all_controls          the sum of the three       "pooled negative-control features"
  neg_probe_rate        probe / (transcripts + 1)  NOT a clean null (its denominator is an
                                                   N5 column) -- reported, never headlined

Usage: python3 code/h1_a7_controls.py [--calls tierA_p95,cdkn1a_pos] [--n-jobs N]
Writes results/phase9_h1/a7_control_probe_{fits,provenance,curves}.csv
"""
import os, sys, time, argparse
import numpy as np, pandas as pd
sys.path.insert(0, "/workspace/code")
import h1_sec                      # noqa: F401  repoint the cache
from joblib import Parallel, delayed
import sasp_phase3 as P
import phase3_core as C
import run_phase3_nulls as R
import h1_common as H

RESPONSES = ["neg_control_probe", "neg_control_codeword", "genomic_control",
             "all_controls", "neg_probe_rate"]
COLS = {"neg_control_probe": "control_probe_counts",
        "neg_control_codeword": "control_codeword_counts",
        "genomic_control": "genomic_control_counts"}
BINS = np.arange(0.0, R.WINDOW_UM + 1e-9, 5.0)


def control_matrix(sec):
    c = H.cells_table(sec.name).set_index("cell_id")
    c = c.reindex(pd.Index(sec.cell_id.astype(str)))
    assert c[COLS["neg_control_probe"]].notna().all(), "cache cell_id not in cells.parquet"
    raw = {k: c[v].to_numpy(float) for k, v in COLS.items()}
    raw["all_controls"] = sum(raw[k] for k in COLS)
    raw["neg_probe_rate"] = raw["neg_control_probe"] / (c["transcript_counts"].to_numpy(float) + 1.0)
    Y = np.column_stack([raw[k] for k in RESPONSES])
    prov = [dict(section=sec.name, response=k, n_cells=len(c),
                 mean_per_cell=float(np.mean(raw[k])), sd_per_cell=float(np.std(raw[k])),
                 frac_cells_nonzero=float(np.mean(raw[k] > 0)),
                 max_per_cell=float(np.max(raw[k]))) for k in RESPONSES]
    sd = Y.std(0); sd[sd < 1e-12] = 1.0
    return (Y - Y.mean(0)) / sd, pd.DataFrame(prov)


def _curves(sf, celltype):
    idx = sf.receivers(celltype)
    if idx.sum() < R.MIN_RECEIVERS:
        return pd.DataFrame()
    d = sf.d_obs[idx]
    k = np.clip(np.digitize(d, BINS) - 1, 0, len(BINS) - 2)
    rows = []
    for j, nm in enumerate(RESPONSES):
        y = sf.Y[idx, j].astype(float)
        X1, _, _, _ = R._designs(sf, idx, j)
        resid = P.FixedLambdaFitter(X1, y[:, None]).Y[:, 0]
        for b in range(len(BINS) - 1):
            m = k == b
            if m.sum() < 50:
                continue
            rows.append(dict(section=sf.sec.name, call=sf.call, response=nm,
                             celltype=celltype, bin_lo=BINS[b], bin_hi=BINS[b + 1],
                             bin_mid=0.5 * (BINS[b] + BINS[b + 1]), n=int(m.sum()),
                             mean_raw=float(y[m].mean()),
                             sem_raw=float(y[m].std(ddof=1) / np.sqrt(m.sum())),
                             mean_resid=float(resid[m].mean()),
                             sem_resid=float(resid[m].std(ddof=1) / np.sqrt(m.sum()))))
    return pd.DataFrame(rows)


def _job(sample, call, seed):
    # loky workers are fresh interpreters; repoint the cache inside the worker too.
    import h1_sec                      # noqa: F401
    t0 = time.time()
    sf = R.SectionFit(sample, call, seed)
    Y, prov = control_matrix(sf.sec)
    sf.Y = Y
    sf.NB = C.neighbour_baseline(sf.sec, Y, sf.sender)
    old = P.MODULES
    P.MODULES = RESPONSES
    try:
        rows = []
        types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
        big = max(types, key=lambda t: int(sf.receivers(t).sum()))
        for t in types:
            if sf.receivers(t).sum() < R.MIN_RECEIVERS:
                continue
            for j in range(len(RESPONSES)):
                rows.append(R.fit_cell(sf, t, j, seed + 7 * j, tag="A7"))
        cur = _curves(sf, big)
    finally:
        P.MODULES = old
    print("[a7] %s %s %d rows %.0fs" % (sample, call, len(rows), time.time() - t0), flush=True)
    return rows, prov, cur


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", default="tierA_p95,cdkn1a_pos")
    ap.add_argument("--sections", default="all")
    ap.add_argument("--n-jobs", type=int, default=7)
    a = ap.parse_args()
    secs = list(H.ALL_SECTIONS) if a.sections == "all" else a.sections.split(",")
    calls = a.calls.split(",")
    jobs = [(s, c, P.MASTER_SEED + 1000 * i + j)
            for i, s in enumerate(secs) for j, c in enumerate(calls)]
    out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
        delayed(_job)(s, c, sd) for s, c, sd in jobs)
    df = pd.DataFrame([r for rs, _, _ in out for r in rs]).rename(columns={"module": "response"})
    df.to_csv(H.RESULTS + "/a7_control_probe_fits.csv", index=False)
    pd.concat([p for _, p, _ in out]).drop_duplicates(["section", "response"]) \
      .to_csv(H.RESULTS + "/a7_control_probe_provenance.csv", index=False)
    cu = [c for _, _, c in out if len(c)]
    if cu:
        pd.concat(cu).to_csv(H.RESULTS + "/a7_control_probe_curves.csv", index=False)
    print(df.shape, "->", H.RESULTS + "/a7_control_probe_fits.csv")


if __name__ == "__main__":
    main()
