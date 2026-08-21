#!/usr/bin/env python3
"""
Phase 3 T2 — combined estimate under N2 + N5 + N6, with a DONOR-level bootstrap
over SBR animals (Section 24.1).

One section per mouse, four SBR mice.  The donor bootstrap is implemented by
making the *donor* the block of `sasp_estimators.BlockProfiler`: a bootstrap
replicate is then a multiplicity vector over four blocks, which is exactly a
resample of animals.  With four donors there are only 35 distinct resamples, so
the interval is lumpy and the result is labelled a CASE STUDY, per Section 24.1.

Responses are z-scored WITHIN section before pooling, so beta is in
within-section sd units and a section-level batch shift in the module score
cannot masquerade as a kernel amplitude.
"""
from __future__ import annotations

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/code")
import sasp_estimators as E
import sasp_phase3 as P
import run_phase3_nulls as RN

RES = P.RESULTS
N_DONOR_BOOT = 2000


_SF_CACHE = {}


def _sf(sample, call, si):
    k = (sample, call)
    if k not in _SF_CACHE:
        _SF_CACHE[k] = RN.SectionFit(sample, call, P.MASTER_SEED + 31 * si,
                                     types=P.CANON_TYPES_MERGED)
    return _SF_CACHE[k]


def pooled(sections, call, celltype, j, zone=None):
    parts = []
    for si, s in enumerate(sections):
        sf = _sf(s, call, si)
        extra = None
        if zone is not None:
            extra = sf.sec.compartment == zone
        idx = sf.receivers(celltype, extra)
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        y = sf.Y[idx, j].astype(float)
        y = (y - y.mean()) / (y.std() + 1e-12)
        X1, X2, X3, pp = RN._designs(sf, idx, j)
        parts.append(dict(si=si, name=s, y=y, d_s=sf.d_obs[idx],
                          d_d=sf.d_dec[idx], X1=X1, pp=pp, n=int(idx.sum()),
                          smd=sf.match["max_smd_after"]))
    if len(parts) < 2:
        return None
    nsec = len(parts)
    # section dummies replace the single intercept
    p_extra = nsec - 1
    Xs = []
    for a, pr in enumerate(parts):
        D = np.zeros((pr["n"], p_extra))
        if a > 0:
            D[:, a - 1] = 1.0
        Xs.append(np.column_stack([pr["X1"][:, :1], D, pr["X1"][:, 1:]]))
    X = np.vstack(Xs)
    y = np.concatenate([pr["y"] for pr in parts])
    d_s = np.concatenate([pr["d_s"] for pr in parts])
    d_d = np.concatenate([pr["d_d"] for pr in parts])
    bid = np.concatenate([np.full(pr["n"], a) for a, pr in enumerate(parts)])
    pp = parts[0]["pp"]
    pb = 1 + p_extra
    pn6 = pb + 1
    pfull = X.shape[1]
    lam = RN.lam_grid()
    prof = E.BlockProfiler(d_s, d_d, y, X, bid, nsec, lam, lam)
    one = np.ones(nsec)

    base = prof.fit1(one, pb)
    t0 = base["t"]
    full = prof.fit1(one, pfull)
    r2 = prof.fit2_shared(one, pfull)
    out = dict(call=call, celltype=celltype, module=P.MODULES[j],
               zone=zone or "all", n_donors=nsec, n=int(y.size),
               sections=";".join(str(parts[a]["name"]) for a in range(nsec)),
               lam_naive=base["lam"], beta_naive=base["beta"],
               lam_naive_railed=int(t0 in (0, lam.size - 1)),
               lam_full_profiled=full["lam"], beta_full_profiled=full["beta"],
               lam_n2n5n6_profiled=r2["lam"], beta_n2n5n6_profiled=r2["beta"],
               beta_full_fixedlam=prof.beta_at(one, pfull, t0)[0],
               beta_n2n5n6_fixedlam=prof.beta2_at(one, pfull, t0)[0],
               max_smd_after=float(np.max([pr["smd"] for pr in parts])))
    out["sf_n5n6"] = out["beta_full_fixedlam"] / out["beta_naive"]
    out["sf_n2n5n6"] = out["beta_n2n5n6_fixedlam"] / out["beta_naive"]

    rng = np.random.default_rng(P.MASTER_SEED + j)
    keep = {k: [] for k in ("lam_naive", "beta_naive", "lam_full",
                            "beta_full", "beta_n2n5n6", "sf_n2n5n6")}
    for _ in range(N_DONOR_BOOT):
        m = rng.multinomial(nsec, np.full(nsec, 1.0 / nsec)).astype(float)
        try:
            b0 = prof.fit1(m, pb)
            bf = prof.fit1(m, pfull)
            b2 = prof.beta2_at(m, pfull, b0["t"])[0]
        except Exception:
            continue
        keep["lam_naive"].append(b0["lam"]); keep["beta_naive"].append(b0["beta"])
        keep["lam_full"].append(bf["lam"]); keep["beta_full"].append(bf["beta"])
        keep["beta_n2n5n6"].append(b2)
        keep["sf_n2n5n6"].append(b2 / b0["beta"] if b0["beta"] else np.nan)
    for k, v in keep.items():
        v = np.asarray(v, float); v = v[np.isfinite(v)]
        if v.size > 20:
            out[f"{k}_donor_lo"] = float(np.quantile(v, .025))
            out[f"{k}_donor_hi"] = float(np.quantile(v, .975))
            out[f"{k}_donor_med"] = float(np.median(v))
    return out


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--celltypes", default="Hepatocytes,Biliary/ductular,"
                                           "Endothelial,Macrophages,"
                                           "Mesenchymal,T/NK cells")
    ap.add_argument("--sections", default="inband")
    ap.add_argument("--call", default=RN.PRIMARY_CALL)
    a = ap.parse_args()
    SECS = {"inband": P.IN_BAND, "sbr": P.SBR, "all": P.ALL_SECTIONS}[a.sections]
    rows = []
    for ct in a.celltypes.split(","):
        for j in range(len(P.MODULES)):
            r = pooled(SECS, a.call, ct, j)
            if r:
                rows.append(r)
                print(f"{ct:24s} {P.MODULES[j]:22s} lam={r['lam_naive']:6.1f} "
                      f"beta={r['beta_naive']:+.4f} "
                      f"SF(N2+N5+N6)={r['sf_n2n5n6']:+.3f}", flush=True)
    if "Hepatocytes" in a.celltypes:
        for zone in ("periportal", "midzonal", "pericentral"):
            for j in range(len(P.MODULES)):
                r = pooled(SECS, a.call, "Hepatocytes", j, zone=zone)
                if r:
                    rows.append(r)
                    print(f"[{zone:12s}] {P.MODULES[j]:22s} "
                          f"lam={r['lam_naive']:6.1f} beta={r['beta_naive']:+.4f} "
                          f"SF={r['sf_n2n5n6']:+.3f}", flush=True)
    pd.DataFrame(rows).to_csv(f"{RES}/combined_donor.csv", index=False)
    print("wrote combined_donor.csv")
