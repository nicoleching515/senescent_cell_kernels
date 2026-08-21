#!/usr/bin/env python3
"""
Phase 3 — N8: gene-set disjointness, the scrambled-response control, and a
quantification of the DeepScence/Tier B circularity Bio flagged (BIO_PHASE2 §4.2).

Three products:
  1. `n8_disjointness.csv`  Tier A x Tier B and CoreScence x Tier B overlaps.
  2. `n8_scrambled.csv`     beta for the real module vs beta for 200
     expression-matched random gene sets (Bio's Tier E3), same scorer, same
     cells, same lambda.  Surviving fraction = (b_obs - mean b_rand)/b_obs.
  3. `n8_circularity.csv`   beta with a CoreScence-derived sender call, with and
     without the CoreScence genes stripped out of the response module.  The
     drop is the part of a DeepScence-sender / Tier B-readout fit that is the
     same genes on both sides.

Scoring: mean(set) - mean(expression-matched control), i.e. the same estimator
scanpy's `score_genes` uses (which is how Bio built modules_*.csv), implemented
as one sparse matrix product so 1,400 gene sets cost one pass over the data.
"""
from __future__ import annotations

import os
import sys
import glob
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, "/workspace/code")
import sasp_real as R
import sasp_phase3 as P
import run_phase3_nulls as RN

RES = P.RESULTS
os.makedirs(RES, exist_ok=True)
N_RAND = 200
N_BINS = 25
CTRL_SIZE = 50
CORE_GS = "/usr/local/lib/python3.11/dist-packages/DeepScence/data/coreGS_v2.csv"
ORTHO = "/workspace/genesets/mouse_human_orthologs_MGI.csv"


def geneset(name):
    return [l.strip() for l in open(f"/workspace/genesets/{name}.txt") if l.strip()]


def corescence_mouse(panel, occ_min=5):
    """CoreScence v2 at occurrence >= 5, mapped human -> mouse the same way
    Bio did (the MGI HOM_MouseHumanSequence report the source authors used)."""
    cs = pd.read_csv(CORE_GS)
    cs = cs[cs["occurrence"] >= occ_min]
    o = pd.read_csv(ORTHO)
    cols = {c.lower(): c for c in o.columns}
    hm = [c for c in o.columns if "human" in c.lower()]
    mm = [c for c in o.columns if "mouse" in c.lower()]
    if hm and mm:
        m = dict(zip(o[hm[0]].astype(str).str.upper(), o[mm[0]].astype(str)))
    else:
        m = {}
    out, direction = [], {}
    for _, r in cs.iterrows():
        g = m.get(str(r["gene_symbol"]).upper())
        if g is None:
            g = str(r["gene_symbol"]).capitalize()
        if g in panel:
            out.append(g)
            direction[g] = r["direction"]
    return sorted(set(out)), direction


class Scorer:
    """score(S) = rowmean(L[:, S]) - rowmean(L[:, ctrl(S)]).

    Control genes are drawn from the same expression bins as the set, which is
    what `sc.tl.score_genes` does and what Bio's module scores use.
    """

    def __init__(self, L, gene_names, seed=P.MASTER_SEED):
        self.L = L.tocsc()
        self.g = np.asarray(gene_names)
        self.idx = {g: i for i, g in enumerate(self.g)}
        mu = np.asarray(L.mean(axis=0)).ravel()
        order = np.argsort(mu)
        self.bin = np.empty(mu.size, int)
        self.bin[order] = (np.arange(mu.size) * N_BINS) // mu.size
        self.rng = np.random.default_rng(seed)
        self.bin_members = [np.flatnonzero(self.bin == b) for b in range(N_BINS)]

    def cols(self, genes):
        return np.array([self.idx[g] for g in genes if g in self.idx], int)

    def control(self, cols):
        out = []
        for b, c in zip(*np.unique(self.bin[cols], return_counts=True)):
            pool = self.bin_members[b]
            out.append(self.rng.choice(pool, size=min(CTRL_SIZE * c, pool.size),
                                       replace=False))
        return np.concatenate(out)

    def weights(self, sets):
        """Sparse (n_genes, n_sets) weight matrix; score = L @ W."""
        rows, cs, vals = [], [], []
        for k, genes in enumerate(sets):
            c = self.cols(genes)
            if c.size == 0:
                continue
            ct = self.control(c)
            rows.append(c); cs.append(np.full(c.size, k)); vals.append(np.full(c.size, 1.0 / c.size))
            rows.append(ct); cs.append(np.full(ct.size, k)); vals.append(np.full(ct.size, -1.0 / ct.size))
        W = sp.csc_matrix((np.concatenate(vals),
                           (np.concatenate(rows), np.concatenate(cs))),
                          shape=(self.g.size, len(sets)))
        return W

    def score(self, sets, chunk=400):
        out = np.empty((self.L.shape[0], len(sets)), np.float32)
        for a in range(0, len(sets), chunk):
            W = self.weights(sets[a:a + chunk])
            out[:, a:a + W.shape[1]] = np.asarray((self.L @ W).todense(),
                                                  np.float32)
        return out


def _betas(d, Yr, lam, idx):
    """Standardised beta of exp(-d/lam) for each column of Yr, intercept only."""
    k = np.exp(-d[idx] / lam)
    kc = k - k.mean()
    kk = float(kc @ kc)
    Y = Yr[idx]
    num = kc @ Y
    sd = Y.std(axis=0)
    sd[sd < 1e-12] = np.nan
    return (num / kk) / sd


def run_section(sample, n_rand=N_RAND):
    sec = P.Sec(sample)
    M, gene_names, barcodes = R.load_expression(sample)
    keep = np.isin(barcodes, sec.cell_id)
    M = M[keep]
    bc = barcodes[keep]
    assert np.array_equal(bc, sec.cell_id), "cell order mismatch"
    L = R.normalize_counts(M)
    S = Scorer(L, gene_names)
    panel = set(gene_names)

    tierA = geneset("A_SENDER_FINAL_strict")
    bmods = {m: geneset(f"B_{m}") for m in P.MODULES}
    core, _ = corescence_mouse(panel)

    # ---- 1. disjointness -------------------------------------------------
    rows = []
    for m, gs in bmods.items():
        on = [g for g in gs if g in panel]
        rows.append(dict(section=sample, module=m, module_n=len(gs),
                         module_on_panel=len(on),
                         tierA_on_panel=len([g for g in tierA if g in panel]),
                         overlap_tierA=len(set(on) & set(tierA)),
                         corescence_on_panel=len(core),
                         overlap_corescence=len(set(on) & set(core)),
                         frac_corescence=len(set(on) & set(core)) / max(len(on), 1)))
    pd.DataFrame(rows).to_csv(f"{RES}/n8_disjointness_{sample}.csv", index=False)

    # ---- 2. scrambled response (Tier E3) ---------------------------------
    rand_sets, rand_key = [], []
    for j, m in enumerate(P.MODULES):
        tsv = pd.read_csv(f"/workspace/genesets/E3_random_matched/{m}.tsv",
                          sep="\t", header=None, names=["mod", "i", "genes"])
        for _, r in tsv.head(n_rand).iterrows():
            rand_sets.append(str(r["genes"]).split(","))
            rand_key.append((j, int(r["i"])))
    real_sets = [bmods[m] for m in P.MODULES]
    stripped = [[g for g in bmods[m] if g not in set(core)] for m in P.MODULES]
    Y = S.score(real_sets + stripped + rand_sets)
    Y_real = Y[:, :len(real_sets)]
    Y_strip = Y[:, len(real_sets):2 * len(real_sets)]
    Y_rand = Y[:, 2 * len(real_sets):]

    sf = RN.SectionFit(sample, RN.PRIMARY_CALL, P.MASTER_SEED)
    types = [t for t in sorted(set(sec.celltype)) if t not in P.EXCLUDE_TYPES]
    out = []
    for t in types:
        idx = sf.receivers(t)
        if idx.sum() < RN.MIN_RECEIVERS:
            continue
        ii = np.flatnonzero(idx)
        for j, m in enumerate(P.MODULES):
            y = Y_real[ii, j].astype(float)
            lam, beta, _, _ = P.profile_lambda(sf.d_obs[ii],
                                               np.ones((ii.size, 1)), y, sf.lam)
            bstd = _betas(sf.d_obs, Y_real[:, [j]].astype(float), lam, ii)[0]
            sel = [i for i, (jj, _) in enumerate(rand_key) if jj == j]
            br = _betas(sf.d_obs, Y_rand[:, sel].astype(float), lam, ii)
            br = br[np.isfinite(br)]
            out.append(dict(section=sample, arm=sec.meta["condition"],
                            celltype=t, module=m, n=int(ii.size), lam=lam,
                            beta_obs_std=float(bstd), n_rand=int(br.size),
                            rand_mean=float(br.mean()), rand_sd=float(br.std()),
                            rand_lo=float(np.quantile(br, .025)),
                            rand_hi=float(np.quantile(br, .975)),
                            rand_absmean=float(np.abs(br).mean()),
                            pct_rand_ge=float((br >= bstd).mean()),
                            sf_n8=float((bstd - br.mean()) / bstd)))
    pd.DataFrame(out).to_csv(f"{RES}/n8_scrambled_{sample}.csv", index=False)

    # ---- 3. DeepScence/CoreScence circularity ----------------------------
    corescore = S.score([core])[:, 0].astype(float)
    ok = ~np.isin(sec.celltype, P.EXCLUDE_TYPES + P.EXCLUDE_FROM_SENDERS)
    cs_mask = np.zeros(sec.n, bool)
    for c in np.unique(sec.celltype):
        s = (sec.celltype == c) & ok
        if s.sum() < 100:
            continue
        cs_mask[s] = corescore[s] > np.percentile(corescore[s], 95)
    d_core = P.dist_to_senders(sf.coords, cs_mask)

    crows = []
    for t in types:
        base = sf.receivers(t)
        if base.sum() < RN.MIN_RECEIVERS:
            continue
        for sender_name, dd, smask in (("corescence_p95", d_core, cs_mask),
                                       ("tierA_p95", sf.d_obs, sf.sender)):
            idx = base & (~smask) & (dd <= RN.WINDOW_UM)
            ii = np.flatnonzero(idx)
            if ii.size < RN.MIN_RECEIVERS:
                continue
            for j, m in enumerate(P.MODULES):
                y = Y_real[ii, j].astype(float)
                lam, _, _, _ = P.profile_lambda(dd[ii], np.ones((ii.size, 1)),
                                                y, sf.lam)
                b_full = _betas(dd, Y_real[:, [j]].astype(float), lam, ii)[0]
                b_strip = _betas(dd, Y_strip[:, [j]].astype(float), lam, ii)[0]
                crows.append(dict(section=sample, celltype=t, module=m,
                                  sender=sender_name, n=int(ii.size), lam=lam,
                                  beta_full_std=float(b_full),
                                  beta_stripped_std=float(b_strip),
                                  ratio_stripped=float(b_strip / b_full),
                                  n_stripped_genes=len(bmods[m]) - len(stripped[j])))
    pd.DataFrame(crows).to_csv(f"{RES}/n8_circularity_{sample}.csv", index=False)
    print(f"[n8] {sample} done", flush=True)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--sections", default="sbr")
    ap.add_argument("--n-rand", type=int, default=N_RAND)
    a = ap.parse_args()
    secs = {"sbr": P.SBR, "sham": P.SHAM, "all": P.ALL_SECTIONS,
            "inband": P.IN_BAND}.get(a.sections, a.sections.split(","))
    missing = [s for s in secs if not os.path.exists(
        os.path.join(P.CACHE3, f"{s}.npz"))]
    if missing:
        raise SystemExit(f"no Phase 3 cache for: {missing}")
    for s in secs:
        if os.path.exists(os.path.join(P.CACHE3, f"{s}.npz")):
            run_section(s, a.n_rand)
