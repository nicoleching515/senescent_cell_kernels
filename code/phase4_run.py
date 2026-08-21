"""Phase 4 runner — every method, on real and on coordinate-shuffled tissue.

  python3 phase4_run.py --method cellchat|spatalk|ncem|commot [--reps 100] [--jobs 16]

Writes one CSV per (method, tile) to /workspace/results/phase4/parts/ so a
crash costs one tile, not the run (Master Plan 18.4).
"""
from __future__ import annotations
import argparse, os, time, warnings
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd
from scipy.spatial import cKDTree

import phase4_data as D
import phase4_methods as M
from phase4_tiles import all_tiles

warnings.filterwarnings("ignore")
OUT = "/workspace/results/phase4/parts"
MASTER_SEED = 20260820
NULLS = ["N3_lig", "N4_lig", "N3_type", "N0_perm"]
# Added after the D7 audit (sec 2.1 B9).  Kept OUT of NULLS on purpose: the job
# seeds below are a function of the enumeration order, so appending to NULLS
# would silently change every seed in the original run.  Reached via --only-null.
COORD_NULLS = ("N3_type", "N0_perm", "N0_type")
ALPHA = 0.05
MIN_LIG = 20          # ligand+ cells needed before a ligand-field null is meaningful
NCEM_RADII = [10., 15., 20., 30., 40., 50., 75., 100.]
CELLCHAT_PERM = 100      # CellChat default nboot
SPATALK_PERM = 200       # reduced from the SpaTalk default of 1000 for runtime;
                         # p-resolution 0.005, ample for alpha = 0.05


def lig_masks(tile):
    return {p["pair"]: (tile.gc(p["ligand"]) > 0) for p in D.LR_PAIRS}


def pair_vectors(tile):
    return [(tile.g(p["ligand"]), tile.receptor(p["receptors"])) for p in D.LR_PAIRS]


def _emit(rows, tile, method, cond, rep, pnames, scores, pvals, extra=None):
    K = tile.K
    for t, pn in enumerate(pnames):
        S, P = scores[t], pvals[t]
        for i in range(K):
            for j in range(K):
                r = dict(tile=tile.name, section=tile.section, arm=tile.meta["condition"],
                         week=tile.meta["week"], method=method, pair=pn,
                         sender=tile.types[i], receiver=tile.types[j],
                         cond=cond, rep=rep, score=float(S[i, j]), p=float(P[i, j]))
                if extra:
                    r.update(extra)
                rows.append(r)


# ---------------------------------------------------------------------------

COMMOT_DIS_THR = D.WINDOW_UM      # 100 um, the project-wide window (CS_PHASE3 s2)


def _commot_one(tile, xy, pdef, seed):
    """One COMMOT run: published software, one pathway, one coordinate set."""
    import _shims.np2_compat            # noqa: F401  (must precede commot)
    import anndata as ad, commot as ct
    A = ad.AnnData(np.asarray(tile.X, dtype=np.float64),
                   obs=pd.DataFrame({"celltype": tile.celltype},
                                    index=[str(i) for i in range(tile.n)]),
                   var=pd.DataFrame(index=tile.genes))
    A.obsm["spatial"] = np.asarray(xy, dtype=np.float64)
    df = pd.DataFrame([[pdef["ligand"], r, pdef["pathway"]] for r in pdef["receptors"]])
    ct.tl.spatial_communication(A, database_name="sasp", df_ligrec=df,
                                dis_thr=COMMOT_DIS_THR, heteromeric=False,
                                pathway_sum=True)
    ct.tl.cluster_communication(A, database_name="sasp",
                                pathway_name=pdef["pathway"],
                                clustering="celltype", n_permutations=100,
                                random_seed=int(seed) % (2 ** 31 - 1))
    r = A.uns["commot_cluster-celltype-sasp-" + pdef["pathway"]]
    names = list(r["communication_matrix"].index)
    S = r["communication_matrix"].reindex(index=tile.types, columns=tile.types).values
    P = r["communication_pvalue"].reindex(index=tile.types, columns=tile.types).values
    assert len(names) == tile.K, (len(names), tile.K)
    return S, P


def _ncem_sweep(tile, rng, lm, pnames):
    """NCEM's own model selection: neighbourhood radius maximising variance
    explained.  Run on real coordinates and again under each null, because a
    method whose LENGTH SCALE is unchanged by coordinate shuffling has not
    measured a length scale."""
    K = tile.K
    sweep = []
    confs = [("real", tile.coords)]
    for null in NULLS:
        xy = (D.null_coords(tile, null, rng) if null in COORD_NULLS
              else D.null_coords(tile, null, rng, lig_mask=lm[pnames[0]]))
        if xy is not None:
            confs.append((null, xy))
    for label, xy in confs:
        for r in NCEM_RADII:
            r2s = [M.ncem_linear(tile.receptor(p["receptors"]), xy, tile.code,
                                 K, r)[2] for p in D.LR_PAIRS]
            sweep.append(dict(tile=tile.name, section=tile.section,
                              cond=label, radius=r, r2=float(np.nanmean(r2s))))
    sw = pd.DataFrame(sweep)
    sw.to_csv(f"{OUT}/ncem_sweep_{tile.name}.csv", index=False)
    real = sw[sw.cond == "real"]
    return float(real.loc[real.r2.idxmax(), "radius"])


def run_job(tile, method, cond, reps, seed, pair_idx=None):
    """One (tile, method, condition[, LR pair]) unit of work."""
    rng = np.random.default_rng(seed)
    pnames = [p["pair"] for p in D.LR_PAIRS]
    lm = lig_masks(tile)
    rows, K, d0 = [], tile.K, tile.med_nn
    sel_all = list(range(len(D.LR_PAIRS))) if pair_idx is None else [pair_idx]
    ncem_radius = _ncem_sweep(tile, np.random.default_rng(seed + 7), lm, pnames) \
        if method == "ncem" else None

    def stat_all(xy, sel):
        pv = pair_vectors(tile)
        if method == "cellchat":
            s, p = M.cellchat_many([pv[i] for i in sel], xy, tile.code, K, d0,
                                   rng, n_perm=CELLCHAT_PERM)
            return s, p, {}
        if method == "spatalk":
            s, p, _ = M.spatalk_many([pv[i] for i in sel], xy, tile.code, K,
                                     rng, n_perm=SPATALK_PERM)
            return s, p, {}
        if method == "ncem":
            S, Q = [], []
            for i in sel:
                B, P, _ = M.ncem_linear(tile.receptor(D.LR_PAIRS[i]["receptors"]),
                                        xy, tile.code, K, ncem_radius)
                S.append(B); Q.append(M.bh(P))
            return S, Q, dict(radius=ncem_radius)
        if method == "commot":
            S, P = [], []
            for i in sel:
                a, b = _commot_one(tile, xy, D.LR_PAIRS[i], rng.integers(1, 2 ** 30))
                S.append(a); P.append(b)
            return S, P, {}
        raise ValueError(method)

    if cond == "real":
        s, p, ex = stat_all(tile.coords, sel_all)
        _emit(rows, tile, method, "real", -1, [pnames[i] for i in sel_all], s, p, ex)
        return pd.DataFrame(rows)

    for rep in range(reps):
        if cond in COORD_NULLS:
            xy = D.null_coords(tile, cond, rng)
            s, p, ex = stat_all(xy, sel_all)
            _emit(rows, tile, method, cond, rep, [pnames[i] for i in sel_all],
                  s, p, ex)
        else:
            for i in sel_all:
                m = lm[pnames[i]]
                if m.sum() < MIN_LIG:
                    continue
                xy = D.null_coords(tile, cond, rng, lig_mask=m)
                s, p, ex = stat_all(xy, [i])
                _emit(rows, tile, method, cond, rep, [pnames[i]], s, p, ex)
    return pd.DataFrame(rows)


def _job(args):
    tile, method, cond, reps, seed, pair_idx = args
    tag = f"{method}_{tile.name}_{cond}" + ("" if pair_idx is None else f"_p{pair_idx}")
    f = f"{OUT}/{tag}.csv"
    if os.path.exists(f):
        return f"[skip] {tag}"
    t0 = time.time()
    df = run_job(tile, method, cond, reps, seed, pair_idx)
    df.to_csv(f, index=False)
    return f"[ok] {tag} n={tile.n} rows={len(df)} {time.time()-t0:.0f}s"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", required=True)
    ap.add_argument("--reps", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--only-null", default=None,
                    help="run a single condition (e.g. N0_type) without "
                         "disturbing the original job enumeration or seeds")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    tiles = all_tiles()
    split_pairs = (a.method == "commot")
    jobs = []
    for i, t in enumerate(tiles):
        for c in ([a.only_null] if a.only_null else ["real"] + NULLS):
            if split_pairs:
                for pi in range(len(D.LR_PAIRS)):
                    jobs.append((t, a.method, c, a.reps,
                                 MASTER_SEED + 1000 * i + 37 * len(jobs), pi))
            else:
                jobs.append((t, a.method, c, a.reps,
                             MASTER_SEED + 1000 * i + 37 * len(jobs), None))
    # Round-robin over SECTIONS by tile index, then largest tile first.  If the
    # run has to be stopped early, every section is covered by its first tile
    # rather than one section being covered three times.
    jobs.sort(key=lambda j: (int(j[0].name.split("_t")[1]), -j[0].n,
                             j[2] in COORD_NULLS or j[2] == "real"))
    print(f"{len(jobs)} jobs, {a.jobs} workers", flush=True)
    if a.jobs == 1:
        for j in jobs:
            print(_job(j), flush=True)
    else:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.jobs, maxtasksperchild=4) as pool:
            for r in pool.imap_unordered(_job, jobs):
                print(r, flush=True)


if __name__ == "__main__":
    main()
