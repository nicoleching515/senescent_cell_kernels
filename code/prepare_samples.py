"""One-time per-sample preparation: parse the h5, score modules, call senders,
build cKDTree geometry, and cache everything to a compact .npz.

Parsing cell_feature_matrix.h5 and building the neighbour statistics costs
~30-75 s per section; the kernel fits themselves cost ~1 s.  Doing the parse
once and caching means the fitting stage (which is run many times, across
modules x sender calls x families x bootstraps) starts instantly, and re-running
after the Bio agent delivers its annotations is a single cheap re-prepare.

Cache is ~20 MB per section (float32), well inside the workspace quota.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import sasp_real as R

CACHE = "/workspace/data/processed/cache"


def prepare(sample: str, force: bool = False) -> str:
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, f"{sample}.npz")
    if os.path.exists(out) and not force:
        return f"[skip] {sample}"
    s = R.load_sample(sample)
    d = dict(
        coords=s.coords.astype(np.float32),
        celltype=np.asarray(s.celltype).astype(str),
        median_nn_um=np.float64(s.median_nn_um),
        cell_id=s.cells["cell_id"].to_numpy().astype(str),
        transcript_counts=s.cells["transcript_counts"].to_numpy(np.float32),
        cell_area=s.cells["cell_area"].to_numpy(np.float32),
        nn1_um=s.cells["nn1_um"].to_numpy(np.float32),
    )
    for c in ("density_25um", "density_50um", "density_100um"):
        d[c] = s.cells[c].to_numpy(np.float32)
    for c in s.cells.columns:
        if c.startswith("anat_") or c == "tierA_score":
            d[c] = s.cells[c].to_numpy(np.float32)
    for k, v in s.modules.items():
        d[f"mod__{k}"] = np.asarray(v, np.float32)
    for k, v in s.sender.items():
        d[f"snd__{k}"] = np.asarray(v, bool)
        # distance to nearest sender for THIS call, precomputed once
        d[f"dist__{k}"] = R.distance_to_set(s.coords, np.asarray(v, bool)
                                            ).astype(np.float32)
    d["meta"] = np.array([f"{k}={v}" for k, v in s.meta.items()])
    d["sources"] = np.array([f"celltype={s.celltype_source}",
                             f"sender={s.sender_source}",
                             f"module={s.module_source}"])
    np.savez_compressed(out, **d)
    return (f"[done] {sample} n={s.n()} medNN={s.median_nn_um:.2f} "
            f"modules={len(s.modules)} senders={len(s.sender)}")


class Cached:
    """Lightweight accessor over the cached .npz."""

    def __init__(self, sample: str):
        self.name = sample
        self.z = np.load(os.path.join(CACHE, f"{sample}.npz"),
                         allow_pickle=False)
        self.meta = dict(x.split("=", 1) for x in self.z["meta"])
        self.sources = dict(x.split("=", 1) for x in self.z["sources"])

    @property
    def coords(self):
        return self.z["coords"].astype(float)

    @property
    def celltype(self):
        return self.z["celltype"]

    @property
    def median_nn_um(self):
        return float(self.z["median_nn_um"])

    def modules(self):
        return [k[5:] for k in self.z.files if k.startswith("mod__")]

    def senders(self):
        return [k[5:] for k in self.z.files if k.startswith("snd__")]

    def module(self, m):
        return self.z[f"mod__{m}"].astype(float)

    def sender(self, k):
        return self.z[f"snd__{k}"]

    def dist(self, k):
        return self.z[f"dist__{k}"].astype(float)

    def col(self, c):
        return self.z[c]


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--n-jobs", type=int, default=4)
    a = ap.parse_args()
    res = Parallel(n_jobs=a.n_jobs, prefer="processes")(
        delayed(prepare)(s, a.force) for s in R.list_samples())
    for r in res:
        print(r, flush=True)
