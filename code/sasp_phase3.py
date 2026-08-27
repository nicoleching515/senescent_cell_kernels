"""
Phase 3 — full null battery N1-N8 on the real SBR arm.

Master Plan Sections 22 Step 3, 23, 24, 25; Bio Phase 2 constraints.

Design decisions carried in from Phases 1-2 (do NOT re-derive):
  * report SURVIVING FRACTIONS, never permutation p-values (Sec 6.5; CS Phase 2 §10);
  * hold lambda FIXED at lambda_hat_naive inside every null, otherwise the null
    compares different models and beta stops being an amplitude (CS Phase 2 §10);
  * N2 uses the SHARED-lambda two-kernel form, not beta_true - beta_decoy;
  * the fitting window is set from the observed distance distribution, not from
    the plan's unreachable 300 um (CS Phase 2 §2);
  * every estimate carries a spatial block bootstrap, never an iid CI.

This module is data prep + the estimator core.  Runners are in
`run_phase3_nulls.py`.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import sasp_real as R

PROC = "/workspace/data/processed"
CACHE3 = os.path.join(PROC, "cache3")
RESULTS = "/workspace/results/phase3"
MASTER_SEED = 20260820

# Bio Phase 2 §6: not receiver cell types, and not sender candidates.
EXCLUDE_TYPES = ("Low_quality", "Unknown", "unknown")
# Bio Phase 2 §4.4: Cdkn1a is induced in cycling cells; not senescence-specific.
EXCLUDE_FROM_SENDERS = ("Proliferating",)

# Receiver stratification label set: "merged" (Bio Phase 3 stable label
# families) or "fine" (Bio Phase 2 labels).  Merged is the default because the
# fine labels are not comparable across sections.
LABELS = "merged"

MODULES = ["downstream_arrest", "emt_ecm", "il6_jak_stat3", "interferon_response",
           "oxidative_stress", "secondary_senescence", "tnfa_nfkb_proximal"]

# Union of Bio's receiver cell types over the six annotated sections.  Used to
# harmonise the k-NN-composition covariate block when sections are POOLED, so
# that the design matrix has the same columns in every section (Phase 3 T2).
CANON_TYPES = ("B-cells", "Biliary/ductular", "Central venous LSECs", "DC",
               "Hepatic stellate cells", "Hepatocytes", "Inflammatory macs",
               "Kupffer cells", "LSECs", "Mesothelial cells",
               "Portal endothelial cells", "Proliferating", "T/NK cells",
               "vSMCs")

# ---------------------------------------------------------------------------
# Section inventory and the Section 8 Test 3 admissibility rule.
#
# Test 3 requires sender prevalence in 1-20 %.  Above the ceiling,
# distance-to-nearest-sender is near zero everywhere and lambda is
# unidentifiable BY CONSTRUCTION; below the floor there are too few senders.
# Measured on Cdkn1a+ hepatocytes across all 11 sections (Bio Phase 3):
#   over ceiling : 7239 45.0 %, 7448 25.6 %, 7361 25.3 %, 7450 22.7 %
#   in band      : 7260 10.5, 7259 9.6, 7001 8.9, 7352 7.2, 7248 4.9, 7435 2.3
#   below floor  : 7250 0.48 %
# The band cuts ACROSS the surgical arms, so the primary analysis is the six
# in-band sections (six animals, both arms) with arm as a contrast, NOT the
# SBR arm alone.  See CS_PHASE3.md.
IN_BAND = ["7259_liver_sbr_Male_26-U1", "7260_liver_sbr_Male_26-U1",
           "7001_liver_sham_Male_52-U1", "7248_liver_sham_Male_26-U1",
           "7352_liver_sham_Male_2-U1", "7435_liver_sham_Male_10-U1"]
OVER_CEILING = ["7239_liver_sbr_Male_52-U1", "7448_liver_sbr_Male_10-U1",
                "7361_liver_sbr_Male_2-U1", "7450_liver_sbr_Male_10-U1"]
BELOW_FLOOR = ["7250_liver_sham_Male_26-U1"]
ALL_SECTIONS = IN_BAND + OVER_CEILING + BELOW_FLOOR

# kept for the Phase-3 sensitivity analysis that reproduces the first,
# SBR-only scoping
# merged-label version of CANON_TYPES (Bio Phase 3 stable label families)
CANON_TYPES_MERGED = ("B-cells", "Biliary/ductular", "DC", "Endothelial",
                      "Hepatocytes", "Low_quality", "Macrophages",
                      "Mesenchymal", "Mesothelial cells", "Proliferating",
                      "T/NK cells", "vSMCs")

SBR = ["7259_liver_sbr_Male_26-U1", "7361_liver_sbr_Male_2-U1",
       "7450_liver_sbr_Male_10-U1", "7239_liver_sbr_Male_52-U1"]
SHAM = ["7250_liver_sham_Male_26-U1", "7435_liver_sham_Male_10-U1"]


# ---------------------------------------------------------------------------
# preparation
# ---------------------------------------------------------------------------

def prep(sample: str, force: bool = False) -> str:
    """Cache everything the null battery needs for one section.

    Analysis cell set = cells that (a) pass the Xenium QC used upstream and
    (b) carry a Bio cell-type label.  Anything else is dropped here rather
    than silently defaulted, so the exclusion is auditable.
    """
    os.makedirs(CACHE3, exist_ok=True)
    out = os.path.join(CACHE3, f"{sample}.npz")
    if os.path.exists(out) and not force:
        return f"[skip] {sample}"

    cells = pd.read_parquet(os.path.join(R.RAW, sample, "cells.parquet"))
    M, gene_names, barcodes = R.load_expression(sample)
    if not np.array_equal(barcodes, cells["cell_id"].to_numpy().astype(str)):
        cells = cells.set_index("cell_id").loc[barcodes].reset_index()

    def _read(kind):
        p = os.path.join(PROC, f"{kind}_{sample}.csv")
        if not os.path.exists(p):
            raise FileNotFoundError(p)
        return pd.read_csv(p).set_index("cell_id")

    ct = _read("celltypes")
    mods = _read("modules")
    sen = _read("senders")
    anat = _read("anatomy")

    cid = cells["cell_id"].to_numpy().astype(str)
    keep = pd.Index(cid).isin(ct.index).astype(bool)
    cells = cells.loc[keep].reset_index(drop=True)
    M = M[keep]
    cid = cid[keep]

    ct = ct.reindex(cid)
    mods = mods.reindex(cid)
    sen = sen.reindex(cid)
    anat = anat.reindex(cid)

    coords = np.column_stack([cells["x_centroid"].to_numpy(float),
                              cells["y_centroid"].to_numpy(float)])
    tree = cKDTree(coords)
    d = {}
    for r in (25.0, 50.0, 100.0):
        d[f"density_{int(r)}um"] = np.asarray(
            tree.query_ball_point(coords, r=r, return_length=True,
                                  workers=-1), np.float32) - 1.0
    nn, nn_idx = tree.query(coords, k=21, workers=-1)
    d["nn1_um"] = nn[:, 1].astype(np.float32)
    d["knn_idx"] = nn_idx[:, 1:].astype(np.int32)
    d["median_nn_um"] = np.float64(np.median(nn[:, 1]))

    genes_det = np.asarray((M > 0).sum(axis=1)).ravel().astype(np.float32)

    segm = cells["segmentation_method"].to_numpy().astype(str)
    seg_levels = np.array(sorted(set(segm)))
    d["seg_levels"] = seg_levels
    d["seg_code"] = np.searchsorted(seg_levels, segm).astype(np.int8)

    d.update(
        cell_id=cid,
        coords=coords.astype(np.float32),
        celltype=ct["cell_type"].to_numpy().astype(str),
        # Bio Phase 3 §1.1: three label families are not stably separable on
        # this panel and the winner flips between sections (7001 lost its
        # entire stellate compartment to `Unknown` on a 0.02 margin).  The
        # merged label set is what cross-section fits must use.
        celltype_merged=(ct["cell_type_merged"].to_numpy().astype(str)
                         if "cell_type_merged" in ct.columns
                         else ct["cell_type"].to_numpy().astype(str)),
        median_depth=np.float64(np.median(
            cells["transcript_counts"].to_numpy(float))),
        ct_conf=ct["cell_type_confidence"].to_numpy(np.float32),
        transcript_counts=cells["transcript_counts"].to_numpy(np.float32),
        genes_detected=genes_det,
        cell_area=cells["cell_area"].to_numpy(np.float32),
        nucleus_area=cells["nucleus_area"].to_numpy(np.float32),
        zonation_score=anat["zonation_score"].to_numpy(np.float32),
        compartment=anat["compartment_label"].to_numpy().astype(str),
        dist_to_boundary_um=anat["dist_to_boundary_um"].to_numpy(np.float32),
        dist_to_portal_triad_um=anat["dist_to_portal_triad_um"].to_numpy(np.float32),
        cdkn1a_counts=sen["cdkn1a_counts"].to_numpy(np.float32),
        tierA_score=sen["tierA_score"].to_numpy(np.float32),
        senepy_score=sen["senepy_score"].to_numpy(np.float32),
    )
    for q in (90, 95, 99):
        d[f"flag_p{q}"] = sen[f"sender_flag_p{q}"].fillna(0).to_numpy().astype(bool)
    for m in MODULES:
        d[f"mod__{m}"] = mods[m].to_numpy(np.float32)
    # Phase 8 / D1: the seven per-module Tier A sensitivity sender sets
    # (`genesets/A_sender_for_<module>.txt`), scored and thresholded by
    # `phase2_downstream.py` exactly as the PRIMARY Tier A set is.  Cached
    # here so `Sec.sender_mask("tierApm_pNN", module=...)` costs nothing.
    # A cache written before Phase 8 simply lacks these keys, and
    # `sender_mask` raises rather than guessing (see `Sec.has_permodule`).
    for m in MODULES:
        col = f"tierA_{m}_score"
        if col not in sen.columns:
            continue
        d[f"tierApm__{m}"] = sen[col].to_numpy(np.float32)
        for q in (90, 95, 99):
            d[f"flag_pm_{m}_p{q}"] = (sen[f"sender_flag_{m}_p{q}"]
                                      .fillna(0).to_numpy().astype(bool))

    np.savez_compressed(out, **d)
    return (f"[done] {sample} n={len(cid)} medNN={float(d['median_nn_um']):.2f} "
            f"seg_levels={len(seg_levels)}")


class Sec:
    """Accessor over the Phase 3 cache for one section."""

    def __init__(self, sample: str, labels: str = None):
        self.name = sample
        self.z = np.load(os.path.join(CACHE3, f"{sample}.npz"), allow_pickle=False)
        self.meta = R.parse_sample(sample)
        self._c = {}
        self.labels = labels or LABELS
        if self.labels == "merged" and "celltype_merged" in self.z.files:
            self._c["celltype"] = self.z["celltype_merged"]
        self.meta["band"] = ("in_band" if sample in IN_BAND else
                             "over_ceiling" if sample in OVER_CEILING
                             else "below_floor")

    def __getattr__(self, k):
        if k in ("z", "_c"):
            raise AttributeError(k)
        if k not in self._c:
            self._c[k] = self.z[k]
        return self._c[k]

    @property
    def n(self):
        return self.coords.shape[0]

    def module(self, m):
        return self.z[f"mod__{m}"].astype(float)

    @property
    def has_permodule(self) -> bool:
        """True when this cache carries the Phase 8 per-module Tier A flags."""
        return all(f"flag_pm_{m}_p95" in self.z.files for m in MODULES)

    # ---- sender calls (N7 axis) ------------------------------------------
    def sender_mask(self, call: str, module: str = None) -> np.ndarray:
        """Sender calls, all excluding Low_quality/Unknown/Proliferating.

        tierA_pNN   : Bio's within-cell-type NNth-percentile Tier A flag
        cdkn1a_pos  : Cdkn1a > 0, the call the source paper uses
        senepy_pNN  : within-cell-type NNth percentile of the SenePy hub score
        tierApm_pNN : the same within-cell-type NNth-percentile rule applied to
                      the PER-MODULE Tier A sensitivity set
                      `genesets/A_sender_for_<module>.txt` (Phase 8 D1).  This
                      call is module-specific, so `module=` is required and the
                      caller must fan out over `MODULES`.
        """
        ok = ~np.isin(self.celltype, EXCLUDE_TYPES + EXCLUDE_FROM_SENDERS)
        if call.startswith("tierApm_p"):
            # checked BEFORE tierA_p*, which is not a prefix of it but is easy
            # to confuse; the two families must never silently swap.
            if module is None:
                raise ValueError(
                    f"{call!r} is a per-module sender call: pass "
                    "module=<one of sasp_phase3.MODULES>.")
            if module not in MODULES:
                raise ValueError(f"unknown module {module!r}")
            key = f"flag_pm_{module}_p{call[9:]}"
            if key not in self.z.files:
                raise ValueError(
                    f"{self.name}: the Phase 3 cache has no {key!r}.  It "
                    "predates Phase 8; re-run phase2_downstream.py and "
                    "sasp_phase3.prep(force=True) before using "
                    f"{call!r}.  (Sec.has_permodule reports this.)")
            m = self.z[key].copy()
        elif call.startswith("tierA_p"):
            m = self.z[f"flag_p{call[7:]}"].copy()
        elif call == "cdkn1a_pos":
            m = self.cdkn1a_counts > 0
        elif call.startswith("senepy_p"):
            q = float(call[8:])
            s = self.senepy_score.astype(float)
            m = np.zeros(self.n, bool)
            for c in np.unique(self.celltype):
                sel = (self.celltype == c) & np.isfinite(s) & ok
                if sel.sum() < 100:
                    continue
                m[sel] = s[sel] > np.percentile(s[sel], q)
        else:
            raise ValueError(call)
        return m & ok


# ---------------------------------------------------------------------------
# geometry helpers
# ---------------------------------------------------------------------------

def dist_to_senders(coords: np.ndarray, mask: np.ndarray) -> np.ndarray:
    if mask.sum() == 0:
        return np.full(coords.shape[0], np.nan)
    dd, _ = cKDTree(coords[mask]).query(coords, k=1, workers=-1)
    return dd


def dist_to_points(coords: np.ndarray, pts: np.ndarray) -> np.ndarray:
    dd, _ = cKDTree(pts).query(coords, k=1, workers=-1)
    return dd


def block_ids(coords: np.ndarray, n_side: int = 10) -> np.ndarray:
    """Quantile-defined spatial blocks (equal-count, so no empty blocks)."""
    qx = np.quantile(coords[:, 0], np.linspace(0, 1, n_side + 1)[1:-1])
    qy = np.quantile(coords[:, 1], np.linspace(0, 1, n_side + 1)[1:-1])
    bx = np.searchsorted(qx, coords[:, 0])
    by = np.searchsorted(qy, coords[:, 1])
    return (bx * n_side + by).astype(np.int64)


# ---------------------------------------------------------------------------
# the estimator core
# ---------------------------------------------------------------------------

def _qr(X):
    Q, _ = np.linalg.qr(X)
    return Q


class FixedLambdaFitter:
    """beta-hat for y = X gamma + beta * exp(-d/lam) + eps, at FIXED lam.

    X is orthonormalised once (thin QR); each null then costs one (n x p)
    matvec pair plus dot products, so 1,000 permutations over 100k cells with
    30 covariates is seconds, not minutes (Sec 18.2).
    """

    def __init__(self, X: np.ndarray, Y: np.ndarray):
        self.Q = _qr(np.asarray(X, float))
        Y = np.asarray(Y, float)
        self.Y = Y - self.Q @ (self.Q.T @ Y)      # residualised responses
        self.p = X.shape[1]

    def beta(self, k: np.ndarray) -> np.ndarray:
        kt = k - self.Q @ (self.Q.T @ k)
        kk = float(kt @ kt)
        if kk <= 1e-12:
            return np.full(self.Y.shape[1], np.nan)
        return (kt @ self.Y) / kk

    def beta_rss(self, k):
        kt = k - self.Q @ (self.Q.T @ k)
        kk = float(kt @ kt)
        b = (kt @ self.Y) / kk
        rss = (self.Y * self.Y).sum(0) - b ** 2 * kk
        return b, rss


def profile_lambda(d, X, y, lam_grid):
    """Profile the exponential kernel's lambda on a grid; exact linear solve
    at each grid point.  Returns (lam_hat, beta_hat, rss curve)."""
    f = FixedLambdaFitter(X, y[:, None])
    rss = np.empty(lam_grid.size)
    bet = np.empty(lam_grid.size)
    for i, lam in enumerate(lam_grid):
        b, r = f.beta_rss(np.exp(-d / lam))
        bet[i] = b[0]
        rss[i] = r[0]
    t = int(np.argmin(rss))
    return float(lam_grid[t]), float(bet[t]), rss, t


def standardize(Z):
    Z = np.asarray(Z, float)
    mu = np.nanmean(Z, axis=0)
    sd = np.nanstd(Z, axis=0)
    sd[sd < 1e-12] = 1.0
    Z = (Z - mu) / sd
    return np.nan_to_num(Z)


def dummies(labels, drop_first=True):
    lev = np.array(sorted(set(labels)))
    if drop_first:
        lev = lev[1:]
    if lev.size == 0:
        return np.zeros((len(labels), 0))
    return np.column_stack([(labels == l).astype(float) for l in lev]), lev


def smd(Z, a_idx, b_idx):
    """Standardized mean difference per covariate (Sec 8 Test 5)."""
    A, B = Z[a_idx], Z[b_idx]
    ma, mb = A.mean(0), B.mean(0)
    pooled = np.sqrt(0.5 * (A.var(0) + B.var(0))) + 1e-12
    return (ma - mb) / pooled
