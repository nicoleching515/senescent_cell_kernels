#!/usr/bin/env python3
"""
Phase 8 / PI decision D15 -- THE COMPOSITION-MATCHED RERUN PROTOCOL AT 5 SEEDS.

`Phase7_Minimal_Human_Replication (1).md` s15 freezes "the composition-matched
rerun protocol at 5 seeds"; s16 run-order step 9 says "Composition-matched
reruns, 5 seeds, both arms"; s17 has a row "Composition surrogate share".
`reports/PREREG_PHASE8.md` s3.8 records that NO code implemented it and asks the
PI to pre-register "its matching variables and its five literal seeds".  This
file is that implementation.  Everything the planning documents left unstated is
listed under SPECIFICATION below with the choice made and the reason.

--------------------------------------------------------------------------
SPECIFICATION
--------------------------------------------------------------------------
FOUND IN THE DOCUMENTS

  * The confound.  Master Plan s6.1 puts k-NN cell-type composition in the
    nuisance vector z_i; s22 Step 3 N2 defines the matched-decoy control as
    "a non-senescent cell matched on cell type, local density (k-NN within
    50 um), and k-NN composition"; s8 Test 5 sets the matching acceptance
    criterion (|SMD| <= 0.1).
  * The size of the confound.  CS_PHASE3 s3.1: receiver cell-type composition
    is 66 % of the naive amplitude (SF 0.344 with receiver-type intercepts
    alone).  CS_PHASE3 s3.2: the composition sub-block of N5 alone leaves
    SF 0.474.  CS_PHASE5 s4: the composition-only surrogate curve reproduces
    76 % of the unstratified contact amplitude, per-module ratios 0.42-1.83.
    Those two numbers are the "66-76 %" of the s17 row.
  * The estimand.  s6.5 / CS_PHASE3: report a SURVIVING FRACTION, lambda held
    fixed at lambda_hat_naive, with a 400-replicate spatial block bootstrap
    over 100 quantile blocks (s15 freezes those two).

NOT FOUND ANYWHERE -- resolved here, and these are the choices to pre-register

  D15.1 Which cells are matched.  The documents say "composition-matched" but
        never name the design.  Resolved: the N2 design of Master Plan s22
        Step 3 -- senders vs non-sender decoys -- with the matching set
        REDUCED TO COMPOSITION.  Reason: N2 is the only sender/non-sender
        matching design in the plan, it is already implemented
        (`phase3_core.match_decoys_section`), and reducing its matching set to
        composition is exactly the "is it composition rather than distance"
        question.  Variant `comp` below.
  D15.2 The matching variables.  Resolved: EXACT stratification on receiver
        cell type (the greedy matcher's stratum) plus 1-1 nearest-neighbour
        propensity matching, without replacement, caliper 0.25 SD, on the
        20-NN cell-type composition vector (`knn_frac_*`, one column per
        cell type present in the section).  Nothing else -- no density, no
        transcript depth, no anatomy.  Reason: the point of the rerun is to
        isolate composition; leaving density and depth in makes it N2 again.
        Variant `full` (the published N2 matching set: log_dens50,
        log_counts, anatomy, plus composition) is run alongside at the same
        five seeds so the isolation is measurable rather than asserted.
  D15.3 The five literal seeds.  Resolved:
        COMPMATCH_SEEDS = (20260901, 20260902, 20260903, 20260904, 20260905).
        Reason: derived from `sasp_phase3.MASTER_SEED = 20260820` by date, and
        chosen outside the range reachable by `run_phase3_nulls._expand`
        (MASTER_SEED + 1000*i + j, and + PM_SEED_OFFSET) so no composition
        rerun reproduces an existing N2 match by accident.
  D15.4 What the seed controls.  Resolved: the seed is passed to
        `match_decoys_section`, i.e. it fixes the order in which senders claim
        decoys in the greedy matcher, AND it seeds the block bootstrap.  Point
        estimates depend only on the matching draw (they are computed before
        the bootstrap), so the between-seed spread reported here is matching
        variability and nothing else.
  D15.5 Receiver scope.  Resolved: BOTH.  `scope == "ALL"` is the pooled,
        unstratified fit -- the one the field reports and the one CS_PHASE3
        s3.1 and CS_PHASE5 s4 measured composition against, so it is the row
        that fills the s17 "Composition surrogate share" cell.  Per-receiver-
        cell-type rows are reported as the secondary, since a within-type fit
        has already removed the between-type composition effect.
  D15.6 Sender calls.  Resolved: `tierA_p95` (A_SENDER_FINAL_strict, 33 genes,
        PRIMARY) and `tierApm_p95` (the seven per-module A_sender_for_*.txt
        sensitivity sets), which are the two pre-registered Tier A variants.
  D15.7 Section set.  Resolved: the six Section 8 Test 3 admissible ("in
        band") sections, `sasp_phase3.IN_BAND`, matching every other Phase 3
        primary.
  D15.9 Matching alone cannot answer the scientific question.  CS_PHASE3 s5
        records that in this dataset N2 matching is nearly inert (SF 0.943)
        while REGRESSING on the same covariates removes 92 %: "matching
        balances the covariates between senders and decoys; it does not remove
        the dependence of the response on those covariates at the receiver,
        which is where the confounding acts."  Resolved: a third, seed-free
        set of variants is computed at the same scopes, same lambda_hat, same
        bootstrap: `comp_adj` (20-NN composition as covariates -- the `comp`
        sub-block of `run_phase3_attribution.py`, extended to the pooled
        scope), `type_adj` (the receiver's own cell-type intercepts, the
        mu_{c_i} of s6.1 -- this is the block CS_PHASE3 s3.1 measured at
        SF 0.344, i.e. the 66 % end of the s17 range) and `typecomp_adj`
        (both).  `type_adj`/`typecomp_adj` are pooled-scope only, being
        degenerate inside a single-cell-type fit.  Without these the protocol
        would report "composition matching removes nothing" and a reader would
        wrongly hear "composition is not the confound".
  D15.8 Reportable population for the summary.  Resolved: identical to
        `summarize_phase3.sf_table` -- fits with beta_naive > 0 AND
        beta_base_lo > 0 -- so the medians here are directly comparable to the
        published SF table.

--------------------------------------------------------------------------
REUSE
--------------------------------------------------------------------------
No estimator, matcher or bootstrap is reimplemented here.
  * `phase3_core.match_decoys_section` / `greedy_ps_match` -- the N2 matcher.
    Only its `Zmatch` argument changes.
  * `run_phase3_nulls.SectionFit` -- subclassed (`CompMatchFit`), which adds
    one method, `rematch`, that swaps the matching set and the seed and
    recomputes the decoy distance field.  Everything else is inherited.
  * `run_phase3_nulls.fit_cell` -- called UNMODIFIED.  Its `n2` columns are
    the shared-lambda two-kernel sender-vs-decoy contrast; under `rematch`
    the decoys are composition-matched, so those columns become the
    composition-matched estimate.  They are renamed `*_matched` on the way
    out so nothing is mistaken for the published N2.
  * `sasp_estimators.BlockProfiler` -- the 400-replicate spatial block
    bootstrap, unchanged, via `fit_cell`.
  * `summarize_phase3.sf_table`'s reportable-population rule, reimplemented in
    two lines rather than imported, because importing that module reads
    `main_fits.csv`, which a concurrent job is rewriting.

--------------------------------------------------------------------------
ARMS
--------------------------------------------------------------------------
The runner is arm-generic: an arm is a section list, a cache directory, a
label set and the two Tier A calls.  M1 runs now.  H1 IS BEHIND THE s15
PRE-REGISTRATION FREEZE and this file REFUSES to touch it: `--arm h1` exits
unless the H1 arm has been populated after the freeze AND
SASP_H1_UNFROZEN=1 is set in the environment.  Nothing here reads
`/workspace/data/raw_h1/`.

Usage
  # the full protocol is two stages plus a merge:
  python3 -u run_phase8_compmatch.py --arm m1 --n-jobs 8 --calls tierA_p95 \
      --variants comp,full,comp_adj,type_adj,typecomp_adj --out-tag _tierA
  python3 -u run_phase8_compmatch.py --arm m1 --n-jobs 8 --calls tierApm_p95 \
      --variants comp,comp_adj,type_adj,typecomp_adj --out-tag _tierApm
  python3 -u run_phase8_compmatch.py --merge \
      results/phase3/compmatch_fits_tierA.csv,results/phase3/compmatch_fits_tierApm.csv

  python3 -u run_phase8_compmatch.py --arm m1 --n-jobs 4
  python3 -u run_phase8_compmatch.py --arm m1 --n-jobs 2 --sections 7259_liver_sbr_Male_26-U1 \
      --calls tierA_p95 --modules il6_jak_stat3 --seeds 20260901 --out-tag smoke

Outputs (never overwrites anything that existed before this file):
  results/phase3/compmatch_reruns<tag>.csv  -- per-seed rows + summary rows
  results/phase3/compmatch_fits<tag>.csv    -- the fit-level audit trail
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, "/workspace/code")
import sasp_estimators as E
import sasp_phase3 as P
import phase3_core as C
import run_phase3_nulls as RN

RES = P.RESULTS

# --- D15.3 the five literal pre-registered seeds --------------------------
COMPMATCH_SEEDS = (20260901, 20260902, 20260903, 20260904, 20260905)

# --- D15.2 the matching sets ----------------------------------------------
# `comp` is the protocol.  `full` is the published N2 matching set, rerun at
# the same seeds so that "composition alone" can be read against it.
VARIANTS = ("comp", "full", "comp_adj", "type_adj", "typecomp_adj")

CALIPER_SD = 0.25          # as in phase3_core.match_decoys_section
SMD_GATE = 0.10            # Master Plan s8 Test 5 / roadmap 9.3


# --------------------------------------------------------------------------
# arms
# --------------------------------------------------------------------------
ARMS = {
    "m1": dict(
        label="M1 mouse liver (GSE310392, Xenium Prime Mouse 5K)",
        cache=P.CACHE3,
        sections=list(P.IN_BAND),            # D15.7
        labels=P.LABELS,
        primary_call="tierA_p95",            # A_SENDER_FINAL_strict, PRIMARY
        permodule_call="tierApm_p95",        # A_sender_for_<module>.txt
        modules=list(P.MODULES),
        frozen=False,
    ),
    # H1.  Populated in Phase 10 (roadmap item 10.2), after Phase 9's A2/A5 gates
    # passed.  It stays `frozen=True`, so SASP_H1_UNFROZEN=1 is still required and
    # running it remains a deliberate act.  The cache is the one Phase 9 built
    # (`code/h1_prep_cache.py`) and extended in Phase 10 (`code/h1_cache_extend.py`).
    # PRIMARY CALL: `tierAmg_p95`, the identical Tier A percentile rule at the MERGED
    # label family -- the family `sasp_phase3.LABELS = "merged"` actually stratifies
    # receivers on.  Declared post-freeze PI decision; deviation H5,
    # `reports/CS_PHASE9_H1_AUDIT.md` §9.4.  The frozen-literal fine-label
    # `tierA_p95` is run alongside as the sensitivity, via `--calls`.
    "h1": dict(
        label="H1 human spleen (GSE326743, Xenium Prime Human 5K + 100 addon)",
        cache="/workspace/data/processed_h1/cache3_h1",
        sections=["SPLN07", "SPLN14", "SPLN21", "SPLN24",
                  "SPLN30", "SPLN43", "SPLN44"],
        labels=P.LABELS,
        primary_call="tierAmg_p95",
        permodule_call="tierApm_p95",
        modules=list(P.MODULES),
        frozen=True,
    ),
}


def arm_config(arm: str) -> dict:
    if arm not in ARMS:
        raise SystemExit(f"unknown arm {arm!r}; known: {sorted(ARMS)}")
    cfg = dict(ARMS[arm])
    if cfg["frozen"]:
        if os.environ.get("SASP_H1_UNFROZEN") != "1":
            raise SystemExit(
                f"REFUSING to run arm {arm!r}.\n"
                "  Phase7 s15: H1 is behind the pre-registration freeze and "
                "Phase 8's tag is not cut.\n"
                "  The protocol is implemented and arm-generic; running it on "
                "H1 is roadmap item 10.2, after Phase 9.\n"
                "  To run it after the freeze: populate ARMS['h1'] "
                "(cache + sections) and set SASP_H1_UNFROZEN=1.")
        if not cfg["sections"]:
            raise SystemExit(
                f"arm {arm!r} unfrozen but ARMS['{arm}']['sections'] is empty: "
                "populate it from Job B before running.")
    return cfg


# --------------------------------------------------------------------------
# the protocol
# --------------------------------------------------------------------------
def comp_col_idx(blocks: dict) -> np.ndarray:
    """Columns of `blocks['Zmatch']` that ARE the k-NN cell-type composition."""
    cols = blocks["zmatch_cols"]
    ix = [i for i, c in enumerate(cols) if c.startswith("knn_frac_")]
    if not ix:
        raise RuntimeError("no knn_frac_* columns in Zmatch; "
                           "phase3_core.build_blocks changed shape")
    return np.asarray(ix, dtype=int)


class CompMatchFit(RN.SectionFit):
    """`run_phase3_nulls.SectionFit` plus a swappable matching set.

    The parent already builds the covariate blocks, the neighbour baseline,
    the block ids, the lambda grid and one N2 match.  `rematch` is the only
    addition: it re-runs the SAME matcher with a different matching set and a
    different seed, and refreshes the decoy distance field that `fit_cell`
    reads.  No estimator code is duplicated.
    """

    def __init__(self, sample, call, seed, types=None, labels=None, module=None):
        super().__init__(sample, call, seed, types=types, labels=labels,
                         module=module)
        self._Zfull = self.blocks["Zmatch"]
        self._comp_ix = comp_col_idx(self.blocks)
        self.match_variant = "full"
        self.match_seed = seed

    def matching_set(self, variant: str):
        if variant == "full":
            return self._Zfull, list(self.blocks["zmatch_cols"])
        if variant == "comp":
            cols = [self.blocks["zmatch_cols"][i] for i in self._comp_ix]
            return self._Zfull[:, self._comp_ix], cols
        raise ValueError(variant)

    def rematch(self, variant: str, seed: int):
        Z, cols = self.matching_set(variant)
        self.match = C.match_decoys_section(self.sec, self.sender, Z, seed,
                                            caliper_sd=CALIPER_SD)
        dec = np.zeros(self.sec.n, bool)
        dec[self.match["decoy_idx"]] = True
        self.decoy = dec
        self.d_dec = (P.dist_to_senders(self.coords, dec) if dec.sum() > 10
                      else np.full(self.sec.n, np.nan))
        self.match_variant = variant
        self.match_seed = int(seed)
        self.match_cols = cols
        return self


# `fit_cell` writes the sender-vs-decoy contrast under `n2`; here the decoys
# are composition-matched, so rename on the way out.  Nothing may be confused
# with the published N2 column of `main_fits.csv`.
RENAME = {
    "beta_n2": "beta_matched", "sf_n2": "sf_matched",
    "beta_n2_lo": "beta_matched_lo", "beta_n2_hi": "beta_matched_hi",
    "sf_n2_lo": "sf_matched_lo", "sf_n2_hi": "sf_matched_hi",
    "lam_n2_profiled": "lam_matched_profiled",
    "beta_n2_profiled": "beta_matched_profiled",
    "beta_n2n5n6": "beta_matched_n5n6", "sf_n2n5n6": "sf_matched_n5n6",
    "beta_n2n5n6_lo": "beta_matched_n5n6_lo",
    "beta_n2n5n6_hi": "beta_matched_n5n6_hi",
    "sf_n2n5n6_lo": "sf_matched_n5n6_lo",
    "sf_n2n5n6_hi": "sf_matched_n5n6_hi",
    "beta_n2zon": "beta_matched_anat", "sf_n2zon": "sf_matched_anat",
    "beta_n2zon_lo": "beta_matched_anat_lo",
    "beta_n2zon_hi": "beta_matched_anat_hi",
    "sf_n2zon_lo": "sf_matched_anat_lo", "sf_n2zon_hi": "sf_matched_anat_hi",
}


def _adjust_designs(sf, ii):
    """Covariate blocks for the regression counterparts, from the SAME
    `phase3_core.build_blocks` output the matcher uses.

      comp_adj      20-NN cell-type composition (`knn_frac_*`) -- the N5 `comp`
                    sub-block of `run_phase3_attribution.py`, and exactly the
                    variables the `comp` matching set matches on.
      type_adj      receiver's OWN cell-type intercepts -- the mu_{c_i} of
                    Master Plan s6.1.  This is the block CS_PHASE3 s3.1
                    measured at SF 0.344, i.e. the "66 %" of the s17 row.
                    Degenerate inside a single-cell-type fit, so pooled only.
      typecomp_adj  both.
    """
    n5, Z5 = sf.blocks["n5_cols"], sf.blocks["N5"]
    ix = [i for i, c in enumerate(n5) if c.startswith("knn_frac_")]
    comp = Z5[np.ix_(ii, ix)]
    ct = sf.sec.celltype[ii]
    lev = np.array(sorted(set(ct)))[1:]
    dm = (np.column_stack([(ct == l).astype(float) for l in lev])
          if lev.size else np.zeros((ii.size, 0)))
    comp_names = "|".join(n5[i] for i in ix)
    type_names = "|".join(f"ct_{l}" for l in lev)
    return {
        "comp_adj": (comp, comp_names, "any"),
        "type_adj": (dm, type_names, "pooled"),
        "typecomp_adj": (np.column_stack([dm, comp]),
                         type_names + "|" + comp_names, "pooled"),
    }


ADJUSTMENTS = ("comp_adj", "type_adj", "typecomp_adj")


def adjust_rows(sf, scopes, js, arm, call, module, bootstrap_seed,
                which=ADJUSTMENTS):
    """D15.9 -- the REGRESSION counterpart of the matched rerun.

    CS_PHASE3 s5 N2 is explicit that in this dataset MATCHING is nearly inert
    (SF 0.943) while REGRESSING ON THE SAME COVARIATES removes 92 %: "matching
    balances the covariates between senders and decoys; it does not remove the
    dependence of the response on those covariates at the receiver, which is
    where the confounding acts."  A matched rerun on its own therefore answers
    "does composition MATCHING remove the kernel", not "does composition
    explain the kernel".  This computes the second question at the same scopes
    and the same lambda_hat_naive, with the same 400-replicate block
    bootstrap.  Deterministic -- no matching draw -- so these rows carry
    seed = -1.

    One `BlockProfiler` per adjustment gives BOTH the naive fit (p = 1, the
    intercept-only design) and the adjusted fit (p = full) for free, because
    `_acc(m, p)` uses `X[:, :p]`.
    """
    rows = []
    for t in scopes:
        idx = sf.receivers(t)
        ii = np.flatnonzero(idx)
        if ii.size < RN.MIN_RECEIVERS:
            continue
        kind = "pooled" if t is None else "by_celltype"
        bid = sf.bid[ii]
        ub, bid = np.unique(bid, return_inverse=True)
        nb = ub.size
        one = np.ones(nb)
        blocks = _adjust_designs(sf, ii)
        for name in which:
            Zb, cols, where = blocks[name]
            if where == "pooled" and kind != "pooled":
                continue
            if Zb.shape[1] == 0:
                continue
            X = np.column_stack([np.ones(ii.size), Zb])
            p_full = X.shape[1]
            for j in js:
                y = sf.Y[ii, j].astype(float)
                pr = E.BlockProfiler(sf.d_obs[ii], None, y, X, bid, nb,
                                     sf.lam, sf.lam)
                b0 = pr.fit1(one, 1)
                t0 = b0["t"]
                if not np.isfinite(b0["beta"]) or b0["beta"] == 0:
                    continue
                b_adj = pr.beta_at(one, p_full, t0)[0]
                rng = np.random.default_rng(int(bootstrap_seed) + 7 * j)
                bb0 = np.full(RN.N_BOOT, np.nan)
                bba = np.full(RN.N_BOOT, np.nan)
                for r in range(RN.N_BOOT):
                    m = rng.multinomial(nb, np.full(nb, 1.0 / nb)).astype(float)
                    try:
                        bb0[r] = pr.beta_at(m, 1, t0)[0]
                        bba[r] = pr.beta_at(m, p_full, t0)[0]
                    except Exception:
                        continue
                with np.errstate(divide="ignore", invalid="ignore"):
                    sfb = bba / bb0
                sfb = sfb[np.isfinite(sfb)]
                v0 = bb0[np.isfinite(bb0)]
                va = bba[np.isfinite(bba)]
                rows.append(dict(
                    section=sf.sec.name, arm=sf.sec.meta["condition"],
                    week=sf.sec.meta["week"], call=call,
                    celltype=(None if t is None else t),
                    module=P.MODULES[j], stratum="all", n=int(ii.size),
                    sender_set=(f"A_sender_for_{module}" if module
                                else "A_SENDER_FINAL_strict"),
                    n_senders=int(sf.sender.sum()),
                    prevalence=float(sf.sender.mean()),
                    lam_naive=b0["lam"], beta_naive=b0["beta"],
                    lam_grid_lo=float(sf.lam[0]),
                    lam_grid_hi=float(sf.lam[-1]),
                    lam_railed=int(t0 == 0 or t0 == sf.lam.size - 1),
                    sd_y=float(y.std()),
                    beta_base=b0["beta"], sf_base=1.0,
                    beta_n2=b_adj, sf_n2=(b_adj / b0["beta"]),
                    beta_base_lo=(float(np.quantile(v0, .025))
                                  if v0.size > 20 else np.nan),
                    beta_base_hi=(float(np.quantile(v0, .975))
                                  if v0.size > 20 else np.nan),
                    beta_n2_lo=(float(np.quantile(va, .025))
                                if va.size > 20 else np.nan),
                    beta_n2_hi=(float(np.quantile(va, .975))
                                if va.size > 20 else np.nan),
                    sf_n2_lo=(float(np.quantile(sfb, .025))
                              if sfb.size > 20 else np.nan),
                    sf_n2_hi=(float(np.quantile(sfb, .975))
                              if sfb.size > 20 else np.nan),
                    match_rate=np.nan, max_smd_before=np.nan,
                    max_smd_after=np.nan,
                    arm_id=arm, variant=name,
                    control_type="covariate_adjustment",
                    seed=-1, seed_index=-1,
                    scope=("ALL" if t is None else t), scope_kind=kind,
                    matched_on=cols, n_match_cols=int(Zb.shape[1]),
                    n_decoys=0, caliper_sd=np.nan,
                    n_boot=RN.N_BOOT, n_blocks=RN.N_BLOCKS_SIDE ** 2,
                    window_um=RN.WINDOW_UM))
    return rows


def _job(arm, cfg, sample, call, module, variants, seeds):
    """One (section, call, sender-module) unit: build the section ONCE, then
    rematch and refit for every (variant, seed, receiver scope, module)."""
    t0 = time.time()
    sf = CompMatchFit(sample, call, seeds[0], module=module)
    types = [t for t in sorted(set(sf.sec.celltype)) if t not in P.EXCLUDE_TYPES]
    scopes = [None] + [t for t in types
                       if int(sf.receivers(t).sum()) >= RN.MIN_RECEIVERS]
    js = RN._js(module)
    rows = []
    adj = [v for v in variants if v in ADJUSTMENTS]
    if adj:
        rows += adjust_rows(sf, scopes, js, arm, call, module,
                            bootstrap_seed=seeds[0], which=adj)
    for variant in [v for v in variants if v not in ADJUSTMENTS]:
        for si, seed in enumerate(seeds):
            sf.rematch(variant, seed)
            n_dec = int(sf.decoy.sum())
            for t in scopes:
                for j in js:
                    r = RN.fit_cell(sf, t, j, int(seed) + 7 * j, tag="all")
                    r["arm_id"] = arm
                    r["variant"] = variant
                    r["control_type"] = "decoy_matching"
                    r["seed"] = int(seed)
                    r["seed_index"] = si
                    r["scope"] = "ALL" if t is None else t
                    r["scope_kind"] = "pooled" if t is None else "by_celltype"
                    r["matched_on"] = "|".join(sf.match_cols)
                    r["n_match_cols"] = len(sf.match_cols)
                    r["n_decoys"] = n_dec
                    r["caliper_sd"] = CALIPER_SD
                    r["n_boot"] = RN.N_BOOT
                    r["n_blocks"] = RN.N_BLOCKS_SIDE ** 2
                    r["window_um"] = RN.WINDOW_UM
                    rows.append(r)
    print(f"[compmatch] {arm} {sample} {call} {module or '-'} "
          f"{len(rows)} rows {time.time()-t0:.0f}s", flush=True)
    return rows


def build_jobs(arm, cfg, sections, calls, modules, variants, seeds):
    jobs = []
    for s in sections:
        for c in calls:
            if RN.is_permodule(c):
                for m in modules:
                    jobs.append((arm, cfg, s, c, m, variants, seeds))
            else:
                jobs.append((arm, cfg, s, c, None, variants, seeds))
    return jobs


# --------------------------------------------------------------------------
# summary
# --------------------------------------------------------------------------
def _q(v, q):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.quantile(v, q)) if v.size else np.nan


def _nanmax(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(v.max()) if v.size else np.nan


def _gate(v):
    """Fraction of matches meeting the Master Plan s8 Test 5 balance gate.
    NaN for the covariate-adjustment variants, which do no matching."""
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.mean(v <= SMD_GATE)) if v.size else np.nan


def _med(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else np.nan


def summarise(fits: pd.DataFrame) -> pd.DataFrame:
    """Per-seed rows, then one across-seed summary row per stratification.

    Reportable population = summarize_phase3.sf_table's rule: a positive naive
    amplitude whose block-bootstrap CI excludes zero.
    """
    d = fits.copy()
    for c in ("beta_naive", "beta_base_lo", "sd_y", "sf_matched",
              "sf_matched_lo", "sf_matched_hi", "sf_matched_n5n6", "sf_n5",
              "sf_n6n5", "beta_matched", "lam_naive", "lam_railed",
              "match_rate", "max_smd_before", "max_smd_after",
              "control_type"):
        if c not in d.columns:
            d[c] = np.nan
    d["reportable"] = (d.beta_naive > 0) & (d.beta_base_lo > 0)
    d["beta_naive_sd"] = d.beta_naive / d.sd_y
    d["comp_share"] = 1.0 - d.sf_matched      # share of the naive amplitude a
                                              # matched decoy already reproduces
    keys = ["arm_id", "call", "sender_set_kind", "variant",
            "control_type", "scope_kind"]
    d["sender_set_kind"] = np.where(d.call.map(RN.is_permodule),
                                    "A_sender_for_<module> (7 sets)",
                                    "A_SENDER_FINAL_strict")
    rows = []
    for k, g in d.groupby(keys, dropna=False):
        rep_all = g[g.reportable]
        per_seed = {}
        for seed, gs in g.groupby("seed"):
            rep = gs[gs.reportable]
            v = rep.sf_matched.to_numpy(float)
            row = dict(zip(keys, k))
            row.update(
                row_type="seed", seed=int(seed),
                seed_index=int(gs.seed_index.iloc[0]),
                n_fits=int(len(gs)), n_reportable=int(len(rep)),
                frac_reportable=float(len(rep) / max(len(gs), 1)),
                median_sf_matched=_med(v),
                sf_matched_q25=_q(v, .25), sf_matched_q75=_q(v, .75),
                frac_sf_le_0=float(np.mean(v[np.isfinite(v)] <= 0))
                if np.isfinite(v).any() else np.nan,
                frac_sf_gt_05=float(np.mean(v[np.isfinite(v)] > .5))
                if np.isfinite(v).any() else np.nan,
                median_comp_share=_med(rep.comp_share),
                median_sf_matched_lo=_med(rep.sf_matched_lo),
                median_sf_matched_hi=_med(rep.sf_matched_hi),
                median_sf_matched_n5n6=_med(rep.sf_matched_n5n6),
                median_sf_n5=_med(rep.sf_n5),
                median_sf_n6n5=_med(rep.sf_n6n5),
                median_beta_naive=_med(rep.beta_naive),
                median_beta_naive_sd=_med(rep.beta_naive_sd),
                median_beta_matched=_med(rep.beta_matched),
                median_lam_naive=_med(rep.lam_naive),
                frac_lam_railed=float(np.nanmean(rep.lam_railed))
                if len(rep) else np.nan,
                median_match_rate=_med(gs.match_rate),
                median_max_smd_before=_med(gs.max_smd_before),
                median_max_smd_after=_med(gs.max_smd_after),
                frac_smd_gate_pass=_gate(gs.max_smd_after),
                n_match_cols=int(gs.n_match_cols.iloc[0]),
                matched_on=gs.matched_on.iloc[0],
                n_sections=int(gs.section.nunique()),
                n_modules=int(gs.module.nunique()),
                n_boot=int(gs.n_boot.iloc[0]),
                n_blocks=int(gs.n_blocks.iloc[0]),
                window_um=float(gs.window_um.iloc[0]),
                caliper_sd=float(gs.caliper_sd.iloc[0]),
            )
            per_seed[int(seed)] = row
            rows.append(row)
        ps = pd.DataFrame(list(per_seed.values()))
        v_all = rep_all.sf_matched.to_numpy(float)
        s = dict(zip(keys, k))
        s.update(
            row_type="summary", seed=-1, seed_index=-1,
            n_seeds=int(ps.shape[0]),
            seeds="|".join(str(x) for x in sorted(per_seed)),
            n_fits=int(len(g)),
            n_reportable=int(len(rep_all)),
            frac_reportable=float(len(rep_all) / max(len(g), 1)),
            # across-seed behaviour of the per-seed medians -- this is the
            # "at 5 seeds" readout
            median_sf_matched=float(np.nanmedian(ps.median_sf_matched)),
            sf_across_seed_mean=float(np.nanmean(ps.median_sf_matched)),
            sf_across_seed_sd=float(np.nanstd(ps.median_sf_matched, ddof=1))
            if ps.shape[0] > 1 else np.nan,
            sf_across_seed_min=float(np.nanmin(ps.median_sf_matched)),
            sf_across_seed_max=float(np.nanmax(ps.median_sf_matched)),
            sf_across_seed_range=float(np.nanmax(ps.median_sf_matched)
                                       - np.nanmin(ps.median_sf_matched)),
            # pooled over every seed's fits (the number to quote in s17)
            sf_matched_q25=_q(v_all, .25), sf_matched_q75=_q(v_all, .75),
            frac_sf_le_0=float(np.mean(v_all[np.isfinite(v_all)] <= 0))
            if np.isfinite(v_all).any() else np.nan,
            frac_sf_gt_05=float(np.mean(v_all[np.isfinite(v_all)] > .5))
            if np.isfinite(v_all).any() else np.nan,
            median_comp_share=_med(rep_all.comp_share),
            comp_share_q25=_q(rep_all.comp_share, .25),
            comp_share_q75=_q(rep_all.comp_share, .75),
            comp_share_across_seed_min=float(np.nanmin(ps.median_comp_share)),
            comp_share_across_seed_max=float(np.nanmax(ps.median_comp_share)),
            median_sf_matched_lo=_med(rep_all.sf_matched_lo),
            median_sf_matched_hi=_med(rep_all.sf_matched_hi),
            median_sf_matched_n5n6=_med(rep_all.sf_matched_n5n6),
            median_sf_n5=_med(rep_all.sf_n5),
            median_sf_n6n5=_med(rep_all.sf_n6n5),
            median_beta_naive=_med(rep_all.beta_naive),
            median_beta_naive_sd=_med(rep_all.beta_naive_sd),
            median_beta_matched=_med(rep_all.beta_matched),
            median_lam_naive=_med(rep_all.lam_naive),
            frac_lam_railed=float(np.nanmean(rep_all.lam_railed))
            if len(rep_all) else np.nan,
            median_match_rate=_med(g.match_rate),
            median_max_smd_before=_med(g.max_smd_before),
            median_max_smd_after=_med(g.max_smd_after),
            max_max_smd_after=_nanmax(g.max_smd_after),
            frac_smd_gate_pass=_gate(g.max_smd_after),
            n_match_cols=int(g.n_match_cols.iloc[0]),
            matched_on=g.matched_on.iloc[0],
            n_sections=int(g.section.nunique()),
            n_modules=int(g.module.nunique()),
            n_boot=int(g.n_boot.iloc[0]),
            n_blocks=int(g.n_blocks.iloc[0]),
            window_um=float(g.window_um.iloc[0]),
            caliper_sd=float(g.caliper_sd.iloc[0]),
        )
        rows.append(s)
    out = pd.DataFrame(rows)
    lead = ["row_type", "arm_id", "call", "sender_set_kind", "variant",
            "control_type", "scope_kind", "seed", "seed_index"]
    return out[lead + [c for c in out.columns if c not in lead]]


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="m1", choices=sorted(ARMS))
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument("--sections", default=None,
                    help="comma list; default = the arm's admissible sections")
    ap.add_argument("--calls", default=None,
                    help="comma list; default = the arm's two Tier A variants")
    ap.add_argument("--modules", default=None,
                    help="comma list; default = all seven")
    ap.add_argument("--variants", default=",".join(VARIANTS))
    ap.add_argument("--seeds", default=",".join(str(s) for s in COMPMATCH_SEEDS))
    ap.add_argument("--out-tag", default="")
    ap.add_argument("--merge", default=None,
                    help="comma list of existing compmatch_fits*.csv to "
                         "concatenate and re-summarise into --out-tag; runs no "
                         "fits and reads no cell data")
    a = ap.parse_args()

    if a.merge:
        parts = [pd.read_csv(x) for x in a.merge.split(",")]
        fits = pd.concat(parts, ignore_index=True)
        fp = f"{RES}/compmatch_fits{a.out_tag}.csv"
        sp = f"{RES}/compmatch_reruns{a.out_tag}.csv"
        for p_ in (fp, sp):
            if os.path.exists(p_):
                raise SystemExit(f"{p_} exists; refusing to overwrite.")
        fits.to_csv(fp, index=False)
        summ = summarise(fits)
        summ.to_csv(sp, index=False)
        print(f"merged {len(parts)} files -> {fp} {fits.shape}")
        print(f"wrote {sp} {summ.shape}")
        show = ["row_type", "call", "variant", "scope_kind", "seed",
                "n_reportable", "median_sf_matched", "median_comp_share",
                "median_match_rate", "median_max_smd_after"]
        with pd.option_context("display.width", 220):
            print(summ[show].to_string(index=False))
        return

    cfg = arm_config(a.arm)
    sections = (a.sections.split(",") if a.sections else cfg["sections"])
    sections = [s for s in sections
                if os.path.exists(os.path.join(cfg["cache"], f"{s}.npz"))]
    if not sections:
        raise SystemExit("no cached sections for this arm")
    calls = (a.calls.split(",") if a.calls
             else [cfg["primary_call"], cfg["permodule_call"]])
    modules = a.modules.split(",") if a.modules else cfg["modules"]
    variants = a.variants.split(",")
    seeds = [int(s) for s in a.seeds.split(",")]

    print(f"arm       : {a.arm}  ({cfg['label']})")
    print(f"sections  : {sections}")
    print(f"calls     : {calls}")
    print(f"modules   : {modules}")
    print(f"variants  : {variants}")
    print(f"seeds     : {seeds}", flush=True)

    jobs = build_jobs(a.arm, cfg, sections, calls, modules, variants, seeds)
    t0 = time.time()
    out = Parallel(n_jobs=a.n_jobs, prefer="processes", verbose=5)(
        delayed(_job)(*j) for j in jobs)
    fits = pd.DataFrame([r for rs in out for r in rs]).rename(columns=RENAME)
    tag = a.out_tag
    fp = f"{RES}/compmatch_fits{tag}.csv"
    sp = f"{RES}/compmatch_reruns{tag}.csv"
    for p in (fp, sp):
        if os.path.exists(p):
            raise SystemExit(f"{p} exists; refusing to overwrite. Use --out-tag")
    fits.to_csv(fp, index=False)
    summ = summarise(fits)
    summ.to_csv(sp, index=False)
    print(f"\nwrote {fp}  {fits.shape}")
    print(f"wrote {sp}  {summ.shape}")
    print(f"elapsed {time.time()-t0:.0f}s")
    show = ["row_type", "call", "variant", "scope_kind", "seed", "n_reportable",
            "median_sf_matched", "median_comp_share", "median_match_rate",
            "median_max_smd_after"]
    with pd.option_context("display.width", 200):
        print(summ[show].to_string(index=False))


if __name__ == "__main__":
    main()
