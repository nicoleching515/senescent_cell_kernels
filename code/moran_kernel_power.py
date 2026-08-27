#!/usr/bin/env python3
"""How much Moran's I can a distance-to-sender kernel of the size A7 measures
actually move?  The quantitative core of the §29 objection-9 answer.

The A7 model is  y = X gamma + beta * k,  k_i = exp(-d_i / lambda), d_i the
distance from cell i to the nearest sender.  If the technical gradient is the
only structure in y and the rest of y is spatially white, then

    I(y)  =  [ beta^2 Var(k) I(k) + Var(eps) I(eps) ] / Var(y)
          ~=  beta_z^2 * Var(k) * I(k)          (y z-scored, I(eps) ~ 0)

so the *entire* gradient A7 detects contributes a predictable, computable
amount to Moran's I.  This script computes that amount on the real sections,
using the real lambda-hat, the real k, the real I(k) and A7's own beta, and
inverts it two ways:

  dI_pred     what the measured A7 gradient contributes to Moran's I
  beta_min    the smallest distance-to-sender amplitude whose contribution to
              Moran's I would reach 2 SE(I) -- i.e. the smallest gradient a
              Moran's I test could see at all on this data

Writes results/moran/moran_kernel_power.csv.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/workspace/code")
import sasp_phase3 as P                       # noqa: E402
import run_phase3_nulls as R                  # noqa: E402
import run_moran_controls as MC               # noqa: E402

RES = "/workspace/results/phase3/"
CALLS = ["tierA_p95", "cdkn1a_pos"]
RESP = "all_controls"


def main():
    A = pd.read_csv(RES + "a7_control_probe_fits.csv")
    A["e_base"] = A.beta_base / A.sd_y
    rows = []
    for s in P.ALL_SECTIONS:
        sec = P.Sec(s)
        co = sec.coords.astype(float)
        cells = (pd.read_parquet(os.path.join(MC.RAW, s, "cells.parquet"))
                 .set_index("cell_id")
                 .reindex(pd.Index(sec.cell_id.astype(str))))
        y_all = (cells["control_probe_counts"] + cells["control_codeword_counts"]
                 + cells["genomic_control_counts"]).to_numpy(float)
        for call in CALLS:
            sub = A[(A.section == s) & (A.call == call) & (A.response == RESP)]
            if not len(sub):
                continue
            lam = float(np.median(sub.lam_naive))
            beta_z = float(sub.e_base.mean())
            snd = sec.sender_mask(call)
            d = P.dist_to_senders(co, snd)
            m = ((~np.isin(sec.celltype, P.EXCLUDE_TYPES)) & (~snd)
                 & np.isfinite(d) & (d <= R.WINDOW_UM))
            if m.sum() < R.MIN_RECEIVERS:
                continue
            k = np.exp(-d[m] / lam)
            y = y_all[m]
            y = (y - y.mean()) / (y.std() + 1e-12)
            W = MC.knn_weights(co[m], MC.KNN_PRIMARY)
            S0, S1, S2 = MC.w_moments(W)
            Y = np.column_stack([y, k])
            I = MC.moran_I(W, S0, Y)
            EI, zn, pn, zr, pr = MC.moran_inference(W, S0, S1, S2, Y, I)
            se_I = abs((I[0] - EI) / zr[0]) if np.isfinite(zr[0]) else np.nan
            vark = float(k.var())
            dI = beta_z ** 2 * vark * I[1]
            beta_min = float(np.sqrt(2 * se_I / (vark * I[1]))) if I[1] > 0 else np.nan
            rows.append(dict(
                section=s, call=call, n_receivers=int(m.sum()),
                lam_naive=lam, a7_beta_z=beta_z,
                I_controls=I[0], SE_I=se_I, I_kernel=I[1], var_k=vark,
                dI_from_a7_gradient=dI,
                dI_as_frac_of_I_controls=dI / I[0] if I[0] else np.nan,
                beta_min_visible_to_moran=beta_min))
            print("[power] %s %-11s lam=%4.1f beta_z=%+.4f I_ctrl=%+.5f "
                  "I_k=%+.4f dI=%.2e (%.3f%% of I) beta_min=%.3f"
                  % (s, call, lam, beta_z, I[0], I[1], dI,
                     100 * dI / I[0] if I[0] else np.nan, beta_min), flush=True)
    D = pd.DataFrame(rows)
    os.makedirs(MC.OUT, exist_ok=True)
    D.to_csv(os.path.join(MC.OUT, "moran_kernel_power.csv"), index=False)
    print()
    print("MEDIANS over %d section x call cells:" % len(D))
    for c in ["lam_naive", "a7_beta_z", "I_controls", "SE_I", "I_kernel",
              "var_k", "dI_from_a7_gradient", "dI_as_frac_of_I_controls",
              "beta_min_visible_to_moran"]:
        print("  %-28s %.6g" % (c, D[c].median()))


if __name__ == "__main__":
    main()
