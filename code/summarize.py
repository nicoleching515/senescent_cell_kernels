"""Aggregate the sweep into the tables that go in the phase report."""
import numpy as np, pandas as pd
pd.set_option("display.width", 220)

RES = "/workspace/results"
LAM = 30.0
EST = ["naive", "naive_bin", "nuis", "decoy", "decoyS", "decoyS_nuis"]


def ci95(x):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if x.size < 2: return (np.nan, np.nan)
    se = x.std(ddof=1) / np.sqrt(x.size)
    return (x.mean() - 1.96 * se, x.mean() + 1.96 * se)


def fmt_ci(x):
    lo, hi = ci95(x); return f"{np.nanmean(x):.2f} [{lo:.2f}, {hi:.2f}]"


def main():
    df = pd.read_csv(f"{RES}/sweep_all.csv")
    # pandas parses the literal string "null" as NaN, which silently drops the
    # beta_true = 0 sweep.  Restore the label.
    df["sweep"] = df["sweep"].fillna("null")
    out = []
    W = out.append

    W("### 1. Easy / clean regime sanity check (does recovery work at all?)\n")
    cl = df[df.sweep == "clean"]
    rows = []
    for reg, g in cl.groupby("regime"):
        r = dict(regime=reg, n_seeds=len(g))
        for e in ["naive", "nuis", "decoyS"]:
            r[f"lam_{e}"] = fmt_ci(g[f"lam_{e}"])
            r[f"beta_{e}"] = fmt_ci(g[f"beta_{e}"])
        r["cover_naive_iid"] = round(g.cover_lam_naive_iid.mean(), 2)
        r["cover_naive_blk"] = round(g.cover_lam_naive_blk.mean(), 2)
        r["maxSMD_after"] = round(g.max_smd_after.mean(), 3)
        rows.append(r)
    W(pd.DataFrame(rows).to_string(index=False))

    W("\n\n### 2. Main grid: clustering x autocorrelation length\n")
    m = df[df.sweep == "main"]
    for col, name in [("lam_naive", "lambda_hat naive"),
                      ("lam_decoyS", "lambda_hat matched-decoy"),
                      ("lam_nuis", "lambda_hat nuisance-cond"),
                      ("cover_lam_naive_iid", "coverage naive iid CI"),
                      ("cover_lam_naive_blk", "coverage naive block CI"),
                      ("cover_lam_decoyS_blk", "coverage decoy block CI"),
                      ("beta_naive", "beta_hat naive"),
                      ("beta_decoyS", "beta_hat matched-decoy"),
                      ("se_ratio_naive", "block sd / iid SE"),
                      ("max_smd_after", "max SMD after matching"),
                      ("ripley50", "realized sender Ripley ratio @50um")]:
        p = m.pivot_table(index="clustering", columns="ell_over_lambda",
                          values=col, aggfunc="mean").round(2)
        W(f"\n{name}  (rows: kappa, cols: ell/lambda_true)\n{p.to_string()}")

    W("\n\n### 3. Sender prevalence (Section 8 Test 3)\n")
    pv = df[df.sweep == "prev"]
    W(pv.groupby("prevalence").agg(
        n_senders=("n_senders", "mean"), med_d=("med_d_sender", "mean"),
        lam_naive=("lam_naive", "mean"), lam_nuis=("lam_nuis", "mean"),
        lam_decoyS=("lam_decoyS", "mean"),
        sd_lam_naive=("lam_naive", "std"),
        cov_iid=("cover_lam_naive_iid", "mean"),
        cov_blk=("cover_lam_naive_blk", "mean"),
        b_naive=("beta_naive", "mean")).round(2).to_string())

    W("\n\n### 4. Confounder strength\n")
    cf = df[df.sweep == "conf"]
    W(cf.groupby("conf_strength").agg(
        lam_naive=("lam_naive", "mean"), lam_decoyS=("lam_decoyS", "mean"),
        lam_nuis=("lam_nuis", "mean"), b_naive=("beta_naive", "mean"),
        b_decoyS=("beta_decoyS", "mean"), b_nuis=("beta_nuis", "mean"),
        cov_iid=("cover_lam_naive_iid", "mean"),
        u_s=("u_sender_mean", "mean"), u_d=("u_decoy_mean", "mean"),
        u_all=("u_all_mean", "mean"), smd=("max_smd_after", "mean")
    ).round(3).to_string())

    W("\n\n### 5. Null (beta_true = 0): does any estimator manufacture signal?\n")
    nl = df[df.sweep == "null"]
    W(nl.groupby(["clustering", "conf_strength", "autocorr_len_um"]).agg(
        b_naive=("beta_naive", "mean"), b_decoyS=("beta_decoyS", "mean"),
        b_nuis=("beta_nuis", "mean"), lam_naive=("lam_naive", "mean"),
        sd_b_naive=("beta_naive", "std"),
        frac_beta_gt_0p2=("beta_naive", lambda s: float((s.abs() > 0.2).mean())),
        cov_beta_naive=("cover_beta_naive_blk", "mean"),
        cov_beta_decoyS=("cover_beta_decoyS_blk", "mean"),
    ).round(3).to_string())

    W("\n\n### 6. N (window size at fixed density)\n")
    ns = df[df.sweep == "nsize"]
    W(ns.groupby("window_um").agg(
        n_cells=("n_cells", "mean"), lam_naive=("lam_naive", "mean"),
        sd=("lam_naive", "std"), cov_iid=("cover_lam_naive_iid", "mean"),
        cov_blk=("cover_lam_naive_blk", "mean")).round(2).to_string())

    W("\n\n### 7. Decoy control cost in TRUE signal (Section 29, objection 6)\n")
    rows = []
    for lbl, g in [("clean, no confounder", cl[cl.regime == "easy"]),
                   ("clean, no nuisance at all", cl[cl.regime == "pure"]),
                   ("ell/lam = 0.25", m[np.isclose(m.ell_over_lambda, .25)]),
                   ("ell/lam = 1", m[np.isclose(m.ell_over_lambda, 1.)]),
                   ("ell/lam = 4", m[np.isclose(m.ell_over_lambda, 4.)])]:
        rows.append(dict(
            regime=lbl,
            beta_naive=fmt_ci(g.beta_naive),
            beta_decoyS=fmt_ci(g.beta_decoyS),
            retained_vs_naive=f"{np.nanmean(g.beta_decoyS / g.beta_naive):.3f}",
            retained_vs_truth=fmt_ci(g.signal_retained_decoyS),
            plan_literal_bt_minus_bd=fmt_ci(g.beta_true_minus_decoy),
            u_gap_closed=f"{np.nanmean((g.u_decoy_mean-g.u_all_mean)/(g.u_sender_mean-g.u_all_mean)):.2f}"))
    W(pd.DataFrame(rows).to_string(index=False))

    txt = "\n".join(out)
    open(f"{RES}/summary_tables.txt", "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
