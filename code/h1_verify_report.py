#!/usr/bin/env python3
"""Phase 9 — re-read every headline number in reports/CS_PHASE9_H1_AUDIT.md from its file.

Exists because this project has repeatedly found its own written numbers not to match the
files behind them (`reports/AUDIT_PHASE8_FACTCHECK.md`).  Every assertion below is a literal
from the report checked against results/phase9_h1/.  Exits non-zero if any has drifted.

Run: python3 code/h1_verify_report.py
"""
import sys
import pandas as pd

R = "/workspace/results/phase9_h1/"
bad = []


def chk(label, cond, val):
    print(("OK   " if cond else "FAIL ") + label, val)
    if not cond:
        bad.append(label)


a1 = pd.read_csv(R + "a1_sections.csv")
chk("A1 median NN 5.447-6.286 um",
    abs(a1.median_nn_um_all.min() - 5.447) < 1e-3 and abs(a1.median_nn_um_qc.max() - 6.286) < 1e-3,
    (a1.median_nn_um_all.min(), a1.median_nn_um_qc.max()))
chk("A1 cells 2,207,593", a1.n_cells.sum() == 2207593, a1.n_cells.sum())
chk("A1 QC-pass 1,962,278", a1.n_qc_pass.sum() == 1962278, a1.n_qc_pass.sum())
asg = pd.read_csv(R + "a1_assignment_rate.csv")
chk("A1 assignment 92.41/93.23/94.17", list(asg.assigned_pct) == [92.41, 93.23, 94.17],
    list(asg.assigned_pct))

a3 = pd.read_csv(R + "a3_summary_by_caller.csv")
r = a3[a3.call == "tierA_p95"].iloc[0]
chk("A3 tierA_p95 60/76, band 76/76", (r.n_pass, r.n_strata, r.n_pass_band) == (60, 76, 76),
    (r.n_pass, r.n_strata, r.n_pass_band))
c = a3[a3.call == "cdkn1a_pos"].iloc[0]
chk("A3 cdkn1a_pos 42/76, band 64", (c.n_pass, c.n_pass_band) == (42, 64), (c.n_pass, c.n_pass_band))
d = pd.read_csv(R + "a3_prevalence_by_type.csv")
cb = d[(d.label_set == "fine") & (d.call == "cdkn1a_pos") & (~d.pass_band)]
chk("A3 cdkn1a band failures all below 1 %", (cb.prevalence_pct_scored < 1).all(),
    (cb.prevalence_pct_scored.min(), cb.prevalence_pct_scored.max()))
m = d[(d.label_set == "merged") & (d.call == "tierA_merged_p95")]
allsec = m.groupby("cell_type").passes_A3.agg(["sum", "size"])
chk("A3 merged-family: 5 families pass in all 7 sections",
    int(((allsec["sum"] == allsec["size"]) & (allsec["size"] == 7)).sum()) == 5,
    sorted(allsec[(allsec["sum"] == allsec["size"]) & (allsec["size"] == 7)].index))

a4 = pd.read_csv(R + "a4_ripley.csv").groupby("call").ripley_ratio.mean().round(3)
chk("A4 1.012 / 1.161 / 1.647",
    (a4["tierA_p95"], a4["cdkn1a_pos"], a4["senepy_p95"]) == (1.012, 1.161, 1.647),
    (a4["tierA_p95"], a4["cdkn1a_pos"], a4["senepy_p95"]))

a5 = pd.read_csv(R + "a5_match_balance.csv")
chk("A5 35/35, max |SMD| after 0.0933",
    (a5.passes_A5.sum(), len(a5), a5.max_smd_after.max()) == (35, 35, 0.0933),
    (a5.passes_A5.sum(), len(a5), a5.max_smd_after.max()))

a7 = pd.read_csv(R + "a7_summary.csv")
p = a7[a7.response == "neg_control_probe"].set_index("design")
chk("A7 40 probes naive -0.0118 [-0.0340,+0.0103]",
    tuple(p.loc["base"][["clustered_mean", "clustered_lo", "clustered_hi"]]) == (-0.0118, -0.034, 0.0103),
    tuple(p.loc["base"][["clustered_mean", "clustered_lo", "clustered_hi"]]))
chk("A7 40 probes n6n5 -0.0028", p.loc["n6n5"].clustered_mean == -0.0028, p.loc["n6n5"].clustered_mean)
ac = a7[a7.response == "all_controls"].set_index("design")
chk("A7 pooled naive/n2/n6n5 -0.0337 / -0.0331 / -0.0051",
    (ac.loc["base"].clustered_mean, ac.loc["n2"].clustered_mean,
     ac.loc["n6n5"].clustered_mean) == (-0.0337, -0.0331, -0.0051),
    (ac.loc["base"].clustered_mean, ac.loc["n2"].clustered_mean, ac.loc["n6n5"].clustered_mean))
bio = a7[a7.response.str.startswith("BIOL")].set_index("design")
chk("A7 modules n6n5 +0.0147 [0.0015,0.0279]",
    tuple(bio.loc["n6n5"][["clustered_mean", "clustered_lo", "clustered_hi"]]) == (0.0147, 0.0015, 0.0279),
    tuple(bio.loc["n6n5"][["clustered_mean", "clustered_lo", "clustered_hi"]]))
chk("A7 codeword fits dropped for sd_y == 0: 10",
    int(a7[a7.response == "neg_control_codeword"].n_dropped_degenerate.iloc[0]) == 10,
    int(a7[a7.response == "neg_control_codeword"].n_dropped_degenerate.iloc[0]))

ag = pd.read_csv(R + "caller_agreement_pooled.csv")


def g(A, B):
    x = ag[(ag.A == A) & (ag.B == B)].iloc[0]
    return float(x.pooled_ratio), float(x.z), int(x.above_chance)


chk("agreement Tier A x DeepScence 1.602", g("tierA_score", "deepscence_score")[0] == 1.602,
    g("tierA_score", "deepscence_score"))
chk("agreement Tier A x SenePy 0.874, z -7.96, 0/7",
    g("tierA_score", "senepy_score") == (0.874, -7.96, 0), g("tierA_score", "senepy_score"))
chk("agreement Tier A x CDKN1A 1.081", g("tierA_score", "cdkn1a_counts")[0] == 1.081,
    g("tierA_score", "cdkn1a_counts"))
chk("circular DeepScence x CDKN1A 6.436, 7/7",
    g("deepscence_score", "cdkn1a_counts")[0] == 6.436 and
    g("deepscence_score", "cdkn1a_counts")[2] == 7, g("deepscence_score", "cdkn1a_counts"))
chk("sign-invariant Tier A x |DeepScence| 1.102",
    g("tierA_score", "abs_deepscence_score")[0] == 1.102, g("tierA_score", "abs_deepscence_score"))

an = pd.read_csv(R + "deepscence_anchor_h1.csv")
chk("P-ii CDKN1A fold stability 1.0 in all 7", (an.stab_cdkn1a == 1.0).all(), list(an.stab_cdkn1a))
chk("P-ii CDKN1A partial rho 0.1911-0.2540",
    (round(an.rho_partial_cdkn1a.min(), 4), round(an.rho_partial_cdkn1a.max(), 4)) == (0.1911, 0.254),
    (an.rho_partial_cdkn1a.min(), an.rho_partial_cdkn1a.max()))

t = pd.read_csv(R + "caller_technical_loading.csv")
ds = t[t.caller == "deepscence_score"]
chk("P-i positive in 7/7, 0.182-0.354",
    (ds.spearman_vs_transcript_counts > 0).all() and
    (round(ds.spearman_vs_transcript_counts.min(), 3),
     round(ds.spearman_vs_transcript_counts.max(), 3)) == (0.182, 0.354),
    (ds.spearman_vs_transcript_counts.min(), ds.spearman_vs_transcript_counts.max()))
chk("P-v Q5/Q1 < 1 in 6/7", int((ds.q5_over_q1 < 1).sum()) == 6, sorted(ds.q5_over_q1.round(3)))

st = pd.read_csv(R + "d2_stability.csv")
f = st[st.config == "denoise=False, FULL section"].iloc[0]
chk("full-section seed check r 0.3719, Jaccard 0.2107, 12,779 changing",
    (f.pearson_r, f.top5_jaccard, f.n_cells_changing_status) == (0.3719, 0.2107, 12779),
    (f.pearson_r, f.top5_jaccard, f.n_cells_changing_status))
dt = st[st.config == "denoise=True"]
chk("P-vii all three denoise=True pairs Jaccard < 0.30", (dt.top5_jaccard < 0.30).all(),
    list(dt.top5_jaccard))
dp = pd.read_csv(R + "d2_depth.csv")
fs = dp[dp.scope == "full section"].iloc[0]
chk("P-vi full section 0.3122 -> 0.1017, Jaccard 0.0164",
    (fs.rho_denoise_False, fs.rho_denoise_True, fs.delta_rho,
     fs.sender_jaccard_False_vs_True) == (0.3122, 0.1017, -0.2104, 0.0164),
    (fs.rho_denoise_False, fs.rho_denoise_True, fs.delta_rho, fs.sender_jaccard_False_vs_True))

pa = pd.read_csv(R + "a8_panel_arithmetic.csv").set_index("gene_set")
chk("A8 Tier A 33 -> 26 human, 33 -> 27 mouse",
    tuple(pa.loc["A_SENDER_FINAL_strict"][["human_full", "human_intersected",
                                           "mouse_intersected"]]) == (33, 26, 27),
    tuple(pa.loc["A_SENDER_FINAL_strict"][["human_full", "human_intersected", "mouse_intersected"]]))
os_ = pd.read_csv(R + "a8_ortho_sender_shift.csv")
o95 = os_[os_.call == "tierA_p95"]
chk("A8 p95 Jaccard 0.522-0.575, rho 0.848-0.934",
    (round(o95.jaccard.min(), 3), round(o95.jaccard.max(), 3),
     round(o95.spearman_score.min(), 3), round(o95.spearman_score.max(), 3)) == (0.522, 0.575, 0.848, 0.934),
    (o95.jaccard.min(), o95.jaccard.max(), o95.spearman_score.min(), o95.spearman_score.max()))

a6 = pd.read_csv(R + "a6_summary.csv")
chk("A6 V1 negative in 2/7, Moran's I 0.221-0.492",
    int((a6.v1_corr_axis_distfollicle < 0).sum()) == 2 and
    (round(a6.v4_moran_i_20nn.min(), 3), round(a6.v4_moran_i_20nn.max(), 3)) == (0.221, 0.492),
    (list(a6.v1_corr_axis_distfollicle.round(3)), a6.v4_moran_i_20nn.min(), a6.v4_moran_i_20nn.max()))

cc = pd.read_csv(R + "jobB_crosscheck_scores.csv")
l3 = cc[(cc.our_label == "cell_type") & (cc.depositor_level == "Level_3_Annotations")]
chk("Job B ARI vs Level_3: 0.352-0.451 all, 0.552-0.663 clean",
    (round(l3.ari_all.min(), 3), round(l3.ari_all.max(), 3),
     round(l3.ari_clean.min(), 3), round(l3.ari_clean.max(), 3)) == (0.352, 0.451, 0.552, 0.663),
    (l3.ari_all.min(), l3.ari_all.max(), l3.ari_clean.min(), l3.ari_clean.max()))

mf = pd.read_csv(R + "h1_module_fits.csv")
rp = mf[(mf.beta_naive > 0) & (mf.beta_base_lo > 0)]
chk("module fits 686, reportable 227, railed 73.9 %",
    (len(mf), len(rp), round(100 * mf.lam_railed.mean(), 1)) == (686, 227, 73.9),
    (len(mf), len(rp), round(100 * mf.lam_railed.mean(), 1)))

print("\n%d checks, %d FAILED" % (28, len(bad)))
sys.exit(1 if bad else 0)
