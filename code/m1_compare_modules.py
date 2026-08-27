#!/usr/bin/env python3
"""Phase 8 / 8.7 -- per-module and per-cell-type pre/post comparison of the
Phase 3 main fits, so movement can be attributed to the gene sets that changed.

  python3 m1_compare_modules.py OLD_PHASE3_DIR NEW_PHASE3_DIR
"""
import sys
import numpy as np, pandas as pd

OLD, NEW = sys.argv[1], sys.argv[2]
IN_BAND = ["7259_liver_sbr_Male_26-U1", "7260_liver_sbr_Male_26-U1",
           "7001_liver_sham_Male_52-U1", "7248_liver_sham_Male_26-U1",
           "7352_liver_sham_Male_2-U1", "7435_liver_sham_Male_10-U1"]


def load(d, call="tierA_p95"):
    m = pd.read_csv(f"{d}/main_fits.csv")
    m = m[m.section.isin(IN_BAND) & (m.call == call) & (m.stratum == "all")]
    return m


def rep(m):
    return m[(m.beta_naive > 0) & (m.beta_base_lo > 0)]


o, n = load(OLD), load(NEW)
print("fits  old %d  new %d ; reportable old %d  new %d"
      % (len(o), len(n), len(rep(o)), len(rep(n))))

print("\n### per module, PRIMARY call tierA_p95, in-band, reportable fits ###")
hdr = ("%-22s %5s %5s | %7s %7s | %7s %7s | %7s %7s | %7s %7s"
       % ("module", "n_o", "n_n", "b/sd_o", "b/sd_n", "N2_o", "N2_n",
          "N5N6_o", "N5N6_n", "N1_o", "N1_n"))
print(hdr); print("-" * len(hdr))
ro, rn = rep(o), rep(n)
for mod in sorted(set(o.module) | set(n.module)):
    a, b = ro[ro.module == mod], rn[rn.module == mod]
    f = lambda d, c: (d[c].median() if len(d) else np.nan)
    print("%-22s %5d %5d | %7.3f %7.3f | %7.3f %7.3f | %7.3f %7.3f | %7s %7s"
          % (mod, len(a), len(b),
             (a.beta_naive / a.sd_y).median() if len(a) else np.nan,
             (b.beta_naive / b.sd_y).median() if len(b) else np.nan,
             f(a, "sf_n2"), f(b, "sf_n2"), f(a, "sf_n6n5"), f(b, "sf_n6n5"),
             "-", "-"))

print("\n### per receiver cell type ###")
hdr = ("%-22s %5s %5s | %7s %7s | %7s %7s | %7s %7s"
       % ("celltype", "n_o", "n_n", "b/sd_o", "b/sd_n", "N2_o", "N2_n",
          "N5N6_o", "N5N6_n"))
print(hdr); print("-" * len(hdr))
for ct in sorted(set(o.celltype) | set(n.celltype)):
    a, b = ro[ro.celltype == ct], rn[rn.celltype == ct]
    f = lambda d, c: (d[c].median() if len(d) else np.nan)
    print("%-22s %5d %5d | %7.3f %7.3f | %7.3f %7.3f | %7.3f %7.3f"
          % (ct, len(a), len(b),
             (a.beta_naive / a.sd_y).median() if len(a) else np.nan,
             (b.beta_naive / b.sd_y).median() if len(b) else np.nan,
             f(a, "sf_n2"), f(b, "sf_n2"), f(a, "sf_n6n5"), f(b, "sf_n6n5")))

print("\n### the N7 sender axis (all calls present in both) ###")
hdr = ("%-14s %5s %5s | %7s %7s | %7s %7s | %7s %7s"
       % ("call", "rep_o", "rep_n", "b/sd_o", "b/sd_n", "N5_o", "N5_n",
          "N2N5N6_o", "N2N5N6_n"))
print(hdr); print("-" * len(hdr))
mo = pd.read_csv(f"{OLD}/main_fits.csv"); mn = pd.read_csv(f"{NEW}/main_fits.csv")
mo = mo[mo.section.isin(IN_BAND) & (mo.stratum == "all")]
mn = mn[mn.section.isin(IN_BAND) & (mn.stratum == "all")]
for c in sorted(set(mo.call) | set(mn.call)):
    a = rep(mo[mo.call == c]); b = rep(mn[mn.call == c])
    f = lambda d, col: (d[col].median() if len(d) else np.nan)
    print("%-14s %5d %5d | %7.3f %7.3f | %7.3f %7.3f | %7.3f %7.3f"
          % (c, len(a), len(b),
             (a.beta_naive / a.sd_y).median() if len(a) else np.nan,
             (b.beta_naive / b.sd_y).median() if len(b) else np.nan,
             f(a, "sf_n5"), f(b, "sf_n5"),
             f(a, "sf_n2n5n6"), f(b, "sf_n2n5n6")))
