#!/usr/bin/env python3
"""
CoreScence circularity, MOUSE arm -- derived from files, never typed in.

Why this file exists
--------------------
`reports/PREREG_PHASE8_genesets.md` sec 7b promises that the gene-set figures are built
"from the CSVs above -- nothing is recomputed or typed in".  The mouse CoreScence anchor
was the one exception: `24/35 = 69%` was a string literal in
`code/make_figure_genesets.py` and in `code/gate_disjointness_human.py`, produced by no
script.  The 2026-08-27 fact-check (`reports/AUDIT_PHASE8_FACTCHECK.md`, item M1) flagged
it.  Re-derivation shows the numerator 24 is right for the STRICT MGI convention but the
denominator 35 is reproducible under no convention at all.  The project's own committed
evidence says 33:

    results/phase3/n8_disjointness_*.csv  ->  corescence_on_panel = 33   (committed at HEAD)
    logs/ds_smoke.log                     ->  "on our ortholog-mapped panel: 31" (strict MGI)
    logs/caller2.log                      ->  "17 up / 14 down on mouse panel"  (= 31, strict)

so the published `24/35 = 69%` understated the mouse arm's circularity by ~10 points.
This script computes the number instead.

Two mapping conventions are reported, both from disk:

  strict    CoreScence human symbol must be the image of a mouse panel gene under the
            pinned 1:1 MGI map (`genesets/mouse_human_orthologs_MGI.csv`).  -> 31 on panel.
  fallback  as strict, plus the project's documented Title-case fallback for symbols the
            pinned map simply lacks a row for (`CDKN2B`, `CXCL1`; the same map-gap class
            `genesets/human/_symbol_resolutions.csv` records and that the Tier C asymmetry
            table already corrects for).  -> 33 on panel.

**`fallback` is the convention to cite.**  It is what `code/run_phase3_n8.py::corescence_mouse`
implements, i.e. what the mouse arm's published Phase 3 circularity numbers were actually
computed under, and it matches the committed `corescence_on_panel = 33`.

Products: `results/phase7_jobA/corescence_circularity_mouse.json`, plus a printed table.
Run: python3 /workspace/code/corescence_circularity.py
"""
from __future__ import annotations

import csv
import gzip
import json
import os
import subprocess

W = "/workspace"
RES = W + "/results/phase7_jobA"
CORE_GS = "/usr/local/lib/python3.11/dist-packages/DeepScence/data/coreGS_v2.csv"
ORTHO_CSV = W + "/genesets/mouse_human_orthologs_MGI.csv"
MODS = ["tnfa_nfkb_proximal", "il6_jak_stat3", "interferon_response", "downstream_arrest",
        "emt_ecm", "oxidative_stress", "secondary_senescence"]
PRE_C6_TAG = "pre-c6-genesets"


def mouse_panel():
    """The 5,097-gene AUTHORITATIVE mouse panel (deviation D17): both panel files,
    minus the 9 genotyping probes.  Verified elsewhere to equal the .h5 feature list."""
    mp = {r["gene_name"] for r in
          csv.DictReader(open(W + "/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv"))}
    p100 = [r["Gene"] for r in
            csv.DictReader(gzip.open(W + "/GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz", "rt"))]
    geno = {g for g in p100 if ("_WT" in g or "_ALT" in g or "_del_" in g or "_splice_" in g)}
    panel = (mp | set(p100)) - geno
    assert len(panel) == 5097, len(panel)
    return panel


def corescence(occ_min=5):
    """CoreScence v2 at DeepScence's own occurrence threshold.  HUMAN symbols."""
    rows = list(csv.DictReader(open(CORE_GS)))
    return sorted({r["gene_symbol"] for r in rows
                   if r["occurrence"] and float(r["occurrence"]) >= occ_min})


def _read_pre_c6(name):
    r = subprocess.run(["git", "-C", W, "show", "%s:genesets/%s.txt" % (PRE_C6_TAG, name)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit("cannot read pre-C6 %s from tag %s: %s"
                         % (name, PRE_C6_TAG, r.stderr.strip()))
    return {l.strip() for l in r.stdout.splitlines() if l.strip()}


def _read_now(name):
    return {l.strip() for l in open("%s/genesets/%s.txt" % (W, name)) if l.strip()}


def derive():
    panel = mouse_panel()
    core = corescence()
    ortho = {r["mouse_symbol"]: r["human_symbol"] for r in csv.DictReader(open(ORTHO_CSV))}

    # human symbol -> mouse panel gene(s), the direction the DeepScence mouse run uses
    # (`code/run_deepscence.py` renames mouse symbols to their 1:1 human orthologs).
    h2m = {}
    for g in panel:
        if g in ortho:
            h2m.setdefault(ortho[g], set()).add(g)

    strict = sorted(c for c in core if c in h2m)
    fallback_extra = sorted(c for c in core if c not in h2m and c.capitalize() in panel)
    fallback = sorted(set(strict) | set(fallback_extra))
    off_panel = sorted(set(core) - set(fallback))

    def mouse_forms(c):
        s = set(h2m.get(c, ()))
        if not s and c.capitalize() in panel:
            s = {c.capitalize()}
        return s

    out = dict(source=CORE_GS, n_occ_ge5=len(core), mouse_panel=len(panel),
               convention_to_cite="fallback",
               n_on_panel_strict=len(strict), n_on_panel_fallback=len(fallback),
               fallback_extra=fallback_extra, off_mouse_panel=off_panel,
               configurations={})

    for tag, rd in (("pre_C6", _read_pre_c6), ("C6_promoted", _read_now)):
        B = {m: rd("B_" + m) & panel for m in MODS}
        U = set().union(*B.values())
        cfg = {}
        for conv, onp in (("strict", strict), ("fallback", fallback)):
            in_any = sorted(c for c in onp if mouse_forms(c) & U)
            cfg[conv] = dict(n_in_any_B=len(in_any), n_on_panel=len(onp),
                             frac=round(len(in_any) / len(onp), 4),
                             genes_in_any_B=in_any,
                             per_module={m: sorted(c for c in onp if mouse_forms(c) & B[m])
                                         for m in MODS})
        out["configurations"][tag] = cfg

    # self-check against the project's own committed Phase 3 evidence
    out["crosscheck"] = dict(
        n8_disjointness_corescence_on_panel=33,
        n8_disjointness_per_module_preC6={
            "downstream_arrest": 10, "emt_ecm": 9, "il6_jak_stat3": 5,
            "interferon_response": 5, "oxidative_stress": 0,
            "secondary_senescence": 14, "tnfa_nfkb_proximal": 8},
        source="results/phase3/n8_disjointness_*.csv @ HEAD (committed)")
    assert out["n_on_panel_fallback"] == 33, out["n_on_panel_fallback"]
    pm = {m: len(v) for m, v in out["configurations"]["pre_C6"]["fallback"]["per_module"].items()}
    assert pm == out["crosscheck"]["n8_disjointness_per_module_preC6"], pm
    return out


def reference_string(d=None, config="pre_C6", conv=None):
    """The one-line mouse-arm reference other scripts quote, e.g. '26/33 = 79%'."""
    d = d if d is not None else derive()
    conv = conv or d["convention_to_cite"]
    c = d["configurations"][config][conv]
    return "%d/%d = %.0f%%" % (c["n_in_any_B"], c["n_on_panel"], 100 * c["frac"])


def load(path=None):
    """Read the derived JSON, or derive it if it has not been written yet."""
    path = path or (RES + "/corescence_circularity_mouse.json")
    if os.path.exists(path):
        return json.load(open(path))
    return derive()


if __name__ == "__main__":
    d = derive()
    os.makedirs(RES, exist_ok=True)
    print("CoreScence v2 occurrence >= 5 : %d human symbols" % d["n_occ_ge5"])
    print("mouse panel (D17 authoritative): %d" % d["mouse_panel"])
    print("on the ortholog-mapped mouse panel: strict %d, with the documented Title-case "
          "fallback %d (adds %s)"
          % (d["n_on_panel_strict"], d["n_on_panel_fallback"], " ".join(d["fallback_extra"])))
    print("not on the mouse panel under either convention: %s" % " ".join(d["off_mouse_panel"]))
    print()
    print("%-14s %-9s %-16s" % ("configuration", "mapping", "in >=1 Tier B module"))
    for tag in ("pre_C6", "C6_promoted"):
        for conv in ("strict", "fallback"):
            c = d["configurations"][tag][conv]
            print("%-14s %-9s %d/%d = %.1f%%"
                  % (tag, conv, c["n_in_any_B"], c["n_on_panel"], 100 * c["frac"]))
    print()
    print("CITE (%s): pre-C6 %s ; C6-promoted %s"
          % (d["convention_to_cite"], reference_string(d, "pre_C6"),
             reference_string(d, "C6_promoted")))
    print("per module, C6-promoted:", {m: len(v) for m, v in
                                       d["configurations"]["C6_promoted"]["fallback"]
                                       ["per_module"].items()})
    json.dump(d, open(RES + "/corescence_circularity_mouse.json", "w"), indent=1)
    print("\nWritten to %s/corescence_circularity_mouse.json" % RES)
