# Commit plan — for the PI on return

**Written 2026-08-27 during the autonomous window. Nothing here has been executed:
no commits, no pushes, no tags were made, per standing instruction.**

You return to **590 changed paths** (107 modified, 483 untracked once directories
are expanded). This is the order to deal with them and the one hazard to clear first.

---

## 0. CLEAR THE HAZARD FIRST — before any `git checkout`

**`git checkout -- code/` is currently destructive.**

`git show HEAD:code/build_genesets.py` contains **none** of the guards added today.
Restoring `code/` from HEAD reinstates a script whose `SCRATCH` constant points at
a dead per-session `/tmp` path — it then globs **zero** MSigDB JSONs and
**silently overwrites `genesets/*.txt` with EMPTY Tier B modules, exiting 0**.

Two working-tree guards prevent that, and neither exists at HEAD:
- it now refuses to run when the archived MSigDB pin is missing
- it now refuses to revert the promoted C6 sets unless `ALLOW_OVERWRITE_C6=1`

**Commit `code/` early.** Until you do, the repo can destroy its own gene sets with
one ordinary command.

Also at risk: **`results/phase3_pre_c6/` is 98 files, 0 tracked** — the sole copy of
the baseline that `reports/CORRECTIONS.md` compares every number against.

---

## 1. What is pending, by category

| Category | Modified | Untracked | Note |
|---|---|---|---|
| `results/` (evidence) | 67 | 255 | Backs C1, A7, D2, compmatch, Moran, the caller gate |
| `genesets/` (**frozen sets**) | 3 | 104 | `human/` and `mouse_c6/` are **0 tracked** |
| `code/` (producers + guards) | 11 | 72 | **Contains the hazard fix — commit early** |
| `figures/` | 16 | 25 | The 16 are 8.7's sanctioned regeneration pass |
| `reports/` | 5 | 23 | The analysis record |
| other | 5 | 4 | incl. the Phase 7 plan's superseded banner |

## 2. Suggested commit order

Small, reviewable commits beat one large one — the pre-registration should be able
to cite specific hashes.

1. **`code/`** — producers and the four silent-failure guards. **Clears §0's hazard.**
2. **`genesets/`** — the frozen Tier A–E sets, both arms, plus `mouse_c6/` and the
   MSigDB/CellMarker pins. This is what `phase8-frozen` must point at.
3. **`results/`** — the evidence base. **Include `results/phase3_pre_c6/`**; without
   it the corrections ledger references a baseline that exists nowhere.
4. **`figures/`** — the 8.7 regeneration pass plus the new Phase 8 figures.
   `python3 code/check_figures_guard.py` passes at **52/52**; re-run after staging.
5. **`reports/`** — including `PREREG_PHASE8.md` and `CORRECTIONS.md`.
6. **THEN tag `phase8-frozen`**, and fill the 7 `TBD` hash fields in
   `PREREG_PHASE8.md` §1 — they are placeholders precisely because the tag does not
   exist yet.

## 3. Verify before tagging

```bash
python3 code/check_figures_guard.py          # expect: OK, all 52 committed figures match
python3 code/gate_genesets_guard.py          # expect: rc=0, both arms pass
md5sum results/phase3/{perm_nulls.csv,sf_summary.csv,summary_phase3.txt}
#   d906394958dbe1b99981756290c511fa  perm_nulls.csv
#   a5ccc9b0e81f4c335e8039e975ec1975  sf_summary.csv
#   dc92ddc6605eef52f6359aeab4e16fd7  summary_phase3.txt
```

`pre-c6-genesets` already exists and captures the pre-C6 state; `git show
pre-c6-genesets:<path>` recovers any pre-C6 baseline.

## 4. What is deliberately NOT to be committed

Already in `.gitignore`, verified:
- `envs/` — 6.6 GB isolated venv on the network volume
- `data/raw_h1/` — 525 MB, re-fetchable via `code/fetch_h1_geo.sh`
- `figures/revised_candidates/` — superseded figure candidates, kept for comparison
- `__pycache__/` — 2,425 directories were removed today

## 5. Two things that do not wait for the commit

1. **`reports/SUBMISSION_PATCH_2026-08-29.md`** — the manuscript is **not in this
   repo**, so this must be applied **by hand**. Deadline is in two days. Its §1–§3
   were re-derived against settled numbers; §4a and §9 carry corrections that
   reverse what earlier drafts said.
2. **The venue question** — evidence is in `SASP_Kernel_Master_Plan.md` §29's
   flagged subsection. The Primary/Secondary paragraphs were left **byte-unchanged**
   because that call is yours.

---

**Read first:** `reports/COMPLETED_TASKS.md` (81 entries, evidence file per row,
including every correction made to the coordinator's own statements), then
`reports/CORRECTIONS.md` (what moved and why).

**The headline did not move:** controlled amplitude **0.0288** against an 80 %-power
detectable bound of **0.1833**; SF under N2+N5+N6 **0.0885**. §18 outcome **A**
stands, and the bound *tightened* from 0.203.
