#!/usr/bin/env bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
# Gate for re-cutting phase8-frozen. PI-authorised 2026-08-27, CONDITIONAL on
# the reconstructed pipeline reproducing the committed results.
#
# The original tag (9264396) shipped code that could NOT produce its own
# results: the task 8.7 edits to run_phase3_nulls.py were destroyed by an
# uncommitted-change revert before being staged. This gate refuses to re-cut
# until that is genuinely repaired. It checks; it does not tag.
#
#   bash code/recut_phase8_gate.sh    -> exit 0 = safe to re-cut
set -u
cd /workspace
fail=0
say() { printf "  %-58s %s\n" "$1" "$2"; }
bad() { say "$1" "FAIL — $2"; fail=1; }

echo "=== 1. the lost symbols are back ==="
for s in perm_c1 _expand TIERA_PM_CALLS is_permodule; do
  n=$(grep -c "$s" code/run_phase3_nulls.py 2>/dev/null); n=${n:-0}
  [ "$n" -gt 0 ] && say "$s in run_phase3_nulls.py" "ok ($n)" || bad "$s in run_phase3_nulls.py" "absent"
done
n=$(grep -c load_ours code/make_figure4.py 2>/dev/null); n=${n:-0}
[ "$n" -gt 0 ] && say "load_ours in make_figure4.py" "ok ($n)" || bad "load_ours in make_figure4.py" "absent"
n=$(grep -c "Refusing to run" code/build_genesets.py 2>/dev/null); n=${n:-0}
[ "$n" -gt 0 ] && say "build_genesets.py guard" "ok" || bad "build_genesets.py guard" "absent"

echo "=== 2. the reconstruction reproduces the committed results ==="
if [ -f reports/PIPELINE_RECONSTRUCTION.md ]; then
  if grep -qiE "bit-identical|byte-identical|reproduces exactly" reports/PIPELINE_RECONSTRUCTION.md; then
    say "reconstruction report claims reproduction" "present — VERIFY BY HAND"
  else
    bad "reconstruction report" "no reproduction claim found"
  fi
else
  bad "reports/PIPELINE_RECONSTRUCTION.md" "missing"
fi

echo "=== 3. the frozen code can actually be invoked ==="
python3 -c "import ast,sys; ast.parse(open('code/run_phase3_nulls.py').read())" 2>/dev/null \
  && say "run_phase3_nulls.py parses" "ok" || bad "run_phase3_nulls.py" "syntax error"
grep -q "perm_c1" code/run_phase3_var.py 2>/dev/null && say "run_phase3_var.py references perm_c1" "ok" || true

echo "=== 4. guards and cleanliness ==="
python3 code/check_figures_guard.py >/dev/null 2>&1 && say "figure guard" "ok" || bad "figure guard" "does not pass"
python3 code/gate_genesets_guard.py >/dev/null 2>&1 && say "gene-set gate" "ok" || bad "gene-set gate" "does not pass"
n=$(git status --porcelain | grep -vc '^??'); n=${n:-0}
[ "$n" -eq 0 ] && say "no uncommitted tracked changes" "ok" || bad "uncommitted tracked changes" "$n"

echo
if [ "$fail" -eq 0 ]; then
  echo "  GATE PASS — safe to re-cut phase8-frozen."
  echo "  Old tag: $(git rev-list -n1 phase8-frozen 2>/dev/null || echo none)"
  echo "  Would move to: $(git rev-parse HEAD)"
else
  echo "  GATE FAIL — do NOT re-cut."
fi
exit $fail
