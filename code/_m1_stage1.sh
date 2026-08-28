#!/usr/bin/env bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
#
# THE HEAD OF THE M1 PIPELINE -- the piece that was missing from the tag.
#
# Every other M1 driver WAITS for these three stages: _m1_chain.sh blocks until
# logs/m1_main.log has 99 `[main]` lines and logs/m1_perm{,_c1}.log have 6 each;
# _m1_driver_{a,b}.sh block on pgrep; _m1_driver_b2.sh blocks on a PID.  Nothing
# started them.  On a clean machine pgrep matches nothing and the log files do not
# exist, so every one of those waits falls straight through and stages 2-5 run against
# absent or stale inputs (AUDIT_REPRODUCIBILITY B6).  This script is stage 1, and it
# writes exactly the three log files the rest of the chain waits on.
#
# The invocations are recovered from two independent records that agree:
#   * logs/m1_chain.log (a `set -x` trace of stages 2+), logs/m1_main.log,
#     logs/m1_perm.log, logs/m1_perm_c1.log -- section list, call list, n_jobs, n_perm;
#   * reports/WRITING_PACK.md:378 and reports/CS_PHASE8_M1_RERUN.md:478,91.
#
# SEED ORDER MATTERS (AUDIT_REPRODUCIBILITY B5).  run_phase3_nulls.py seeds each job as
# base + step_i*i + step_j*j where i is the POSITION of the section and j the POSITION of
# the call in the --calls list.  A different call order is a different permutation null,
# silently.  The keyword forms (`--calls all9`, `all`, `tierA_pm`) expand from constants
# in the script and are therefore order-stable; a hand-written comma list is not.  Use the
# keywords, and if you must pass a list, pass the exact string recorded here.
#
# Usage:
#   bash code/_m1_stage1.sh            # preflight + print what would run (default)
#   bash code/_m1_stage1.sh --run      # actually launch stage 1, then hand over to
#                                      # _m1_chain.sh, which runs stages 2-5
#
# Cost, from the original run: --stage main ~17 min at 22 jobs / ~25 min at 66;
# --stage perm and --stage perm_c1 are 1,000 permutations over 6 in-band sections and
# took 10-70 minutes PER (section, call) job.  Budget hours, not minutes, and expect
# tens of GB of RAM: the container cgroup, not `free`, is the ceiling.
set -u
R=/workspace
L=$R/logs
PY="$SASP_PYTHON"
RUN=0
[ "${1:-}" = "--run" ] && RUN=1

fail() { echo "PREFLIGHT FAILED: $*" >&2; exit 1; }

# ---------------------------------------------------------------- preflight
echo "interpreter: $PY  ($("$PY" -V 2>&1))"

# 1. the Phase 3 cache.  run_phase3_nulls.py now refuses when the filter empties the
#    section list, but a PARTIAL cache is still silently accepted -- and --stage main
#    writes main_fits.csv with no tag and no merge, so a partial run REPLACES the
#    nine-call table with a smaller one.  Check for all eleven here.
SECS=$("$PY" - <<'PYEOF'
import sys; sys.path.insert(0, '/workspace/code')
import sasp_phase3 as P
print(' '.join(P.ALL_SECTIONS))   # what `--sections all` resolves to
PYEOF
) || fail "cannot import sasp_phase3"
miss=""
for s in $SECS; do
  [ -f "$R/data/processed/cache3/$s.npz" ] || miss="$miss $s"
done
if [ -n "$miss" ]; then
  fail "no Phase 3 cache for:$miss
  data/processed/cache3/ is gitignored (166 MB, regenerable) so a fresh clone has none.
  Build it before running stage 1 -- see sasp_phase3.prep / code/README.md."
fi
echo "cache3: all $(echo $SECS | wc -w) sections present"

# 2. the features stage 1 needs.  The tag phase8-frozen shipped a run_phase3_nulls.py
#    that had none of these, so it could not produce the results attributed to it
#    (AUDIT_REPRODUCIBILITY B1/B2).  Check rather than assume: a silent regression here
#    is exactly how the last three days were lost.
"$PY" - <<'PYEOF' || exit 1
import sys; sys.path.insert(0, '/workspace/code')
import run_phase3_nulls as RN, sasp_phase3 as P, inspect
bad = []
src = inspect.getsource(RN)
for need in (('perm_c1',), ('"--tag"', "'--tag'"), ('all9',), ('is_permodule',), ('_expand',)):
    if not any(n in src for n in need):
        bad.append(need[0])
try:
    P.Sec.sender_mask
    if 'tierApm_p' not in inspect.getsource(P.Sec.sender_mask):
        bad.append('sender_mask tierApm_p')
except Exception as e:
    bad.append('sender_mask (%s)' % e)
if bad:
    print('PREFLIGHT FAILED: run_phase3_nulls.py / sasp_phase3.py are missing: %s'
          % ', '.join(bad), file=sys.stderr)
    print('This is the phase8-frozen revision, which cannot produce perm_nulls_c1.csv,\n'
          'the _pm variant, or the three tierApm_* calls in main_fits.csv.', file=sys.stderr)
    raise SystemExit(1)
print('run_phase3_nulls.py: perm_c1, --tag, all9, is_permodule, _expand all present')
print('sasp_phase3.py: sender_mask resolves tierApm_p*')
PYEOF

mkdir -p "$L"

MAIN="$PY -u run_phase3_nulls.py --stage main    --sections all    --calls all9      --n-jobs 24"
PERM="$PY -u run_phase3_nulls.py --stage perm    --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6"
PC1="$PY -u run_phase3_nulls.py --stage perm_c1 --sections inband --calls tierA_p95 --n-perm 1000 --n-jobs 6"

echo
echo "stage 1 would run, from $R/code:"
echo "  $MAIN  > $L/m1_main.log"
echo "  $PERM  > $L/m1_perm.log"
echo "  $PC1  > $L/m1_perm_c1.log"
echo "then: bash _m1_chain.sh   (stages 2-5; it waits on those three logs)"

if [ "$RUN" != 1 ]; then
  echo
  echo "dry run -- pass --run to launch.  These overwrite results/phase3/main_fits.csv,"
  echo "perm_nulls.csv, perm_curves.csv and perm_nulls_c1.csv, three of which carry the"
  echo "checksums in the phase8-frozen tag annotation.  Commit or copy them first."
  exit 0
fi

cd "$R/code"
$MAIN > "$L/m1_main.log" 2>&1 &
$PERM > "$L/m1_perm.log" 2>&1 &
$PC1  > "$L/m1_perm_c1.log" 2>&1 &
echo "stage 1 launched: main=$! (and two perm stages); logs in $L"
wait
echo "M1_STAGE1_DONE $(date -u)"
bash "$R/code/_m1_chain.sh"
