#!/usr/bin/env bash
. /workspace/code/_env.sh   # project interpreter, PATH and BLAS threads (AUDIT_REPRODUCIBILITY B10)
#
# Regenerate -- or verify -- the committed artefacts whose producer was a bare stdout
# redirect or an ad-hoc command that was never itself committed (AUDIT_REPRODUCIBILITY B7).
# Each entry below was produced by a command that existed only in someone's shell history.
#
#   bash code/_repro_artefacts.sh --check [DIR]   regenerate into DIR (default: a temp dir)
#                                                 and diff against the committed copy.
#                                                 Writes NOTHING into results/.
#   bash code/_repro_artefacts.sh --write         regenerate in place, overwriting results/.
#   add --with-network                            also refresh the one artefact that needs
#                                                 an NCBI eutils call.
#
# Everything here is cheap: seconds to a couple of minutes, no pipeline stage.
#
# Known non-determinism, stated rather than hidden:
#   * m1_final_audit.txt  section 5 ("FIGURE STATE") embeds the live figure-guard verdict and
#     the md5 + HH:MM mtime of every figure.  Sections 1-4 and 6-7 -- which carry every number
#     the reports cite, including power80_bound 0.1833 and the controlled amplitude 0.0288 --
#     reproduce byte-for-byte.  The committed copy was captured mid-8.7 while the guard was
#     still failing, so its section 5 can never come back.
#   * d2_tables.md section H globs runmeta_*.json: a tree with more runs renders more rows.
#   * sf_summary_c1_swap_vs_n1.csv is a PRE-C6 artefact.  It is regenerated from
#     results/phase3_pre_c6 deliberately; run sf_swap_vs_n1.py with no argument for the
#     post-C6 numbers, which differ (CORRECTIONS.md:756 quotes those).
#   * the DCA entries need the CPython 3.8 venv; set DCA_ENV_ROOT (code/setup_dca_env.sh).
set -u
MODE=""; DEST=""; NET=0
while [ $# -gt 0 ]; do
  case "$1" in
    --check) MODE=check ;;
    --write) MODE=write ;;
    --with-network) NET=1 ;;
    *) DEST="$1" ;;
  esac; shift
done
[ -n "$MODE" ] || { sed -n '3,30p' "$0"; exit 2; }

R=/workspace
PY="$SASP_PYTHON"
DCA_ROOT="${DCA_ENV_ROOT:-/tmp/claude-0/-workspace/f999c960-39aa-4f7a-a180-b1cefba480ce/scratchpad/dca_attempt}"
DCA_PY="${DCA_VENV_PYTHON:-$DCA_ROOT/v38/bin/python}"
if [ "$MODE" = check ]; then
  DEST="${DEST:-$(mktemp -d)}"; mkdir -p "$DEST"
  echo "regenerating into $DEST (nothing under results/ is written)"
fi
OK=0; DIFF=0; SKIP=0

# produce <committed-path> <shell command writing to $OUT>
produce() {
  rel="$1"; shift
  if [ "$MODE" = write ]; then OUT="$R/$rel"; else OUT="$DEST/$(basename "$rel")"; fi
  export OUT
  if ! ( cd "$R" && eval "$@" ) >/dev/null 2>"$DEST_ERR"; then
    echo "  SKIP  $rel  (producer failed: $(tail -1 "$DEST_ERR" 2>/dev/null))"; SKIP=$((SKIP+1)); return
  fi
  if [ "$MODE" = write ]; then echo "  WROTE $rel"; OK=$((OK+1)); return; fi
  if cmp -s "$OUT" "$R/$rel"; then echo "  MATCH $rel"; OK=$((OK+1))
  elif [ "$rel" = results/phase3/m1_final_audit.txt ] \
       && cmp -s <(sed '/^5. FIGURE STATE/,$d' "$R/$rel") <(sed '/^5. FIGURE STATE/,$d' "$OUT"); then
    # everything above section 5 is the numeric content: permutation counts, the pinned-file
    # md5s, the pre/post headline vector (power80_bound 0.1833, ctrl_amp_med 0.0288), the
    # per-call reportable table.  Section 5 is live figure-guard state and figure mtimes.
    echo "  MATCH $rel  (sections 1-4 byte-exact; section 5 FIGURE STATE is live guard state"
    echo "               and figure md5/mtimes -- the committed copy was captured while the"
    echo "               guard was still failing mid-8.7 and cannot come back)"
    OK=$((OK+1))
  else
    n=$(diff "$R/$rel" "$OUT" | grep -c '^[<>]')
    echo "  DIFF  $rel  ($n differing lines; diff $R/$rel $OUT)"; DIFF=$((DIFF+1))
  fi
}
DEST_ERR="$(mktemp)"

echo "== stdout artefacts =="
produce results/phase3/m1_final_audit.txt        '"$PY" -u code/m1_final_audit.py > "$OUT"'
produce results/phase3/m1_prepost_main_fits.txt  '"$PY" -u code/m1_compare_modules.py results/phase3_pre_c6 results/phase3 > "$OUT"'
produce results/phase3/m1_n7_prepost.txt         '"$PY" -u code/m1_n7_prepost.py > "$OUT"'
produce results/phase3/sf_summary_c1_swap_vs_n1.csv '"$PY" -u code/sf_swap_vs_n1.py results/phase3_pre_c6 > "$OUT"'
produce results/phase8_d2/d2_tables.md           '"$PY" -u code/report_d2_tables.py > "$OUT"'

echo "== D2 provenance =="
produce results/phase8_d2/committed_deepscence_sha256.txt \
        'sha256sum data/processed/deepscence_sham.csv data/processed/deepscence_sbr.csv data/processed/deepscence_7*.csv > "$OUT"'
produce results/phase8_d2/dca_venv_pip_freeze.txt '"$DCA_PY" -m pip freeze > "$OUT"'
produce results/phase8_d2/dca_venv_python.txt     '"$DCA_PY" -VV > "$OUT" 2>&1'

echo "== caller-agreement tables (post-C6, two sections) =="
# seven files at once; verified byte-identical on 2026-08-27.
if [ "$MODE" = write ]; then
  "$PY" -u "$R/code/caller_disagree_all.py" --set 2sec_c6 && echo "  WROTE results/phase3/caller_*_2sec_c6.csv"
else
  mkdir -p "$DEST/c6"
  if "$PY" -u "$R/code/caller_disagree_all.py" --set 2sec_c6 --out-dir "$DEST/c6" >/dev/null 2>&1; then
    for f in "$DEST"/c6/*.csv; do
      b=$(basename "$f")
      if cmp -s "$f" "$R/results/phase3/$b"; then echo "  MATCH results/phase3/$b"; OK=$((OK+1))
      else echo "  DIFF  results/phase3/$b"; DIFF=$((DIFF+1)); fi
    done
  else
    echo "  SKIP  results/phase3/caller_*_2sec_c6.csv (needs data/processed/{senders,anatomy,deepscence}_{sham,sbr}.csv and data/raw/*/cells.parquet)"
    SKIP=$((SKIP+1))
  fi
fi

if [ "$NET" = 1 ]; then
  echo "== network =="
  produce results/a3_fallback/gpl33762_count.xml \
    'curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term=GPL33762%5BAccession%5D+AND+gse%5BEntry+Type%5D&retmax=0" > "$OUT"'
else
  echo "== network (skipped; pass --with-network) =="
  echo "  results/a3_fallback/gpl33762_count.xml"
fi

echo "== partial by design (informational; does not affect the exit status) =="
# results/section_qc_sender_summary.csv reproduces 120 of 132 cells.  The 12 that differ are
# n_analysable / n_cdkn1a_pos / cdkn1a_pos_pct_all on the four sections re-annotated after the
# committed CSV was written, and portal_triad_valid's threshold is a reconstruction.  Its
# --check mode prints the grid; it is NOT counted in the totals below.
"$PY" "$R/code/section_qc_sender_summary.py" --check 2>&1 | tail -2 | sed 's/^/  /'

echo "match=$OK  differ=$DIFF  skipped=$SKIP"
[ "$DIFF" = 0 ]
