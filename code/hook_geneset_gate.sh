#!/usr/bin/env bash
# Claude Code PostToolUse hook -- Phase 8 hazard B6.
#
# Fires the section-11 disjointness gate as soon as any agent writes a gene set
# or a panel file, so a change cannot sit on disk unchecked until freeze time.
# `B_oxidative_stress` clears the >=30 floor by one gene on the mouse arm, so a
# gene-set edit is never a safe no-op. Reads the hook payload (JSON) on stdin
# and exits 0 immediately unless the edited path is a watched one.
# Exit 2 = blocking feedback to the model.
set -u
ROOT=/workspace
PY="${SASP_PYTHON:-}"
if [ -z "$PY" ]; then
  # Prefer the interpreter that carries the scientific stack. The gate itself is
  # stdlib-only, but preferring a full interpreter means it keeps working if the
  # gate ever grows a numpy dependency. envs/sasp311 is the persistent fallback.
  if python3 -c "import sys" >/dev/null 2>&1; then PY="python3";
  elif [ -x "$ROOT/envs/sasp311/bin/python" ]; then PY="$ROOT/envs/sasp311/bin/python";
  else PY="python"; fi
fi

FILE="$("$PY" -c '
import json,sys
try: d = json.load(sys.stdin)
except Exception: sys.exit(0)
ti = d.get("tool_input") or {}
tr = d.get("tool_response") or {}
p = ti.get("file_path") or (tr.get("filePath") if isinstance(tr, dict) else None) or ""
print(p)
' 2>/dev/null)"
[ -n "$FILE" ] || exit 0

# Two watched families, two different gates.
#
#   reports/, results/, README.md  -> the section-10 prohibition checker.  These paths used to
#     `exit 0` here unconditionally, which meant the only always-on hook deliberately skipped
#     every path where a section-10 violation would actually be written
#     (reports/AUDIT_PREREG_VS_CODE.md PART 4, "Coverage hole in the one always-on hook").
#   genesets/, panel files        -> the section-11 disjointness gate, as before.
case "$FILE" in
  */reports/*.md|*/reports/*.txt|*/results/*.md|*/results/*.txt|*/README.md)
    OUT="$("$PY" "$ROOT/code/check_prohibitions.py" "$FILE" 2>&1)"
    if [ $? -ne 0 ]; then
      printf 'SECTION 10 PROHIBITION VIOLATION in %s\n\n%s\n' "$FILE" "$OUT" >&2
      exit 2
    fi
    exit 0 ;;
  */results/*|*/reports/*) exit 0 ;;
  */genesets/*|*/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv|*/GSE310392_Q6VTXC_mMulti_100g_gene_list.csv.gz) ;;
  *) exit 0 ;;
esac

OUT="$("$PY" "$ROOT/code/gate_genesets_guard.py" --quiet 2>&1)"
RC=$?
if [ $RC -ne 0 ]; then
  printf 'SECTION 11 DISJOINTNESS GATE FAILED after writing %s\n\n%s\n' "$FILE" "$OUT" >&2
  exit 2
fi
printf 'section-11 gate re-run after %s: PASS\n' "$FILE"
exit 0
