#!/usr/bin/env bash
# The project interpreter, in one place.  Source this from every driver:
#
#     . "$(dirname "${BASH_SOURCE[0]}")/_env.sh"
#
# WHY THIS FILE EXISTS
# --------------------
# Every tracked driver used to call bare `python3`.  On this container that resolves to
# /usr/bin/python3 (3.11.10) whose sys.path is /usr/local/lib/python3.11/dist-packages -- a
# SECOND, complete, 249-package scientific stack that no manifest in this repository
# describes, and which has been wiped twice.  The environment `requirements.txt` actually
# describes is a 144-package venv at /workspace/envs/sasp311, and until now that path was
# named nowhere in the repo.  So a reader who followed README.md ("pip install -r
# requirements.txt" into a venv) and then ran any driver silently got the OTHER stack --
# working today only by coincidence, since all 21 pinned versions happen to match in both.
# (reports/AUDIT_REPRODUCIBILITY.md B10, D2.)
#
# WHAT IT DOES
#   * SASP_PYTHON   -- absolute path to the project interpreter.  Override to use another.
#   * PATH          -- the interpreter's bin/ is prepended, so bare `python3` inside a driver
#                      resolves to the project interpreter too.  Drivers therefore need no
#                      other edit, and a stray `python3` cannot silently reach the overlay.
#   * fails LOUDLY (exit 1) if the interpreter is missing -- never falls back to /usr/bin.
#
# Rebuild the environment with:
#     python3.11 -m venv /workspace/envs/sasp311
#     /workspace/envs/sasp311/bin/pip install -r /workspace/requirements.txt
# The DCA arm is a SEPARATE CPython 3.8 venv; see code/setup_dca_env.sh and DCA_ENV_ROOT.
set -u

SASP_ROOT="${SASP_ROOT:-/workspace}"
SASP_PYTHON="${SASP_PYTHON:-$SASP_ROOT/envs/sasp311/bin/python}"

if [ ! -x "$SASP_PYTHON" ]; then
  cat >&2 <<MSG
FATAL: the project interpreter is missing.

  expected: $SASP_PYTHON

This project does NOT run on the system python3: /usr/bin/python3 loads
/usr/local/lib/python3.11/dist-packages, a different stack that no manifest here
describes.  Build the project environment first:

  python3.11 -m venv $SASP_ROOT/envs/sasp311
  $SASP_ROOT/envs/sasp311/bin/pip install -r $SASP_ROOT/requirements.txt

or point SASP_PYTHON at an interpreter that satisfies requirements.txt.
MSG
  exit 1
fi

PATH="$(cd "$(dirname "$SASP_PYTHON")" && pwd):$PATH"
export SASP_PYTHON SASP_ROOT PATH
# Keep BLAS single-threaded unless a driver says otherwise: every parallel stage in this
# project fans out with joblib and oversubscribes badly without this.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}" \
       MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}" \
       OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
