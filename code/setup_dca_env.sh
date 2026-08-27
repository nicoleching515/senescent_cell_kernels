#!/bin/bash
# Phase 8, task 8.5 / C7-D2, §6 path 1 -- build the ISOLATED environment that runs DCA.
#
# THIS SCRIPT MUST NEVER BE RUN AGAINST THE MAIN ENVIRONMENT.  It builds a self-contained
# Python 3.8 interpreter and venv under $DCA_ENV_ROOT and installs TensorFlow there and
# only there.  The pinned Python 3.11 stack in requirements.txt is not touched, does not
# gain a TensorFlow, and never imports one: code/_shims_dca_bridge talks to this venv over
# a subprocess with scipy .npz / numpy .npy files.
#
# Why a downloaded interpreter.  DCA 0.3.4 requires tensorflow>=2.0,<2.5 and
# keras>=2.4,<2.6.  Every TensorFlow release below 2.5 ships wheels for cp36/cp37/cp38
# only (verified against the PyPI JSON API: tensorflow 2.4.4 has exactly nine files, all
# cp36/cp37/cp38).  This box has Python 3.10 and 3.11; Ubuntu 22.04 has no python3.8
# package, and /usr/bin/python3.10 has no ensurepip.  A python-build-standalone CPython
# 3.8.19 is therefore downloaded and checksum-verified.
#
# One deviation from DCA's own dependency list, forced: DCA imports dca.train ->
# dca.hyper -> kopt unconditionally, and kopt 0.1.0 calls yaml.load() without a Loader,
# which PyYAML 6 removed.  PyYAML is pinned to 5.4.1 in this venv.  It affects only kopt's
# config parsing, not the autoencoder.
set -euo pipefail
DCA_ENV_ROOT="${DCA_ENV_ROOT:?set DCA_ENV_ROOT to a scratch directory OUTSIDE the repo}"
PY_URL='https://github.com/astral-sh/python-build-standalone/releases/download/20240726/cpython-3.8.19%2B20240726-x86_64-unknown-linux-gnu-install_only.tar.gz'
PY_SHA256='e81ea4dd16e6057c8121bdbcb7b64e2956068ca019f244c814bc3ad907cb2765'

mkdir -p "$DCA_ENV_ROOT"; cd "$DCA_ENV_ROOT"
curl -sSL -o py38.tar.gz "$PY_URL"
echo "$PY_SHA256  py38.tar.gz" | sha256sum -c -
tar xzf py38.tar.gz
./python/bin/python3.8 -m venv v38
./v38/bin/pip install --upgrade "pip<25" "setuptools<70" wheel
./v38/bin/pip install "dca==0.3.4"      # pulls tensorflow 2.4.4 + keras 2.4.3
./v38/bin/pip install "pyyaml==5.4.1"   # see the kopt note above
./v38/bin/python -c "import warnings; warnings.filterwarnings('ignore'); \
import dca.api, tensorflow as tf, keras; print('DCA OK under TF', tf.__version__, 'keras', keras.__version__)"
echo
echo "Interpreter for DCA_VENV_PYTHON:  $DCA_ENV_ROOT/v38/bin/python"
