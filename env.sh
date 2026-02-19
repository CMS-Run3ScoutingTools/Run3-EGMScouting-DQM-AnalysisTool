#!/usr/bin/env bash

# This file is meant to be sourced:
#   source env.sh
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
  echo "Please source this file instead of executing it:"
  echo "  source env.sh"
  exit 1
fi

export DQM_PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# LCG runtime (ROOT, PyROOT, scientific stack)
LCG_SETUP="/cvmfs/sft.cern.ch/lcg/views/LCG_105b/x86_64-el9-gcc13-opt/setup.sh"
if [[ ! -f "${LCG_SETUP}" ]]; then
  echo "[env][ERROR] LCG setup not found: ${LCG_SETUP}"
  return 1
fi
source "${LCG_SETUP}"

# Allow 'python -m' imports from scripts/ package layout.
if [[ -n "${PYTHONPATH:-}" ]]; then
  export PYTHONPATH="${DQM_PIPELINE_ROOT}/scripts:${PYTHONPATH}"
else
  export PYTHONPATH="${DQM_PIPELINE_ROOT}/scripts"
fi

export DQM_PIPELINE_ENV_LOADED=1

echo "[env] loaded"
echo "[env] DQM_PIPELINE_ROOT=${DQM_PIPELINE_ROOT}"
echo "[env] PYTHONPATH=${PYTHONPATH}"
