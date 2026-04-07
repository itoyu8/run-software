#!/bin/bash
# Usage: bash util/python-container/run_python.sh script.py [args...]
# Usage: bash util/python-container/run_python.sh -c "import pandas; print(pandas.__version__)"
#
# Wrapper to run Python 3 via Singularity container on HPC.
# Passes all arguments directly to python3.

set -euo pipefail

CONTAINER="/home/itoyu8/singularity/python3_0.1.0.sif"

if [ ! -f "${CONTAINER}" ]; then
    echo "Error: Container not found: ${CONTAINER}" >&2
    exit 1
fi

singularity exec \
    --bind /home/itoyu8/:/home/itoyu8/,/lustre1/:/lustre1/ \
    "${CONTAINER}" \
    python3 "$@"
