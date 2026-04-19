#!/usr/bin/env bash
# install.sh — creates the 'remap' conda environment with all pipeline dependencies.
# Run once before using run_pipeline.sh.
#
# Usage:
#   bash install.sh
#
# Requires: conda or mamba (mamba recommended for speed).
set -euo pipefail

ENV_NAME="remap"

# Prefer mamba if available
if command -v mamba &>/dev/null; then
    CONDA_CMD="mamba"
else
    CONDA_CMD="conda"
    echo "[install] mamba not found, falling back to conda (slower)"
fi

echo "[install] Creating environment '${ENV_NAME}'..."
${CONDA_CMD} create -y -n "${ENV_NAME}" -c bioconda -c conda-forge \
    minimap2 pysam samtools bcftools pandas \
    ucsc-liftover pip

echo "[install] Installing CrossMap via pip (bioconda's crossmap has Python-version pins that conflict with current conda-forge Python)..."
conda run -n "${ENV_NAME}" pip install --quiet CrossMap

echo ""
echo "[install] Done. Activate the environment before running the pipeline:"
echo "    conda activate ${ENV_NAME}"
echo "    bash run_pipeline.sh -i <manifest.csv> -r <reference.fa> -a <assembly_name> -o <output_dir>"
