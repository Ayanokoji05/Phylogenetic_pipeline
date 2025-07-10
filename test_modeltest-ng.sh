#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
INPUT_DIR="${DATA_DIR}/guidance_outputs/guidance_cds_output"
OUTDIR="${DATA_DIR}/modeltest_cds_output"

mkdir -p "$OUTDIR"
cd "$INPUT_DIR" || exit 1

modeltest-ng \
    -i MSA.MAFFT.Without_low_SP_Col.With_Names \
    -d aa \
    -o "${OUTDIR}/modeltest_results" \
    -T raxml \
    -p 4