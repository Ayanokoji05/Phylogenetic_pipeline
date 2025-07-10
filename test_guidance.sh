#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
OUTDIR="${DATA_DIR}/guidance_outputs/guidance_cds_output"

cd "${DATA_DIR}/guidance_input" || exit 1

# Clean headers
sed -E 's/^>([^ ]+).*/>\1/' combined_filtered_cds.fa > cleaned.fa

# Validate
SEQ_COUNT=$(seqkit stats cleaned.fa | awk 'NR==2 {print $4}')
[ "$SEQ_COUNT" -lt 4 ] && { echo "❌ Need ≥4 sequences"; exit 1; }

# Run GUIDANCE
mkdir -p "$OUTDIR"
guidance.pl \
    --program GUIDANCE \
    --seqFile cleaned.fa \
    --msaProgram MAFFT \
    --seqType aa \
    --outDir "$OUTDIR" \
    --bootstraps 10 \
    --proc_num 4

# Validate basic output
if ! [[ -s "$OUTDIR/MSA.MAFFT.aln" ]]; then
    echo "❌ GUIDANCE2 failed or did not produce alignment."
    exit 1
fi