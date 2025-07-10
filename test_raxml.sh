#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
INPUT_DIR="${DATA_DIR}/guidance_outputs/guidance_cds_output"
OUTDIR="${DATA_DIR}/raxml_cds_output"
CLEANED_DIR="${DATA_DIR}/cleaned_alignment"

mkdir -p "$OUTDIR" "$CLEANED_DIR"
cd "$INPUT_DIR" || exit 1

# ----------------------------------------
# Step 1: Clean and deduplicate the MSA
# ----------------------------------------
echo "🔧 Trimming and deduplicating MSA..."
TRIMMED_MSA="${CLEANED_DIR}/MSA.trimmed.fa"
DEDUP_MSA="${CLEANED_DIR}/MSA.trimmed.dedup.fa"

trimal -in MSA.MAFFT.Without_low_SP_Col.With_Names -out "$TRIMMED_MSA" -gt 0.8
seqkit rmdup "$TRIMMED_MSA" -o "$DEDUP_MSA"

# ----------------------------------------
# Step 2: Get best-fit model
# ----------------------------------------
MODEL_FILE="${DATA_DIR}/modeltest_cds_output/modeltest_results.out"
BEST_MODEL=$(awk '/Best model according to BIC/{flag=1; next} flag && /^Model:/ {print $2; exit}' "$MODEL_FILE")
FALLBACK_MODEL="LG+G4"

echo "📦 Best model according to BIC: $BEST_MODEL"

# ----------------------------------------
# Step 3: Run RAxML-NG with best-fit model
# ----------------------------------------
echo "🌲 Running RAxML-NG with best model..."
if ! raxml-ng \
    --msa "$DEDUP_MSA" \
    --model "$BEST_MODEL" \
    --prefix "${OUTDIR}/raxml_output" \
    --threads 4 \
    --seed 12345 \
    --bs-trees 100 \
    --all; then

    echo "❌ RAxML-NG failed with model $BEST_MODEL"
    echo "🔁 Retrying with fallback model: $FALLBACK_MODEL"

    # Try fallback model
    raxml-ng \
        --msa "$DEDUP_MSA" \
        --model "$FALLBACK_MODEL" \
        --prefix "${OUTDIR}/raxml_output_fallback" \
        --threads 4 \
        --seed 12345 \
        --bs-trees 100 \
        --all

    echo "✅ Tree built using fallback model: $FALLBACK_MODEL"
else
    echo "✅ Tree built successfully using best model: $BEST_MODEL"
fi