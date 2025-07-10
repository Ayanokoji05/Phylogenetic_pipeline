#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data"
CDS_DIR="${DATA_DIR}/cds_filtered"
FILTERED_DIR="${DATA_DIR}/guidance_input"

# Create output directory if not exists
mkdir -p "$FILTERED_DIR"

# Check for CDS input
if [ ! -d "$CDS_DIR" ] || [ -z "$(ls -A "$CDS_DIR")" ]; then
    echo "❌ Error: No CDS files found in $CDS_DIR" >&2
    exit 1
fi

cd "$CDS_DIR" || exit 1

# Clean filenames with spaces
for f in *\ *; do
    [ -e "$f" ] || continue
    mv -- "$f" "${f// /_}"
done

# Initialize combined filtered output
> "${FILTERED_DIR}/combined_filtered_cds.fa"

# Process files ending with _cds.fa
for file in *_cds.fa; do
    [ -e "$file" ] || continue  # Skip if no matching files
    base="${file%_cds.fa}"
    out="${FILTERED_DIR}/${base}_filtered.fa"

    # Skip if already filtered
    if [ -f "$out" ]; then
        echo "ℹ️ Skipping already filtered: $file"
        cat "$out" >> "${FILTERED_DIR}/combined_filtered_cds.fa"
        continue
    fi

    echo "🔧 Filtering: $file"
    if seqkit seq -g -m 300 "$file" -o "$out"; then
        if [ -s "$out" ]; then
            cat "$out" >> "${FILTERED_DIR}/combined_filtered_cds.fa"
        else
            echo "⚠️ No sequences passed filter for: $file"
            rm -f "$out"
        fi
    else
        echo "❌ Error filtering: $file" >&2
    fi
done

echo "✅ CDS filtering completed successfully"