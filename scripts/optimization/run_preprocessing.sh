#!/bin/bash

mkdir -p inputs/cleaned

ANTIBODIES=("ma1" "ma2" "ma3")

echo "--- Starting Batch Preprocessing ---"

for file in inputs/*.csv; do
    
    if [[ "$file" == *"quant"* ]] || [[ "$file" == *"_cleaned"* ]]; then
        continue
    fi

    filename=$(basename "$file" .csv)
    echo ">>> Processing: $filename"
    
    CHAIN_FLAG=""
    
    for ab in "${ANTIBODIES[@]}"; do
        if [[ "$filename" == "$ab" ]]; then
            CHAIN_FLAG="--chain light"
            break
        fi
    done
    
    python -m instanexus.preprocessing \
        --input-csv "$file" \
        --metadata-json json/sample_metadata.json \
        --contaminants-fasta fasta/contaminants.fasta \
        --output-csv-path "inputs/cleaned/${filename}_cleaned.csv" \
        --fdr 1.0 \
        --conf 0.0 \
        $CHAIN_FLAG

done

echo "--- Preprocessing Complete. Files are in inputs/cleaned/ ---"