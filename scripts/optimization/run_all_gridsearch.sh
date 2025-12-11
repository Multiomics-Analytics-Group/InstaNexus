#!/bin/bash

MODES=("greedy" "dbg_weighted" "multimodal_dbg")
METADATA="json/sample_metadata.json"
GRID_PARAMS="json/gridsearch_params.json"
WORKERS=6

# antibody samples (Light and Heavy chains)
#ANTIBODIES=("ma1" "ma2" "ma3") 

ANTIBODIES=("ma1") 

# single chain (Nanobodies, Binders, BSA, etc.)
#OTHERS=("nb1" "nb2" "nb3" "nb4" "nb5" "nb6" "nb7" "nb8" "nb9" "nb10" "bsa" "bind1" "bind2" "bind3")
OTHERS=("nb1" "bsa" "bind1")

run_grid() {
    local run=$1
    local chain=$2
    local mode=$3
    
    local input_file="inputs/cleaned/${run}_cleaned.csv"
    
    local output_dir="outputs/_grid_search/${mode}"  

    if [ -f "$input_file" ]; then
        echo ">>> [${mode}] Launching: ${run} (Chain: '${chain}')"
        
        python scripts/optimization/grid_search.py \
            --input-csv "$input_file" \
            --metadata-json "$METADATA" \
            --grid-json "$GRID_PARAMS" \
            --mode "$mode" \
            --chain "$chain" \
            --workers "$WORKERS" \
            > "logs/grid_${mode}_${run}_${chain}.log" 2>&1
            
    else
        echo "WARNING: Input file for ${run} not found. Skipping."
    fi
}

mkdir -p logs

for mode in "${MODES[@]}"; do
    echo "========================================"
    echo " PROCESSING MODE: $mode "
    echo "========================================"

    for sample in "${ANTIBODIES[@]}"; do
        run_grid "$sample" "light" "$mode"
        run_grid "$sample" "heavy" "$mode"
    done

    for sample in "${OTHERS[@]}"; do
        run_grid "$sample" "" "$mode"
    done
    
done

echo "----------------------------------------"
echo "All grid search jobs completed."