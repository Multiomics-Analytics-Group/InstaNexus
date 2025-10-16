process preprocess {

    input:
    path input_csv

    output:
    path "${output_dir}/preprocessed"

    script:
    """
    echo "Preprocessing input CSV: ${input_csv}"
    python ../src/preprocessing.py --input_csv ${input_csv} --chain ${chain} --folder_outputs ${output_dir}/preprocessed --reference ${reference}
    """
}