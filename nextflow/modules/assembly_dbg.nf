process assembly_dbg {

    tag "assembly_dbg"

    input:
    path input_csv
    val output_dir
    val mode
    val chain

    output:
    path "${output_dir}"

    script:
    """
    mkdir -p ${output_dir}
    echo ">>> Running De Bruijn assembly for \$input_csv"
    python ../src/script_dbg.py --input_csv \$input_csv --folder_outputs ${output_dir} \
        ${mode == 'reference' ? '--reference' : ''} \
        ${chain ? "--chain ${chain}" : ''}
    """
}
