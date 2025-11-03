#!/usr/bin/env nextflow

// Pipeline test for single input CSV file bsa.csv

process Preprocessing {

    publishDir "${projectDir}/results/preprocessing", mode: 'copy'

    input:
    path input_csv

    script:
    """
    python ${projectDir}/src/instanexus/preprocessing.py \
      --input-csv ${projectDir}/inputs/bsa.csv \
      --folder-outputs ${projectDir}/outputs \
      --assembly-mode ${params.assembly_mode} \
      --conf ${params.conf} \
      --kmer-size ${params.kmer_size} \
      --size-threshold ${params.size_threshold} \
      --min-overlap ${params.min_overlap}
    """
}

process Assembly {
    publishDir "${projectDir}/results/assembly", mode: 'copy'

    input:
    path cleaned_csv

    script:
    """
    cleaned_dir=\$(dirname \$(dirname "${cleaned_csv}"))
    python ${projectDir}/src/instanexus/assembly.py \
      --input-data "${cleaned_csv}" \
      --output-folder "\$cleaned_dir" \
      --assembly-mode greedy \
      --min-overlap ${params.min_overlap} \
      --size-threshold ${params.size_threshold}
    """
}


params.input = "${projectDir}/inputs/bsa.csv"
params.assembly_mode  = "dbg"
params.conf           = 0.9
params.kmer_size      = 7
params.size_threshold = 12
params.min_overlap    = 5


workflow {
    preprocessing_ch = Channel.fromPath(params.input)
    Preprocessing(preprocessing_ch)

    cleaned_ch = Channel.fromPath("${projectDir}/outputs/*/*/cleaned/cleaned_data.csv")
    Assembly(cleaned_ch)
}
