#! /usr/bin/env nextflow

/*
 * InstaNexus pipeline (v0.1)
 * Author: Marco Reverenna
 * Description: example for InstaNexus pipeline using NextFlow.
 */

nextflow.enable.dsl=2

workflow {

    params.input_csv = params.input_csv ?: "inputs/bsa.csv"
    params.chain = params.chain ?: ""
    params.output_dir = params.output_dir ?: "outputs"
    params.reference = params.reference ?: ""

    log.info "Running InstaNexus pipeline with the following parameters:"
    log.info "Input CSV: ${params.input_csv}"
    log.info "Chain: ${params.chain}"
    log.info "Output Directory: ${params.output_dir}"
    log.info "Reference: ${params.reference}"

    // Define the main process
    preprocess(params.input_csv)
    assembly_dbg(params.input_csv, params.chain, params.output_dir, params.reference)
    postprocess(params.output_dir)
}