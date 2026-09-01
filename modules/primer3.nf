process DESIGN_PRIMERS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    container 'quay.io/biocontainers/primer3:2.6.1--pl5321h503566f_7'

    publishDir mode: 'copy', path: "${params.outdir}/candidate_markers/${meta.ID}/primers/"

    input:
    tuple val(meta), path(filtered_fasta)

    output:
    tuple val(meta), path(primers_tsv),     emit: primers
    tuple val(meta), path(no_primers_tsv),  emit: no_primers

    script:
    primers_tsv = "${meta.ID}_primers.tsv"
    no_primers_tsv = "${meta.ID}_no_primers.tsv"
    """
    ${moduleDir}/../bin/design_primers.py \\
        ${filtered_fasta} \\
        --label ${meta.ID} \\
        --out-dir . \\
        --product-size-range ${params.primer3_product_size_range} \\
        --num-return ${params.primer3_num_return} \\
        --opt-size ${params.primer3_opt_size} \\
        --min-size ${params.primer3_min_size} \\
        --max-size ${params.primer3_max_size} \\
        --opt-tm ${params.primer3_opt_tm} \\
        --min-tm ${params.primer3_min_tm} \\
        --max-tm ${params.primer3_max_tm} \\
        --min-gc ${params.primer3_min_gc} \\
        --max-gc ${params.primer3_max_gc}
    """
}
