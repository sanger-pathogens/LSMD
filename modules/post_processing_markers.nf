process POST_PROCESS_MARKERS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/${meta.species}/post_processed_markers/"

    input:
    tuple val(meta), path(unitigs_fasta)

    output:
    tuple val(meta), path(filtered_fasta), emit: filtered
    tuple val(meta), path(rejected_fasta), emit: rejected, optional: true
    tuple val(meta), path(plot_png),       emit: plot,     optional: true

    script:
    def id = "${meta.species}_${meta.ID}"
    filtered_fasta = "${id}_markers.fasta"
    rejected_fasta = "${id}_rejected_markers.fasta"
    plot_png = "${id}_marker_analysis.png"
    def reject_flag = params.marker_write_rejected ? "-r" : ""
    def plot_flag = params.marker_plot ? "--plot" : ""
    """
    ${moduleDir}/../bin/post_processing_unitigs.py \\
        ${unitigs_fasta} \\
        -o . \\
        -l ${params.marker_min_length} \\
        -g ${params.marker_gc_min} \\
        -G ${params.marker_gc_max} \\
        -w ${params.marker_window_size} \\
        ${reject_flag} \\
        ${plot_flag}

    mv filtered_unitigs.fasta ${filtered_fasta}
    [ -f rejected_unitigs.fasta ] && mv rejected_unitigs.fasta ${rejected_fasta}
    [ -f unitig_analysis.png ] && mv unitig_analysis.png ${plot_png}
    true
    """
}
