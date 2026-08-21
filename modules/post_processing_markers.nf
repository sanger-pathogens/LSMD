process POST_PROCESS_MARKERS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    // Same image as color_mapping.nf/lineage_index_filtering.nf -- already has
    // Biopython (1.87) and matplotlib (3.10.9), no separate container needed.
    container 'quay.io/sangerpathogens/pandas:2.2.1'

    publishDir mode: 'copy', path: "${params.outdir}/candidate_markers/${meta.ID}/"

    input:
    tuple val(meta), path(unitigs_fasta)

    output:
    tuple val(meta), path(filtered_fasta), emit: filtered
    tuple val(meta), path(rejected_fasta), emit: rejected, optional: true
    tuple val(meta), path(plot_png),       emit: plot,     optional: true

    script:
    filtered_fasta = "${meta.ID}_filtered_markers.fasta"
    rejected_fasta = "${meta.ID}_rejected_markers.fasta"
    plot_png = "${meta.ID}_marker_analysis.png"
    def reject_flag = params.primer_write_rejected ? "-r ${rejected_fasta}" : ""
    def plot_flag = params.primer_plot ? "--plot -p ${plot_png}" : ""
    """
    ${moduleDir}/../bin/post_processing_unitigs.py \\
        ${unitigs_fasta} \\
        -o ${filtered_fasta} \\
        -l ${params.primer_min_length} \\
        -g ${params.primer_gc_min} \\
        -G ${params.primer_gc_max} \\
        -w ${params.primer_window_size} \\
        ${reject_flag} \\
        ${plot_flag}
    """
}
