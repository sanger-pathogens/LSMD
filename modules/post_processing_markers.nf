process POST_PROCESS_MARKERS {
    tag "${meta.ID}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_30m'

    // NOTE: previously used quay.io/sangerpathogens/pandas:2.2.1 on the assumption
    // it already had Biopython -- it doesn't (confirmed via ModuleNotFoundError: No
    // module named 'Bio' on 2026-08-22, test_3/test_4 runs). post_processing_unitigs.py
    // only hard-requires Biopython (Bio.SeqUtils.gc_fraction as of 2026-08-25 -- the
    // old Bio.SeqUtils.GC() this used to call was removed upstream of 1.84); matplotlib
    // is optional and degrades gracefully (try/except) when --plot is requested but
    // unavailable.
    container 'quay.io/biocontainers/biopython:1.84'

    publishDir mode: 'copy', path: "${params.outdir}/Final_markers/${meta.ID}/"

    input:
    tuple val(meta), path(unitigs_fasta)

    output:
    tuple val(meta), path(filtered_fasta), emit: filtered
    tuple val(meta), path(rejected_fasta), emit: rejected, optional: true
    tuple val(meta), path(plot_png),       emit: plot,     optional: true

    script:
    // Nextflow DSL2 won't let a bare output-bound variable (rejected_fasta,
    // plot_png -- needed as bare names for output: path(...) to bind to) be
    // read again inside another expression in the same script block ("already
    // defined in the process scope") -- so the flag strings below are built
    // from their own def locals instead of re-referencing those bare vars.
    //
    // post_processing_unitigs.py's real CLI (confirmed 2026-08-25, it does not
    // match what used to be called here): -o takes an OUTPUT DIRECTORY, not a
    // filename -- fixed names (filtered_unitigs.fasta/rejected_unitigs.fasta/
    // unitig_analysis.png) land inside it. -r/--write-rejected and --plot are
    // both bare boolean flags with no filename argument (there's no -p at all).
    // So: run it against '.' and rename its fixed-name outputs to the
    // meta.ID-prefixed names declared below.
    filtered_fasta = "${meta.ID}_filtered_markers.fasta"
    rejected_fasta = "${meta.ID}_rejected_markers.fasta"
    plot_png = "${meta.ID}_marker_analysis.png"
    def reject_flag = params.primer_write_rejected ? "-r" : ""
    def plot_flag = params.primer_plot ? "--plot" : ""
    """
    ${moduleDir}/../bin/post_processing_unitigs.py \\
        ${unitigs_fasta} \\
        -o . \\
        -l ${params.primer_min_length} \\
        -g ${params.primer_gc_min} \\
        -G ${params.primer_gc_max} \\
        -w ${params.primer_window_size} \\
        ${reject_flag} \\
        ${plot_flag}

    mv filtered_unitigs.fasta ${filtered_fasta}
    [ -f rejected_unitigs.fasta ] && mv rejected_unitigs.fasta ${rejected_fasta}
    [ -f unitig_analysis.png ] && mv unitig_analysis.png ${plot_png}
    true
    """
}
