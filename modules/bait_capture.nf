// Bait-capture panel design (PAT-3586).
//
// Tiles hybridisation-capture baits across the post-processed lineage-specific
// markers using BaitsTools `tilebaits` (https://github.com/campanam/BaitsTools).
//
// This is the bait-capture counterpart to DESIGN_PRIMERS (modules/primer3.nf) on
// the PCR track. It is NOT wired into main.nf yet -- see PAT-3581 (per-assay
// post-processing tracks) for the `--assays pcr,bait` fan-out that will feed it.
//
// Input is the *bait-mode* post-processed marker FASTA: long enough to tile at
// least one bait, low-complexity / repeat windows soft-masked (lowercase) so the
// -K masked-fraction filter can drop baits that sit on them.
//
// BaitsTools applies NO filtration unless the corresponding flag is passed -- the
// "defaults" below are the values BaitsTools uses *when a filter is enabled*, not
// automatic behaviour. We enable GC / masked-fraction / homopolymer / linguistic-
// complexity / complete-length / no-Ns filtering explicitly. Melting-temperature
// filtering is opt-in (params.bait_min_tm / bait_max_tm; null = not applied).
//
// No public BaitsTools conda package or biocontainer exists (checked 2026-09-07:
// not on bioconda, quay.io/biocontainers, Docker Hub or RubyGems.org -- it ships
// only as a GitHub-hosted Ruby gem). containers/baitstools/Dockerfile builds one;
// push it to quay.io/sangerpathogens/baitstools:<ver> and point `container` at
// that (same pattern as themisto2 / sbwt).

process DESIGN_BAITS {
    tag "${meta.ID}"
    label 'cpu_2'
    label 'mem_2'
    label 'time_30m'

    container 'quay.io/sangerpathogens/baitstools:1.8.3'  // TODO(PAT-3586): build + push (containers/baitstools/Dockerfile)

    publishDir mode: 'copy', path: "${params.outdir}/Final_markers/bait/${meta.ID}/"

    input:
    tuple val(meta), path(markers_fasta)

    output:
    tuple val(meta), path("${meta.ID}-baits.fa"),             emit: baits
    tuple val(meta), path("${meta.ID}-filtered-baits.fa"),    emit: filtered,   optional: true
    tuple val(meta), path("${meta.ID}-filtered-params.txt"),  emit: params_tsv, optional: true
    tuple val(meta), path("${meta.ID}-baits.bed"),            emit: bed,        optional: true

    script:
    // Boolean filter flags -- presence is what enables the filter.
    def complete_flag = params.bait_complete      ? "-c" : ""
    def no_ns_flag    = params.bait_no_ns         ? "-N" : ""
    def params_flag   = params.bait_params_table  ? "-w" : ""
    def bed_flag      = params.bait_bed           ? "-B" : ""
    // Melting-temperature filter is opt-in (null = don't pass the flag at all).
    def min_tm_flag   = params.bait_min_tm != null ? "-q ${params.bait_min_tm}" : ""
    def max_tm_flag   = params.bait_max_tm != null ? "-z ${params.bait_max_tm}" : ""
    """
    baitstools tilebaits \\
        -i ${markers_fasta} \\
        -L ${params.bait_length} \\
        -O ${params.bait_offset} \\
        -n ${params.bait_gc_min} \\
        -x ${params.bait_gc_max} \\
        -K ${params.bait_max_mask} \\
        -J ${params.bait_max_homopolymer} \\
        -y ${params.bait_min_complexity} \\
        -T ${params.bait_hyb_type} \\
        -s ${params.bait_na} \\
        -f ${params.bait_formamide} \\
        ${min_tm_flag} \\
        ${max_tm_flag} \\
        ${complete_flag} \\
        ${no_ns_flag} \\
        ${params_flag} \\
        ${bed_flag} \\
        -X ${task.cpus} \\
        -o ${meta.ID} \\
        -Z .
    """
}
