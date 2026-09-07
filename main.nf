#!/usr/bin/env nextflow
// Copyright (C) 2024 Genome Research Ltd.

/*
========================================================================================
    HELP
========================================================================================
*/

def logo = NextflowTool.logo(workflow, params.monochrome_logs)

log.info logo

NextflowTool.commandLineParams(workflow.commandLine, log, params.monochrome_logs)


def printHelp() {
    NextflowTool.help_message("${workflow.ProjectDir}/schema.json",
                               ["${workflow.ProjectDir}/assorted-sub-workflows/themisto2/schema.json"],
    params.monochrome_logs, log)
}

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/
//
// SUBWORKFLOWS
//
include { BUILD_COLOR_INDEX } from './assorted-sub-workflows/themisto2/subworkflows/build_color_index.nf'
include { SET_DIFF_CALCULATIONS } from './assorted-sub-workflows/themisto2/subworkflows/setdiff_filter.nf'
include { SBWT_DUMP_UNITIGS } from './assorted-sub-workflows/themisto2/modules/sbwt.nf'
include { POST_PROCESS_MARKERS } from './modules/post_processing_markers.nf'
include { DESIGN_PRIMERS } from './modules/primer3.nf'
include { CHECKPOINT_COUNT } from './assorted-sub-workflows/themisto2/modules/checkpoint_count.nf'
include { MANIFEST_PARSE } from './subworkflows/manifest_parse.nf'


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {

    if (params.help) {
        printHelp()
        exit 0
    }

    // BUILD_COLOR_INDEX takes one pre-paired item per species:
    //   [ [ID: species, target_groups: <csv>], metadata_file, assembly_input ]
    // built from the --manifest TSV (one row per species).
    if (!params.manifest) {
        exit 1, "ERROR: --manifest is required -- a TSV, one row per species, columns " +
                "species / metadata / assemblies / target_groups. See assets/example_manifest.tsv."
    }
    MANIFEST_PARSE(params.manifest)
    samples_ch = MANIFEST_PARSE.out.samples

    BUILD_COLOR_INDEX(samples_ch)

    // step08 -- set-difference filtering. bg_excl (C) = background - species_index (A);
    // markers (G) = candidate_index (E) - bg_excl. candidate_index (E) is only
    // non-empty for species whose meta.target_groups is set; bg_excl is built inside
    // SET_DIFF_CALCULATIONS from --bg_index / --bg_excl_index.
    SET_DIFF_CALCULATIONS(
        BUILD_COLOR_INDEX.out.sbwt_index,
        BUILD_COLOR_INDEX.out.candidate_index
    )

    // Per-stage count checkpoints -> one funnel TSV of this run's own numbers.
    // Accumulate rows from every stage, then collectFile once at the end.
    checkpoint_rows = BUILD_COLOR_INDEX.out.checkpoints
        .mix(SET_DIFF_CALCULATIONS.out.checkpoints)

    // Candidate marker post-processing (step09) -- off by default (see
    // --primer_post_processing's help_text). markers (G) is one .sbwt per
    // species/lineage combo produced by SET_DIFF_CALCULATIONS; dump each to
    // FASTA, then filter/mask for PCR/primer-design suitability.
    if (params.primer_post_processing) {
        SBWT_DUMP_UNITIGS(SET_DIFF_CALCULATIONS.out.markers)
        POST_PROCESS_MARKERS(SBWT_DUMP_UNITIGS.out.unitigs)

        // Checkpoint the post-processing funnel: G dumped to FASTA -> markers
        // passing / rejected by the length + GC filter.
        SBWT_DUMP_UNITIGS.out.unitigs.map    { meta, f -> [meta, 'G_markers_09_dumped_fasta', 'fasta', f] }
        | mix( POST_PROCESS_MARKERS.out.filtered.map { meta, f -> [meta, 'H_markers_09_postproc_pass', 'fasta', f] } )
        | mix( POST_PROCESS_MARKERS.out.rejected.map { meta, f -> [meta, 'H_markers_09_postproc_reject', 'fasta', f] } )
        | set { postproc_checkpoint_inputs }

        CHECKPOINT_COUNT(postproc_checkpoint_inputs)
        checkpoint_rows = checkpoint_rows.mix(CHECKPOINT_COUNT.out.row)

        // Primer3 design (step10) -- off by default, and only meaningful once
        // POST_PROCESS_MARKERS has actually run (it needs the non_designable
        // coordinates from that step's FASTA headers). Runs on the passed/filtered
        // markers only -- rejected markers (too short/global-GC-out-of-range)
        // are never worth designing primers against.
        if (params.primer3_design) {
            DESIGN_PRIMERS(POST_PROCESS_MARKERS.out.filtered)
        }
    }

    // Funnel TSV of this run's own numbers. Publish it only when intermediates
    // are published -- same gate as CHECKPOINT_COUNT's `publishDir enabled:` --
    // so a default run leaves no checkpoints/ dir. Without storeDir the file
    // still collects, it just stays in the work dir.
    def counts_args = [name: 'pipeline_counts.tsv', keepHeader: true, skip: 1, sort: true]
    if( params.publish_intermediate )
        counts_args.storeDir = "${params.outdir}/checkpoints"

    checkpoint_rows
    | map { meta, row -> row }
    | collectFile(counts_args)
    // baitcapture tool TODO create in location: /data/pam/team230/sm71/scratch/gps_project/lsmd/modules/
}
