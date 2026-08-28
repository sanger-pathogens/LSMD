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

    // Input channels from params. TODO: single metadata/assembly pair only --
    // manifest/samplesheet channel for multi-run is PAT-3553 / PAT-3569.
    metadata_ch = Channel.fromPath(params.metadata)
    assembly_ch = Channel.fromPath(params.assembly_input)

    BUILD_COLOR_INDEX(metadata_ch, assembly_ch)

    // step08 -- set-difference filtering (D/F/G). lineage_index/candidate_index
    // (B/E) are only non-empty when --target_groups was set; bg_excl (C) is built
    // inside SET_DIFF_CALCULATIONS from --bg_index / --bg_excl_index.
    SET_DIFF_CALCULATIONS(
        BUILD_COLOR_INDEX.out.sbwt_index,
        BUILD_COLOR_INDEX.out.lineage_index,
        BUILD_COLOR_INDEX.out.candidate_index
    )

    // Per-stage count checkpoints -> one funnel TSV of this run's own numbers.
    BUILD_COLOR_INDEX.out.checkpoints
    | mix(SET_DIFF_CALCULATIONS.out.checkpoints)
    | map { meta, row -> row }
    | collectFile(
        name: 'pipeline_counts.tsv',
        storeDir: "${params.outdir}/checkpoints",
        keepHeader: true,
        skip: 1,
        sort: true,
    )

    // Candidate marker post-processing (step09) -- off by default (see
    // --primer_post_processing's help_text). markers (G) is one .sbwt per
    // species/lineage combo produced by SET_DIFF_CALCULATIONS; dump each to
    // FASTA, then filter/mask for PCR/primer-design suitability.
    if (params.primer_post_processing) {
        SBWT_DUMP_UNITIGS(SET_DIFF_CALCULATIONS.out.markers)
        POST_PROCESS_MARKERS(SBWT_DUMP_UNITIGS.out.unitigs)

        // Primer3 design (step10) -- off by default, and only meaningful once
        // POST_PROCESS_MARKERS has actually run (it needs the non_designable
        // coordinates from that step's FASTA headers). Runs on the passed/filtered
        // markers only -- rejected markers (too short/global-GC-out-of-range)
        // are never worth designing primers against.
        if (params.primer3_design) {
            DESIGN_PRIMERS(POST_PROCESS_MARKERS.out.filtered)
        }
    }
    // baitcapture tool TODO create in location: /data/pam/team230/sm71/scratch/gps_project/lsmd/modules/
}
