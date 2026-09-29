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
include { BUILD_COLOUR_INDEX                    } from './assorted-sub-workflows/themisto2/subworkflows/build_colour_index.nf'
include { MARKER_FILTERING                      } from './assorted-sub-workflows/themisto2/subworkflows/marker_filtering.nf'
include { POST_PROCESS_MARKERS                  } from './modules/post_processing_markers.nf'
include { VALIDATE_PARAMS                       } from './modules/validate_parameters.nf'
// In progress, not wired in yet: primer design (modules/primer3.nf) and bait design
// (modules/bait_capture.nf).
// include { DESIGN_PRIMERS } from './modules/primer3.nf'
include { CHECKPOINT_FASTA; CHECKPOINT_REPORT   } from './assorted-sub-workflows/themisto2/modules/checkpoint.nf'
include { MANIFEST_PARSE                        } from './subworkflows/manifest_parse.nf'


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

    VALIDATE_PARAMS()
    
    MANIFEST_PARSE(params.manifest)
    samples_ch = MANIFEST_PARSE.out.samples

    BUILD_COLOUR_INDEX(samples_ch)

    // Lineage-specificity filtering, candidate index rebuild, and the ATB cross-species
    // check (replaces the old bg_excl/markers sbwt set-diff) -- all keyed off the manifest's
    // species (as the ATB target) / target_groups / atb_exclude_species columns, kept out of meta
    // upstream so editing either doesn't bust the species index cache (see manifest_parse.nf).
    MARKER_FILTERING(
        BUILD_COLOUR_INDEX.out.species_export,
        MANIFEST_PARSE.out.target_groups,
        MANIFEST_PARSE.out.atb_target_species,
        MANIFEST_PARSE.out.atb_exclude_species
    )

    // Per-stage count checkpoints -> one pipeline_counts.tsv per species.
    // Accumulate rows from every stage, then CHECKPOINT_REPORT once per species at the end.
    checkpoint_rows = BUILD_COLOUR_INDEX.out.checkpoints
        .mix(MARKER_FILTERING.out.checkpoints)

    // Candidate marker post-processing -- off by default (see
    // --marker_post_processing's help_text). MARKER_FILTERING.out.markers is already one
    // FASTA per species/lineage combo (ATB-checked, or unchecked with a warning -- see
    // marker_filtering.nf); filter/mask it for downstream assay design.
    if (params.marker_post_processing) {
        POST_PROCESS_MARKERS(MARKER_FILTERING.out.markers)

        POST_PROCESS_MARKERS.out.filtered.map { meta, f -> [meta, 'markers_postproc_pass', 'fasta', f] }
        | mix( POST_PROCESS_MARKERS.out.rejected.map { meta, f -> [meta, 'markers_postproc_reject', 'fasta', f] } )
        | set { postproc_checkpoint_inputs }

        CHECKPOINT_FASTA(postproc_checkpoint_inputs)
        checkpoint_rows = checkpoint_rows.mix(CHECKPOINT_FASTA.out.row)

        // Primer design -- in progress, not wired in yet. When it is, it runs on the
        // passed/filtered markers only (it needs the soft-masking from this step).
        // if (params.primer3_design) {
        //     DESIGN_PRIMERS(POST_PROCESS_MARKERS.out.filtered)
        // }
    }

    // Rows sort into pipeline order inside CHECKPOINT_REPORT (see checkpoint_steps() in
    // checkpoint.nf); always published to results/<species>/checkpoint/pipeline_counts.tsv.
    checkpoint_rows
    | map { meta, row -> [meta.species ?: meta.ID, row] }
    | groupTuple
    | CHECKPOINT_REPORT
}
