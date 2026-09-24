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
include { BUILD_COLOUR_INDEX } from './assorted-sub-workflows/themisto2/subworkflows/build_colour_index.nf'
include { MARKER_FILTERING } from './assorted-sub-workflows/themisto2/subworkflows/marker_filtering.nf'
include { POST_PROCESS_MARKERS } from './modules/post_processing_markers.nf'
include { DESIGN_PRIMERS } from './modules/primer3.nf'
include { CHECKPOINT_FASTA } from './assorted-sub-workflows/themisto2/modules/checkpoint.nf'
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

    // BUILD_COLOUR_INDEX takes, per species: samples [ [ID: species], metadata_file,
    // assembly_input ], from the --manifest TSV (one row per species), and builds ONLY the
    // species-wide colour index -- no filtering happens there any more (see
    // build_colour_index.nf's header comment).
    if (!params.manifest) {
        exit 1, "ERROR: --manifest is required -- a TSV, one row per species, columns " +
                "species / metadata / assemblies / target_groups / atb_exclude_species. " +
                "See the README's Input section."
    }
    if (params.primer3_design && !params.marker_post_processing) {
        exit 1, "ERROR: --primer3_design needs --marker_post_processing (primer3 runs on the post-processed markers)."
    }
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

    // Per-stage count checkpoints -> one funnel TSV of this run's own numbers.
    // Accumulate rows from every stage, then collectFile once at the end.
    checkpoint_rows = BUILD_COLOUR_INDEX.out.checkpoints
        .mix(MARKER_FILTERING.out.checkpoints)

    // Candidate marker post-processing -- off by default (see
    // --marker_post_processing's help_text). MARKER_FILTERING.out.markers is already one
    // FASTA per species/lineage combo (ATB-checked, or unchecked with a warning -- see
    // marker_filtering.nf); filter/mask it for downstream assay design (primer3 and/or baits).
    if (params.marker_post_processing) {
        POST_PROCESS_MARKERS(MARKER_FILTERING.out.markers)

        POST_PROCESS_MARKERS.out.filtered.map { meta, f -> [meta, 100, 'markers_postproc_pass', 'fasta', f] }
        | mix( POST_PROCESS_MARKERS.out.rejected.map { meta, f -> [meta, 110, 'markers_postproc_reject', 'fasta', f] } )
        | set { postproc_checkpoint_inputs }

        CHECKPOINT_FASTA(postproc_checkpoint_inputs)
        checkpoint_rows = checkpoint_rows.mix(CHECKPOINT_FASTA.out.row)

        // Primer3 design -- off by default, and only meaningful once
        // POST_PROCESS_MARKERS has actually run (it needs the non_designable
        // coordinates from that step's FASTA headers). Runs on the passed/filtered
        // markers only -- rejected markers (too short/global-GC-out-of-range)
        // are never worth designing primers against.
        if (params.primer3_design) {
            DESIGN_PRIMERS(POST_PROCESS_MARKERS.out.filtered)
        }
    }

    def counts_args = [name: 'pipeline_counts.tsv', keepHeader: true, skip: 1, sort: true]
    if( params.publish_intermediate )
        counts_args.storeDir = "${params.outdir}/checkpoints"

    checkpoint_rows
    | map { meta, row -> row }
    | collectFile(counts_args)
    // baitcapture tool TODO create in location: /data/pam/team230/sm71/scratch/gps_project/lsmd/modules/
}
