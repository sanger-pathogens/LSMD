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

    // build the real channels from your params — this is the "COMBINE_IRODS" step for you
    metadata_ch = Channel.fromPath(params.metadata)
    assembly_ch = Channel.fromPath(params.assembly_input)

    // TODO: next step in development -- this only supports a single metadata/assembly
    // pair per run. To support multiple runs, replace these two params-derived channels
    // with a manifest/samplesheet channel (see combined_input.nf's manifest_of_lanes
    // pattern for precedent) and update the call below to pass that single channel in.

    // pipe them in — matches take: ch_metadata / ch_assembly
    BUILD_COLOR_INDEX(metadata_ch, assembly_ch)
}
