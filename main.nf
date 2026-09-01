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

    // Input channels from params. TODO: single metadata/assembly pair only --
    // manifest/samplesheet channel for multi-run is PAT-3553 / PAT-3569.
    metadata_ch = Channel.fromPath(params.metadata)
    assembly_ch = Channel.fromPath(params.assembly_input)

    BUILD_COLOR_INDEX(metadata_ch, assembly_ch)
}
