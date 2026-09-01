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
    // Either from a --manifest TSV (one row per species) or, for a single species,
    // built from the --metadata / --assembly_input / --target_groups params.
    if (params.manifest) {
        MANIFEST_PARSE(params.manifest)
        samples_ch = MANIFEST_PARSE.out.samples
    } else {
        // Single species: ID defaults to the metadata file's basename.
        samples_ch = Channel.of([
            [ID: file(params.metadata).baseName, target_groups: params.target_groups ?: ''],
            file(params.metadata),
            file(params.assembly_input),
        ])
    }

    BUILD_COLOR_INDEX(samples_ch)
}
