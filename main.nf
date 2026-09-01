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

    // Input: either a --manifest TSV (one row per species) or the single-species
    // --metadata / --assembly_input params directly.
    if (params.manifest) {
        MANIFEST_PARSE(params.manifest)
        // TODO (PAT-3569 ASW side): pass MANIFEST_PARSE.out.samples straight through so
        // BUILD_COLOR_INDEX gets meta.id (species -> output folder) and per-species
        // meta.target_groups, instead of splitting it back out into the two channels
        // below and deriving the folder from the metadata basename + a global
        // params.target_groups.
        MANIFEST_PARSE.out.samples
            .multiMap { meta, metadata, assemblies ->
                metadata: metadata
                assembly: assemblies
            }
            .set { mf }
        metadata_ch = mf.metadata
        assembly_ch = mf.assembly
    } else {
        metadata_ch = Channel.fromPath(params.metadata)
        assembly_ch = Channel.fromPath(params.assembly_input)
    }

    BUILD_COLOR_INDEX(metadata_ch, assembly_ch)
}
