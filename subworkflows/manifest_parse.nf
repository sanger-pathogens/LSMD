// Parse the run manifest -- a TSV, one row per species, columns:
//
//   species             name used for the output folder (index_species/<species>/ etc.)
//   metadata            path to that species' metadata table (.tsv or .csv)
//   assemblies          path to a directory of assemblies OR a .txt file listing
//                       one assembly path per line
//   target_groups       optional, comma-separated lineage labels for that species
//                       (e.g. "GPSC1,GPSC2"); empty = every lineage with
//                       >= candidate_min_genome_count genomes (excluding "unclassified")
//   atb_target_species  space-separated ATB colour name(s) for the ATB cross-species
//                       check (marker_filtering.nf). Required whenever that check
//                       runs -- NOT derived from 'species', since ATB uses full
//                       binomial names (vibrio_cholerae) that abbreviated manifest
//                       species keys (v_cholerae) don't match, and some species span
//                       more than one ATB colour (streptococcus_pneumoniae splits into
//                       lettered chunks a/b/c/...). Look the name(s) up in ATB's
//                       color_names.txt yourself, e.g.
//                       grep -i '<species>' <atb_color_names> -- marker_filtering.nf's
//                       expand_target_species() picks up clean lettered splits of
//                       whatever you list here automatically, you don't need to
//                       enumerate those yourself.
//
// Emits three channels, all keyed on the SAME slim meta [ID: <species>]:
//   samples             [ [ID: <species>], <metadata file>, <assemblies path> ]
//   target_groups       [ [ID: <species>], <target_groups string> ]
//   atb_target_species  [ [ID: <species>], <atb_target_species string> ]
//
// target_groups / atb_target_species are kept OUT of meta on purpose: meta rides
// through every species-wide build process (COLOR_MAPPING .. THEMISTO2_EXPORT) as
// part of the task hash, but neither is consumed there -- only marker_filtering.nf
// consumes them. Carrying either in meta would make an edit to it invalidate the
// entire species index on -resume. build_color_index.nf never sees them at all;
// marker_filtering.nf joins them back in on the slim meta key.

def parse_manifest_row(row) {
    def species = (row.species ?: "").trim()
    if (!species) {
        error("manifest: every row needs a non-empty 'species' value")
    }
    if (species ==~ /.*[\/\s].*/) {
        error("manifest: 'species' is used as a directory name, so it can't contain "
            + "whitespace or '/' -- got '${species}'")
    }

    def metadata = file((row.metadata ?: "").trim())
    if (!metadata.exists() || metadata.isDirectory()) {
        error("manifest (${species}): 'metadata' must be an existing file -- got '${row.metadata}'")
    }

    def assemblies = file((row.assemblies ?: "").trim())
    if (!assemblies.exists()) {
        error("manifest (${species}): 'assemblies' path does not exist -- got '${row.assemblies}'")
    }

    def target_groups = (row.target_groups ?: "").trim()
    def atb_target_species = (row.atb_target_species ?: "").trim()

    return [[ID: species], metadata, assemblies, target_groups, atb_target_species]
}

workflow MANIFEST_PARSE {
    take:
    manifest   // path to the manifest TSV

    main:
    Channel.fromPath(manifest, checkIfExists: true)
        | splitCsv(header: true, sep: '\t', strip: true)
        | map { row -> parse_manifest_row(row) }
        | multiMap { meta, metadata, assemblies, tg, atb_tgt ->
            samples:            [meta, metadata, assemblies]
            target_groups:      [meta, tg]
            atb_target_species: [meta, atb_tgt]
          }
        | set { parsed }

    emit:
    samples             = parsed.samples             // [ [ID], metadata, assemblies ]
    target_groups       = parsed.target_groups       // [ [ID], target_groups string ]
    atb_target_species  = parsed.atb_target_species  // [ [ID], atb_target_species string ]
}
