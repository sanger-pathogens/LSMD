// Parse the run manifest -- a TSV, one row per species, columns:
//
//   species        name used for the output folder (index_species/<species>/ etc.)
//   metadata       path to that species' metadata table (.tsv or .csv)
//   assemblies     path to a directory of assemblies OR a .txt file listing
//                  one assembly path per line
//   target_groups  optional, comma-separated lineage labels for that species
//                  (e.g. "GPSC1,GPSC2"); empty = species-wide only
//
// Emits one tuple per row, the shape BUILD_COLOR_INDEX's samples_ch expects:
//   [ [ID: <species>, target_groups: <string>], <metadata file>, <assemblies path> ]

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

    return [[ID: species, target_groups: target_groups], metadata, assemblies]
}

workflow MANIFEST_PARSE {
    take:
    manifest   // path to the manifest TSV

    main:
    Channel.fromPath(manifest, checkIfExists: true)
        | splitCsv(header: true, sep: '\t', strip: true)
        | map { row -> parse_manifest_row(row) }
        | set { samples }

    emit:
    samples   // [ [ID, target_groups], metadata, assemblies ]
}
