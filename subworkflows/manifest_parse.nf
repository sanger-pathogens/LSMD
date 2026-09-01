// Parse the run manifest -- a TSV, one row per species, columns:
//
//   species        name used for the output folder (index_species/<species>/ etc.)
//   metadata       path to that species' metadata table (.tsv or .csv)
//   assemblies     path to a directory of assemblies OR a .txt file listing
//                  one assembly path per line
//   target_groups  optional, comma-separated lineage labels for that species
//                  (e.g. "GPSC1,GPSC2"); empty = species-wide only
//
// Emits one tuple per row:
//   [ [id: <species>, target_groups: <string>], <metadata file>, <assemblies path> ]
//
// TODO (PAT-3569 ASW side): BUILD_COLOR_INDEX still takes two separate channels and
// combine()s them, so it only handles one species per run and derives its folder
// name from the metadata basename, not `id`. Until it takes this pre-paired channel,
// MANIFEST_PARSE enforces a single row -- see the guard in main.nf.

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

    return [[id: species, target_groups: target_groups], metadata, assemblies]
}

workflow MANIFEST_PARSE {
    take:
    manifest   // path to the manifest TSV

    main:
    Channel.fromPath(manifest, checkIfExists: true)
        | splitCsv(header: true, sep: '\t', strip: true)
        | toList
        | map { rows ->
            if (rows.isEmpty()) {
                error("manifest has no data rows: ${manifest}")
            }
            // Single-species guard -- drop once the ASW side takes a per-species channel.
            if (rows.size() > 1) {
                error("manifest has ${rows.size()} rows -- multi-species runs need the "
                    + "PAT-3569 ASW change (BUILD_COLOR_INDEX taking a pre-paired "
                    + "channel). Use a one-row manifest for now.")
            }
            rows
        }
        | flatMap { it }
        | map { row -> parse_manifest_row(row) }
        | set { samples }

    emit:
    samples   // [ [id, target_groups], metadata, assemblies ]
}
