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
//   atb_exclude_species optional, space-separated ATB colour name(s) dropped from the
//                       ATB check's max-outside side entirely (e.g. close relatives
//                       ATB can't reliably separate from the target). Blank = 'unknown'
//                       (ATB's catch-all bucket for unassigned/low-confidence genomes).
//
// Optional group-label cleaning columns (see color_mapping.py / the themisto2 README):
//   label_missing       '|'-separated label values that count as missing (-> unclassified),
//                       case-insensitive. Blank = color_mapping.py's default list.
//   label_multi         keep | smallest | unclassified, for labels containing ';'. Blank = keep.
//   label_map           path to a raw_label -> group TSV. Blank = no map. Checked here, so a
//                       typo stops the run before any jobs start.
//   unclassified_genomes keep | drop -- whether genomes whose label ends up 'unclassified'
//                       stay in the index. Blank = keep.
//
// Emits four channels, all keyed on the SAME slim meta [ID: <species>]:
//   samples             [ [ID: <species>], <metadata file>, <assemblies path>,
//                         label_missing, label_multi, <label_map file>, unclassified_genomes ]
//   target_groups       [ [ID: <species>], <target_groups string> ]
//   atb_target_species  [ [ID: <species>], <atb_target_species string> ]
//   atb_exclude_species [ [ID: <species>], <atb_exclude_species string> ]
//
// The label-cleaning columns ride in samples on purpose: they change which colour each
// genome gets, so editing one must rebuild that species' index.
//
// target_groups / atb_target_species / atb_exclude_species are kept OUT of meta on purpose: meta rides
// through every species-wide build process (COLOR_MAPPING .. THEMISTO2_EXPORT) as
// part of the task hash, but neither is consumed there -- only marker_filtering.nf
// consumes them. Carrying either in meta would make an edit to it invalidate the
// entire species index on -resume. build_color_index.nf never sees them at all;
// marker_filtering.nf joins them back in on the slim meta key.

// Staged for species with no label_map (path inputs can't be empty); COLOR_MAPPING
// recognises it by name and never passes it to the script.
def no_label_map() {
    return "${projectDir}/assorted-sub-workflows/themisto2/assets/NO_LABEL_MAP"
}

def parse_choice(row, species, column, allowed) {
    def value = (row[column] ?: "").trim() ?: allowed[0]
    if (!(value in allowed)) {
        error("manifest (${species}): '${column}' must be one of ${allowed.join(' / ')} (or blank) -- got '${row[column]}'")
    }
    return value
}

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
    def atb_exclude_species = (row.atb_exclude_species ?: "").trim() ?: "unknown"

    def label_missing = (row.label_missing ?: "").trim()
    def label_multi = parse_choice(row, species, 'label_multi', ['keep', 'smallest', 'unclassified'])
    def unclassified = parse_choice(row, species, 'unclassified_genomes', ['keep', 'drop'])
    def label_map_path = (row.label_map ?: "").trim()
    def label_map = file(label_map_path ?: no_label_map())
    if (!label_map.exists() || label_map.isDirectory()) {
        error("manifest (${species}): 'label_map' must be an existing file -- got '${row.label_map}'")
    }

    return [[ID: species], metadata, assemblies, target_groups, atb_target_species, atb_exclude_species,
            label_missing, label_multi, label_map, unclassified]
}

workflow MANIFEST_PARSE {
    take:
    manifest   // path to the manifest TSV

    main:
    Channel.fromPath(manifest, checkIfExists: true)
        | splitCsv(header: true, sep: '\t', strip: true)
        | map { row -> parse_manifest_row(row) }
        | multiMap { meta, metadata, assemblies, tg, atb_tgt, atb_excl, lbl_missing, lbl_multi, lbl_map, uncl ->
            samples:             [meta, metadata, assemblies, lbl_missing, lbl_multi, lbl_map, uncl]
            target_groups:       [meta, tg]
            atb_target_species:  [meta, atb_tgt]
            atb_exclude_species: [meta, atb_excl]
          }
        | set { parsed }

    emit:
    samples             = parsed.samples             // [ [ID], metadata, assemblies, label_missing, label_multi, label_map, unclassified_genomes ]
    target_groups       = parsed.target_groups       // [ [ID], target_groups string ]
    atb_target_species  = parsed.atb_target_species  // [ [ID], atb_target_species string ]
    atb_exclude_species = parsed.atb_exclude_species // [ [ID], atb_exclude_species string ]
}
