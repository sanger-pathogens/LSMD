// Parse the run manifest -- a TSV (tab-separated; a .csv or a comma-separated header
// stops the run), one row per species, exactly these five columns:
//
//   species             the species' ATB colour name, e.g. streptococcus_pneumoniae or
//                       vibrio_cholerae. Used as the output folder/file prefix
//                       (index_species/<species>/ etc.), so no whitespace or '/', and as
//                       the target of the ATB cross-species check (marker_filtering.nf).
//                       ATB's lettered splits of a species (streptococcus_pneumoniaea,
//                       ...b, ...) are picked up automatically by expand_target_species().
//                       A species that isn't in --atb_color_names skips the ATB check with
//                       a warning (and a "did you mean" suggestion, in case it's a typo).
//   metadata            path to that species' metadata table
//   assemblies          path to a directory of assemblies OR a .txt file listing
//                       one assembly path per line
//   target_groups       optional, comma-separated lineage labels for that species
//                       (e.g. "GPSC1,GPSC2"); empty = every lineage with
//                       >= candidate_min_genome_count genomes (excluding "unclassified")
//   atb_exclude_species optional, comma-separated ATB colour name(s) dropped from the
//                       ATB check's max-outside side entirely (e.g. close relatives
//                       ATB can't reliably separate from the target). 'unknown' (ATB's
//                       catch-all bucket for unassigned/low-confidence genomes) is always
//                       excluded; values here are added on top of it.
//
// Emits four channels, all keyed on the SAME slim meta [ID: <species>]:
//   samples             [ [ID: <species>], <metadata file>, <assemblies path> ]
//   target_groups       [ [ID: <species>], <target_groups string> ]
//   atb_target_species  [ [ID: <species>], <species, or '' when it isn't in ATB> ]
//   atb_exclude_species [ [ID: <species>], <space-separated ATB colour names> ]
//
// target_groups / atb_target_species / atb_exclude_species are kept OUT of meta on purpose: meta rides
// through every species-wide build process (COLOR_MAPPING .. THEMISTO2_EXPORT) as
// part of the task hash, but neither is consumed there -- only marker_filtering.nf
// consumes them. Carrying either in meta would make an edit to it invalidate the
// entire species index on -resume. build_color_index.nf never sees them at all;
// marker_filtering.nf joins them back in on the slim meta key.

def manifest_columns() {
    return ['species', 'metadata', 'assemblies', 'target_groups', 'atb_exclude_species']
}

// Stop before any rows are read if the manifest isn't a 5-column TSV. splitCsv on a
// comma-separated file doesn't fail by itself: it reads each line as one column.
def check_manifest_format(manifest_file) {
    if (manifest_file.name.toLowerCase().endsWith('.csv')) {
        error("manifest: must be a tab-separated .tsv file, not a .csv -- got '${manifest_file}'")
    }
    def header = manifest_file.withReader { it.readLine() } ?: ''
    // -1 keeps trailing empty fields, so a trailing tab shows up as a blank column below.
    // trim() also drops a Windows '\r'.
    def columns = header.split('\t', -1).collect { it.trim() }
    def expected = manifest_columns()
    if (columns.size() == 1 && header.contains(',')) {
        error("manifest: the header is comma-separated; the manifest must be tab-separated -- got '${manifest_file}'")
    }
    // splitCsv rejects blank header names with an unhelpful error, so catch them here.
    if (columns.any { !it }) {
        error("manifest: the header has a blank column name (often a trailing tab, e.g. from Excel). "
            + "Remove it -- got '${manifest_file}'")
    }
    def missing = expected - columns
    def extra = columns - expected
    if (missing || extra) {
        error("manifest: the header must have exactly these tab-separated columns: ${expected.join(', ')}."
            + (missing ? " Missing: ${missing.join(', ')}." : '')
            + (extra ? " Not recognised: ${extra.join(', ')}." : ''))
    }
}

// Bare ATB colour names from color_names.txt ('<id>\tper_species_unitigs/<name>-unitigs-k31.fna'),
// unwrapped the same way atb_cross_species_filter.py's _clean_atb_name() does.
def load_atb_names(color_names_path) {
    return file(color_names_path, checkIfExists: true).readLines()
        .findAll { it.trim() }
        .collect { line ->
            def name = line.contains('\t') ? line.split('\t', 2)[1] : line
            name = name.trim()
            if (name.startsWith('per_species_unitigs/')) { name = name - 'per_species_unitigs/' }
            if (name.endsWith('-unitigs-k31.fna')) { name = name[0..<(name.length() - '-unitigs-k31.fna'.length())] }
            name
        } as Set
}

// True if ATB has the species itself or a clean lettered split of it (same rule as
// expand_target_species(): base + exactly one trailing lowercase letter).
def in_atb(species, atb_names) {
    return atb_names.contains(species) || atb_names.any { it.length() == species.length() + 1 && it.startsWith(species) && it[-1] ==~ /[a-z]/ }
}

def edit_distance(String a, String b) {
    def prev = (0..b.length()).toList()
    for (int i = 1; i <= a.length(); i++) {
        def cur = [i]
        for (int j = 1; j <= b.length(); j++) {
            cur << [prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] == b[j - 1] ? 0 : 1)].min()
        }
        prev = cur
    }
    return prev[b.length()]
}

// Suggestions skip ATB's split chunks (a name that's another name plus 1-2 trailing
// letters, e.g. streptococcus_mitisbw), so they offer base species names.
def closest_atb_names(species, atb_names) {
    return atb_names
        .findAll { n -> !(1..2).any { k -> n.length() > k && atb_names.contains(n[0..<(n.length() - k)]) } }
        .collect { [it, edit_distance(species, it)] }
        .findAll { it[1] <= 3 }
        .sort { it[1] }
        .take(3)
        .collect { it[0] }
}

def parse_manifest_row(row, atb_names) {
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

    def atb_target_species = species
    if (!in_atb(species, atb_names)) {
        def suggestions = closest_atb_names(species, atb_names)
        log.warn("manifest (${species}): '${species}' isn't an ATB colour name in ${params.atb_color_names}, so its "
            + "markers will skip the ATB cross-species check UNVERIFIED."
            + (suggestions ? " Did you mean: ${suggestions.join(', ')}?" : ""))
        atb_target_species = ""
    }

    def excluded = (row.atb_exclude_species ?: "").split(',').collect { it.trim() }.findAll { it }
    excluded.findAll { !in_atb(it, atb_names) }.each { name ->
        def suggestions = closest_atb_names(name, atb_names)
        log.warn("manifest (${species}): atb_exclude_species '${name}' couldn't be found in ${params.atb_color_names}, "
            + "so it has no effect. Please check the spelling."
            + (suggestions ? " Closest matches: ${suggestions.join(', ')}." : ""))
    }
    def atb_exclude_species = (['unknown'] + excluded).unique().join(' ')

    return [[ID: species], metadata, assemblies, target_groups, atb_target_species, atb_exclude_species]
}

workflow MANIFEST_PARSE {
    take:
    manifest   // path to the manifest TSV

    main:
    def manifest_file = file(manifest, checkIfExists: true)
    check_manifest_format(manifest_file)
    def atb_names = load_atb_names(params.atb_color_names)

    Channel.fromPath(manifest_file)
        | splitCsv(header: true, sep: '\t', strip: true)
        | map { row -> parse_manifest_row(row, atb_names) }
        | multiMap { meta, metadata, assemblies, tg, atb_tgt, atb_excl ->
            samples:             [meta, metadata, assemblies]
            target_groups:       [meta, tg]
            atb_target_species:  [meta, atb_tgt]
            atb_exclude_species: [meta, atb_excl]
          }
        | set { parsed }

    emit:
    samples             = parsed.samples             // [ [ID], metadata, assemblies ]
    target_groups       = parsed.target_groups       // [ [ID], target_groups string ]
    atb_target_species  = parsed.atb_target_species  // [ [ID], species, or '' when not in ATB ]
    atb_exclude_species = parsed.atb_exclude_species // [ [ID], space-separated ATB colour names ]
}
