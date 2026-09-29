def VALIDATE_PARAMS() {
    // BUILD_COLOUR_INDEX takes, per species: samples [ [ID: species], metadata_file,
    // assembly_input ], from the --manifest TSV (one row per species), and builds ONLY the
    // species-wide colour index -- no filtering happens there any more (see
    // build_colour_index.nf's header comment).
    if (!params.manifest) {
        exit 1, "ERROR: --manifest is required -- a TSV, one row per species, columns " +
                "species / metadata / assemblies / target_groups / atb_exclude_species. " +
                "See the README's Input section."
    }
    // Primer design is still in progress: stop
    // rather than silently ignore the flag.
    if (params.primer3_design) {
        exit 1, "ERROR: --primer3_design isn't available yet -- primer design is still in progress."
    }
    if (params.marker_post_processing) {
        // Thresholds default values are `null`: it is required that the user picks numeric values for their assay. A bare
        // `--marker_gc_min` (no value) arrives as boolean true, so check for a number.
        def missing = ['marker_min_length', 'marker_gc_min', 'marker_gc_max'].findAll { !(params[it] instanceof Number) }
        if (missing) {
            exit 1, "ERROR: --marker_post_processing needs a numeric value for " +
                    missing.collect { "--${it}" }.join(', ') + " (no defaults)."
        }
        if (params.marker_gc_min >= params.marker_gc_max) {
            exit 1, "ERROR: --marker_gc_min (${params.marker_gc_min}) must be below --marker_gc_max (${params.marker_gc_max})."
        }
    }
}

