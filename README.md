# LSMD: Lineage-Specific Marker Discovery

## What is this?

LSMD discovers k-mer markers that distinguish a target group (lineage, sublineage, serotype, etc.) from related groups within a species. It combines the accuracy of whole-genome analysis with the speed of k-mer-based approaches — discovering markers in hours rather than weeks.

The pipeline was developed and tested for **lineage-level grouping** (e.g. GPSC for _S. pneumoniae_, 7PET for _V. cholerae_), but is flexible enough to work with any grouping you define in your metadata — serotypes, clades, resistance phenotypes, or any other categorical column. You provide the genome assemblies and define the groupings; the pipeline discovers markers that distinguish your chosen groups.

Throughout this README, **group** means whatever categorical column you point `--group_label` at (GPSC, lineage, serotype, …). "Lineage" is just the most common example.

[[_TOC_]]

## Why k-mers?

K-mer analysis enables LSMD to achieve what older discovery methods cannot:

- **Alignment-free**: No need to align genomes or call variants — faster and more objective
- **Annotation-independent**: Works with non-coding sequences; doesn't require gene calls or functional annotation
- **Unbiased**: Uses all k-mers in the input genomes, not a subset (unlike sketching/sub-sampling approaches)
- **Highly scalable**: Tested on 40,000+ genomes; computational cost scales with data volume, not complexity
- **Deployment-ready**: Output k-mers are short, concrete sequences suitable for immediate PCR primer design or bait-capture panel construction

The pipeline uses **SBWT** for compact k-mer indexing and **Themisto2** for colour-mapped k-mer indexes and pseudoalignment.

## Quick start

**Prerequisites:** Nextflow, Docker or Singularity

```bash
git clone --recurse-submodules git@gitlab.internal.sanger.ac.uk:sanger-pathogens/pipelines/lsmd.git
cd lsmd
```

Then prepare a manifest and run it: see [Running your own data](#running-your-own-data).

## Running your own data

### Prerequisites

- Nextflow ≥ 21.04.0
- Docker or Singularity
- A manifest TSV (one row per species; see [Input](#input))
- Assembly FASTA files + a metadata CSV with Sample_ID and grouping column (e.g. Lineage, GPSC)

### Steps

1. **Prepare your manifest** (manifest.tsv):

| species         | metadata     | assemblies           | target_groups     | atb_exclude_species |
| --------------- | ------------ | -------------------- | ----------------- | ------------------- |
| vibrio_cholerae | metadata.csv | /path/to/assemblies/ | Lineage1,Lineage2 |                     |

2. **Prepare your metadata table** (metadata.csv):

| Sample_ID | Lineage  | other_columns |
| --------- | -------- | ------------- |
| sample_1  | Lineage1 | data...       |
| sample_2  | Lineage1 | data...       |
| sample_3  | Lineage2 | data...       |

3. **Run the pipeline:**

   ```bash
   nextflow run main.nf \
     -profile sanger_local \
     --manifest manifest.tsv \
     --assembly_suffix .fasta.gz \
     --group_label Lineage \
     --sample_col Sample_ID \
     --outdir my_output
   ```

4. **Inspect results** in `my_output/` (see [Output](#output) for file descriptions).

### On the Sanger farm (LSF)

Load the module and use `bsub.py` to submit:

```bash
module load lsmd bsub.py

bsub.py --threads 32 64 lsmd_run \
  "lsmd -profile sanger_local \
    --manifest manifest.tsv \
    --assembly_suffix .fasta.gz \
    --group_label Lineage \
    --sample_col Sample_ID \
    --outdir my_output"
```

The `--threads 32 64` request (32 CPUs, 64 GB RAM) matches the pipeline's resource requirements.

### Optional: Marker post-processing and primer design

By default, the pipeline stops at the ATB-checked markers FASTA. `--marker_post_processing` filters markers by length/GC% and soft-masks them for downstream assay design (primer3 PCR primers and/or bait-capture tiling); `--primer3_design` then designs primers with primer3. `--primer3_design` requires `--marker_post_processing`; the run stops at launch if it's set on its own.

```bash
nextflow run main.nf \
  -profile sanger_local \
  --manifest manifest.tsv \
  --group_label Lineage \
  --sample_col Sample_ID \
  --outdir my_output \
  --marker_post_processing \
  --primer3_design \
  --marker_min_length 100 \
  --marker_gc_min 35.0 \
  --marker_gc_max 60.0
```

Filtered/soft-masked markers go to `my_output/post_processed_markers/`; designed primers to `my_output/<group>/primers/` (one `_primers.tsv` per group, plus `_no_primers.tsv` listing markers primer3 could not design against).

## Input

### Manifest TSV (--manifest, required)

One row per species. Columns:

| Column              | Required | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| ------------------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| species             | yes      | The species' [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria) (ATB) colour name: the full binomial, lowercase, with an underscore, e.g. `streptococcus_pneumoniae` or `vibrio_cholerae`. Used as the output-file prefix (so no whitespace or `/`) and as the target species of the ATB cross-species check. ATB splits some species into lettered chunks (_S. pneumoniae_ into `a`/`b`/`c`/...); give the base name and the chunks are picked up automatically. Check the name with `grep -i '<species>' <atb_colour_names>`. **Not in ATB:** the run still goes ahead, but that species skips the ATB check (its markers go forward unverified) and the log warns at launch, suggesting the closest ATB names in case it's a typo |
| metadata            | yes      | Path to a CSV file with Sample_ID + grouping column (--group_label); see below                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| assemblies          | yes      | Directory of assembly FASTAs, or a .txt file listing one assembly path per line                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| target_groups       | no       | Comma-separated `--group_label` values to discover markers for (e.g. `GPSC1,GPSC2`). Blank = every group with ≥ `--candidate_min_genome_count` genomes. See the note below.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| atb_exclude_species | no       | Comma-separated ATB colour name(s) left out of the ATB check's "absent from every other species" (`atb_max_outside`) test for this species, e.g. close relatives ATB can't reliably tell apart from your target. `unknown` (ATB's catch-all bucket for unassigned/low-confidence genomes) is **always** excluded; anything listed here is added on top                                                                                                                                                                                                                                                                                                                                                                                               |

The manifest must be tab-separated with exactly these five header columns (any order). A `.csv` file, a comma-separated header, or a missing or unrecognised column stops the run before any jobs start.

**Example:**

| species                  | metadata     | assemblies          | target_groups     | atb_exclude_species |
| ------------------------ | ------------ | ------------------- | ----------------- | ------------------- |
| streptococcus_pneumoniae | metadata.csv | /data/s_pneumoniae/ | GPSC1,GPSC2,GPSC3 |                     |
| vibrio_cholerae          | metadata.csv | /data/v_cholerae/   | 7PET              |                     |

Each row is processed independently; one run can build indexes for multiple species.

**About `target_groups`:**

- **Specify groups explicitly** (e.g. `GPSC1,GPSC2,GPSC3`): Candidate markers are discovered only for the listed groups. Use the exact values from your `--group_label` column in the metadata, after [label cleaning](#cleaning-group-labels)
- **Leave blank**: Candidate markers are discovered for **every** group with ≥ `--candidate_min_genome_count` genomes, after label cleaning. Only the exact label `unclassified` is skipped. Useful for exploring new datasets and discovering markers for all major groups automatically
- **Catch-all labels are real groups, not targets you want.** A label like `Non-7PET_unclassified` is not `unclassified`: it becomes its own group, and would become a target if `target_groups` were blank. For _V. cholerae_, **7PET is the only target**, so list `7PET` explicitly. `Non-7PET_unclassified` still does useful work as a non-target: it counts as an outside group, so markers that also appear in non-7PET genomes are rejected. It isn't dropped like `unclassified` genomes are
- For species with many groups (e.g. _S. pneumoniae_ with 765+ GPSCs), specifying groups explicitly is more efficient than discovering markers for all of them

**Note on folder naming:** `<group>` in the output paths is a placeholder for each value in your `--group_label` column. With `--group_label GPSC` you get `candidate_marker_filtering/streptococcus_pneumoniae_GPSC1_…`, `atb_cross_species/GPSC1/`, …; with `--group_label Lineage` and species `vibrio_cholerae` you get `atb_cross_species/7PET/`.

### Metadata table (per-species)

A CSV (comma-separated) file with one row per genome assembly. Required columns:

| Column                              | Meaning                                                                                                                          |
| ----------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| --sample_col (default: `Sample_ID`) | Genome identifier, matched against assembly filenames after normalisation (e.g. metadata `VC_O1_8` matches file `VC_O1_8.fasta`) |
| --group_label (required)            | The column whose values become the groups (e.g. `GPSC`, `Lineage`, `sublineage`)                                                 |

**Example:**

| Sample_ID | GPSC  | Country | Resistance |
| --------- | ----- | ------- | ---------- |
| sample_1  | GPSC1 | UK      | sensitive  |
| sample_2  | GPSC1 | USA     | resistant  |
| sample_3  | GPSC2 | UK      | sensitive  |
| sample_4  | NA    | UK      | sensitive  |

Sample 4's GPSC is a missing value, so it's labelled `unclassified` and left out of the index (see [Unclassified genomes](#unclassified-genomes)).

### Cleaning group labels

Group labels become colour groups. Messy labels create fake groups: `1215;5` would become its own group instead of GPSC5, and would then reject GPSC5's own markers as an "outside group". Labels are cleaned **before** the index is built, so every later step (targets, genome-count thresholds, specificity filtering) uses the cleaned labels. There's nothing to configure: the rules are fixed.

Metadata is read as plain text: `3` stays `3`, not `3.0`. Headers, sample IDs and labels are trimmed of surrounding spaces.

Each label goes through these rules **in order**:

| Order | Rule            | What it does                                                                                                                                                                                                                                                                                                                                     |
| ----- | --------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| 1     | Missing values  | Blank labels and the missing values below become `unclassified`                                                                                                                                                                                                                                                                                  |
| 2     | GPSC `;` labels | Only when `--group_label` is `GPSC` (any case). A merge-history label becomes its smallest number, keeping its `GPSC` prefix if it has one: `1215;5` → `5`, `GPSC3;28` → `GPSC3`. **Exception:** 235 with 9 in any order or prefix form (`235;9`, `9;235`, `GPSC235;9`, `GPSC9;235`) becomes its own group `235_9` (`GPSC235_9` with the prefix) |

For any other `--group_label`, labels containing `;` are left as written.

**Missing values** (matched case-insensitively):

`""` `NA` `N/A` `#N/A` `NaN` `null` `none` `unknown` `missing` `-` `?` `.` `not applicable` `not available` `not collected` `not provided`

**The colour-mapping step stops if** `--group_label` is `GPSC` and a `;` label has a part that isn't a whole number (e.g. `5;abc`). All bad labels are listed; fix them in the metadata.

**Checking the result:** every changed label is listed in `colour_mapping/<species>_stats.json` under `label_changes` (raw label, new label, genome count, and which rule changed it: `missing_value` or `gpsc_multi`), and printed in the colour-mapping log.

### Unclassified genomes

A genome ends up `unclassified` when:

- its label is blank or a missing value
- its metadata label already reads `unclassified` (any case)
- its assembly file has no metadata row (Sample_ID taken from the filename)

Unclassified genomes are **always left out of the index**. They're listed in `colour_mapping/<species>_dropped_unclassified.tsv` with the reason (`missing_value`, `labelled_unclassified` or `no_metadata_row`). The colour-mapping step stops if no genomes are left.

**Markers are not checked against unclassified genomes.** There's no evidence that markers are absent from them. For _S. pneumoniae_ this includes GPS novel clusters (`NA`), which are real non-targets. Classify unknown genomes where possible, and validate markers against near-neighbours (e.g. with BLAST) as a backstop.

**Also left out:** metadata rows whose assembly FASTA can't be found, and non-existent paths in a `.txt` assembly list.

### Worked examples

**_S. pneumoniae_ (GPSC).** Per the [GPS notes](https://www.pneumogen.net/gps/assigningGPSCs.html), a label like `1215;5` is a **merge history**: the canonical GPSC is the smaller number. `235;9` is the exception, a mixture rather than a merge. `NA` is a **novel cluster** with no GPSC yet. With `--group_label GPSC`:

- `1215;5` → `5` and `1250;156` → `156`
- `235;9` → `235_9`, its own group (it isn't in the current GPS metadata, but is handled if it appears)
- `NA` → `unclassified`, dropped from the index, so markers aren't checked against novel clusters

**_V. cholerae_ (7PET).** 7PET is the only target, so set `target_groups` to `7PET`. Every other lineage, including `Non-7PET_unclassified`, is an outside group that markers are checked against. Genomes with no metadata row are of unknown lineage and may be 7PET, so they're dropped rather than rejecting 7PET markers.

### Run-wide parameters

These apply to all species in the manifest:

- `--sample_col` (default: `Sample_ID`): metadata column matched to filenames
- `--group_label` (required): metadata column whose values become the groups
- `--assembly_suffix` (default: `.fasta`): suffix appended to `--sample_col` to form the expected filename (directory input only)

## Output

Results are written to `--outdir` (default: `./results`), one set per species.

### Directory structure

```
results/
├── colour_mapping/
│   ├── <species>_file_colours_input.txt         # Themisto input
│   ├── <species>_label_mapping.tsv             # Sample_ID → group, ordered by colour ID
│   ├── <species>_stats.json                    # Genome counts + label-cleaning summary
│   └── <species>_dropped_unclassified.tsv      # Unclassified genomes left out of the index
├── themisto2/
│   ├── <species>_build/                        # Species-wide Themisto2 index
│   ├── <species>_export/                       # Exported unitigs & colour sets
│   └── candidate_<group>_build/                # Per-group candidate index (QC gate only)
├── candidate_marker_filtering/
│   ├── <species>_<group>_candidate_unitigs.fasta       # Group-core, group-specific candidates
│   ├── <species>_<group>_specificity.tsv               # Unitig-level specificity scores
│   └── <species>_<group>_stats.txt                     # Group-level filtering stats
├── atb_cross_species/<group>/
│   ├── <group>_atb_check_PASS.fasta            # Final markers (passed the ATB check)
│   ├── <group>_atb_check_FLAG.fasta            # Also found in another ATB species (dropped)
│   ├── <group>_atb_check_ABSENT.fasta          # Not found in the target species at all (dropped)
│   ├── <group>_atb_check_validation.tsv        # Per-marker fractions, verdict and reason
│   ├── <group>_atb_check_summary.txt           # Counts, pass rate, top off-target species
│   └── <group>_atb_pseudoalign.jsonl           # Raw pseudoalignment against ATB
├── post_processed_markers/                     # With --marker_post_processing
│   ├── <species>_<group>_markers.fasta         # Soft-masked (length/GC filtered)
│   ├── <species>_<group>_rejected_markers.fasta        # Rejected candidates
│   └── <species>_<group>_marker_analysis.png           # Length/GC diagnostic plot
├── <group>/primers/                            # With --primer3_design
└── ggcat/, sbwt/, checkpoints/                 # With --publish_intermediate only
```

### Key output files

Paths are relative to `--outdir`.

| File                                                                   | Description                                                                                                                                         |
| ---------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `colour_mapping/<species>_stats.json`                                  | Genomes written and dropped, groups, and every label change                                                                                         |
| `colour_mapping/<species>_dropped_unclassified.tsv`                    | Unclassified genomes left out of the index, with the reason                                                                                         |
| `candidate_marker_filtering/<species>_<group>_candidate_unitigs.fasta` | Candidate markers after group-specificity filtering, before the ATB check                                                                           |
| `candidate_marker_filtering/<species>_<group>_specificity.tsv`         | Per-unitig within-group / max-outside-group presence (diagnostic)                                                                                   |
| `atb_cross_species/<group>/<group>_atb_check_PASS.fasta`               | **Final markers**: candidates that passed the ATB cross-species check. For a species that isn't in ATB, the candidates go forward unchecked instead |
| `atb_cross_species/<group>/<group>_atb_check_validation.tsv`           | Why each candidate passed or failed the ATB check                                                                                                   |
| `post_processed_markers/<species>_<group>_markers.fasta`               | `--marker_post_processing`: markers filtered by length/GC and soft-masked (a subset of the final markers)                                           |
| `<group>/primers/<species>_<group>_primers.tsv`                        | `--primer3_design`: designed primer pairs                                                                                                           |
| `checkpoints/pipeline_counts.tsv`                                      | `--publish_intermediate`: marker counts at each stage                                                                                               |

### stats.json fields

`colour_mapping/<species>_stats.json` summarises how metadata and assemblies were matched up, then how labels were cleaned:

| Field                                | Meaning                                                                                                              |
| ------------------------------------ | -------------------------------------------------------------------------------------------------------------------- |
| `species`                            | Species name (the manifest `species` value)                                                                          |
| `group_label_column`                 | The `--group_label` column used                                                                                      |
| `assemblies_total`                   | `assemblies_written` + `assemblies_dropped_fasta_not_found` + `assemblies_dropped_unclassified`                      |
| `assemblies_written`                 | Genomes included in the index (one colour per line in the colour file)                                               |
| `assemblies_dropped_fasta_not_found` | Metadata rows whose expected assembly FASTA wasn't found on disk                                                     |
| `assemblies_dropped_unclassified`    | Unclassified genomes left out of the index (listed in `_dropped_unclassified.tsv`)                                   |
| `assembly_paths_missing_file`        | Entries in a `.txt` assembly path-list that don't exist on disk                                                      |
| `relabelled_unclassified`            | Metadata rows whose final label is `unclassified` because it was blank or a missing value                            |
| `assemblies_without_metadata_row`    | Assembly files with no metadata row (labelled `unclassified`)                                                        |
| `assemblies_per_group`               | Genome count per group after cleaning, including `unclassified` when kept                                            |
| `missing_values`                     | The labels treated as missing (case-insensitive)                                                                     |
| `label_changes`                      | One entry per changed label: `raw_label`, `new_label`, `genomes`, and `changed_by` (`missing_value` or `gpsc_multi`) |

## How it works

### Pipeline stages

1. **BUILD_COLOUR_INDEX**: clean group labels (see [Cleaning group labels](#cleaning-group-labels)), colour-map assemblies by group, and build a species-wide SBWT/Themisto2 index from all genomes. Changing any label column rebuilds this index, which takes a long time for large species (~42k genomes for GPSC), so set them before a run
2. **Group-specificity filtering**: for each target group (from `target_groups`, or every group with ≥ `--candidate_min_genome_count` genomes when it's left blank), keep only k-mers that are group-core (present in ≥ `--candidate_min_freq` of the group) and group-specific (present in ≤ `--specificity_max_outside` of any single other group)
3. **Rebuild candidate index**: rebuild each group's filtered k-mers into its own SBWT/Themisto2 index (a QC gate), then dump it back to candidate unitigs
4. **ATB cross-species check**: pseudoalign the candidates against the AllTheBacteria species index (`--atb_index`) and keep only markers found in ≥ `--atb_min_within` of the target species' k-mers and ≤ `--atb_max_outside` of any other ATB species → final markers. Species that aren't in ATB skip this step with a warning, and their candidates go forward unchecked
5. **POST_PROCESS_MARKERS** (opt-in, `--marker_post_processing`): filter on length and global GC%; soft-mask (lowercase) local windows outside the GC range
6. **DESIGN_PRIMERS** (opt-in, `--primer3_design`): run primer3 on the soft-masked markers

### Key concepts

- **species index**: all genomes of a species, coloured by group (`--group_label`)
- **candidate index**: k-mers that are group-core and group-specific within the species (per target group)
- **markers**: candidates that also pass the ATB cross-species check, i.e. aren't found in other bacterial species

### Why not a colour index per group?

Building a full colour index for every group doesn't scale (e.g. _S. pneumoniae_ has 765+ GPSCs). Instead, LSMD filters the species-wide Themisto2 export and rebuilds only the survivors into each group's candidate index — much faster and more memory-efficient.

## Parameters

Run `nextflow run main.nf --help` for the full, always-up-to-date list.

### Core parameters

| Option                   | Type    | Default     | Description                                                                                                                           |
| ------------------------ | ------- | ----------- | ------------------------------------------------------------------------------------------------------------------------------------- |
| `--manifest`             | path    | —           | **Required.** Manifest TSV (see [Input](#input))                                                                                      |
| `--group_label`          | string  | —           | **Required.** Metadata column whose values become the groups                                                                          |
| `--sample_col`           | string  | `Sample_ID` | Metadata column matched to assembly filenames                                                                                         |
| `--assembly_suffix`      | string  | `.fasta`    | Suffix appended to `--sample_col` to form the expected filename (directory input only)                                                |
| `--outdir`               | path    | `./results` | Output directory                                                                                                                      |
| `--publish_intermediate` | boolean | `false`     | Also publish GGCAT unitigs, raw SBWT builds, dumped candidate unitigs and per-stage marker counts (`ggcat/`, `sbwt/`, `checkpoints/`) |

### Index building

| Option                     | Type    | Default | Description                                                                          |
| -------------------------- | ------- | ------- | ------------------------------------------------------------------------------------ |
| `--colour_index_kmer_size` | integer | 31      | k-mer size for GGCAT/SBWT/Themisto2                                                  |
| `--temp_dir`               | path    | —       | Scratch root for GGCAT/SBWT temp files (optional; falls back to task-local work dir) |
| `--temp_space`             | integer | 10000   | Temp storage (MB) requested for processes needing it                                 |

### Group-specificity filtering

| Option                         | Type    | Default | Description                                                                                                                            |
| ------------------------------ | ------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `--candidate_min_freq`         | string  | `core`  | Within-group presence cutoff: `core` (≥0.95), `relaxed` (≥0.5), `catchall` (>0), or a literal 0.0–1.0 (e.g. `0.8`)                     |
| `--candidate_min_genome_count` | integer | 5       | Absolute genome-count floor (on top of `--candidate_min_freq`). Also the size below which another group is ignored as an outside group |
| `--specificity_max_outside`    | string  | `0.05`  | Max presence fraction (0.0–1.0) allowed in any single other group; `null` = core-only (no specificity test)                            |

### ATB cross-species check

The defaults point at the ATB species index on the Sanger farm. Off the farm, supply your own.

| Option               | Type  | Default                        | Description                                                            |
| -------------------- | ----- | ------------------------------ | ---------------------------------------------------------------------- |
| `--atb_index`        | path  | ATB-species.thm2 (Sanger farm) | Themisto2 index of AllTheBacteria, one colour per species              |
| `--atb_colour_names` | path  | color_names.txt (Sanger farm)  | Colour ID → ATB species name, matching `--atb_index`                   |
| `--atb_min_within`   | float | 0.95                           | Minimum fraction of a marker's k-mers found in the target species      |
| `--atb_max_outside`  | float | `--specificity_max_outside`    | Maximum fraction of a marker's k-mers allowed in any other ATB species |

### Marker post-processing (--marker_post_processing)

| Option                     | Type    | Default                    | Description                                                                                |
| -------------------------- | ------- | -------------------------- | ------------------------------------------------------------------------------------------ |
| `--marker_post_processing` | boolean | `false`                    | Filter/mask markers for downstream assay design (primer3 and/or bait capture)              |
| `--marker_min_length`      | integer | 100                        | Minimum marker length (bp)                                                                 |
| `--marker_gc_min`          | float   | 35.0                       | Minimum global GC%                                                                         |
| `--marker_gc_max`          | float   | 60.0                       | Maximum global GC%                                                                         |
| `--marker_window_size`     | integer | `--colour_index_kmer_size` | Sliding-window size (bp) for local GC check; out-of-range windows soft-masked (lowercased) |
| `--marker_write_rejected`  | boolean | `true`                     | Write rejected candidates to separate FASTA                                                |
| `--marker_plot`            | boolean | `true`                     | Generate length/GC diagnostic plot                                                         |

### Primer design (--primer3_design)

Requires `--marker_post_processing`; the run stops at launch if it's set on its own.

| Option             | Type    | Default | Description                             |
| ------------------ | ------- | ------- | --------------------------------------- |
| `--primer3_design` | boolean | `false` | Run primer3_core on soft-masked markers |

## Troubleshooting

### Q: My group has no markers (empty FASTA)

**A:** Either no k-mers passed group-specificity filtering, or none passed the ATB check.

Try one of:

- **Check where they were lost:** `candidate_marker_filtering/<species>_<group>_stats.txt` for group-specificity filtering, and `atb_cross_species/<group>/<group>_atb_check_summary.txt` for the ATB check
- **Relax the within-group cutoff:** `--candidate_min_freq relaxed` (≥0.5 instead of ≥0.95)
- **Increase other-group tolerance:** `--specificity_max_outside 0.1` (allow up to 10% presence in other groups)
- **If most candidates are `FLAG` in the ATB check:** look at the top off-target species in the summary. A close relative ATB can't tell apart from your species can be listed in `atb_exclude_species`

### Q: A group I expected is missing, or there's a group I didn't expect

**A:** Check `label_changes` and `assemblies_per_group` in `colour_mapping/<species>_stats.json`. A GPSC label may have been merged (`gpsc_multi`), sent to `unclassified` and dropped (`missing_value`), or kept as written because no rule applied. Fix the label in the metadata.

### Q: A big group is losing markers it should have

**A:** Look at the `outside_lineage` column in `<species>_<group>_specificity.tsv`. If it's a catch-all label (e.g. `Non-7PET_unclassified`) that may contain your target's genomes, classify those genomes. If it's a GPS merge label like `1215;5`, check `--group_label` is `GPSC` so it's merged into its GPSC.

### Q: The ATB check step is queued for a long time

**A:** Each group's ATB pseudoalignment loads the full ATB species index, so it asks for 350 GB (doubling on each retry); under the `standard` profile on LSF that goes to the `hugemem` queue. With many target groups, those jobs queue behind each other. List fewer groups in `target_groups` to reduce the load.

### Q: Can I reuse indexes from a previous run?

**A:** Yes, with `-resume` from the same launch directory and work directory. The species-wide index is reused as long as the metadata and assemblies are unchanged, so changing `--candidate_min_freq`, `target_groups` or the ATB settings only reruns the later steps.

## Known limitations and follow-ups

- **Group labels with `/`, spaces or `;`** end up in output paths and commands and may break them. They aren't checked yet; rename such labels in the metadata (GPSC `;` labels are handled, see [Cleaning group labels](#cleaning-group-labels))
- **Markers aren't checked against unclassified genomes**, since they're always dropped. For GPS this includes the `NA` novel clusters, which are real non-targets. Validate markers against near-neighbours (e.g. with BLAST) as a backstop. A classification step for unknown genomes is planned before the pipeline runs; if popPUNK's `_clusters.csv` is available, novel clusters could be given their real clusters and kept

## Software versions

| Software           | Version                     | Container                                               | Used by                                                                     |
| ------------------ | --------------------------- | ------------------------------------------------------- | --------------------------------------------------------------------------- |
| GGCAT              | 2.2.0                       | `quay.io/biocontainers/ggcat:2.2.0--hf1b6044_0`         | index building                                                              |
| SBWT (sbwt-rs-cli) | 0.4.2 (patched, `-f93d92c`) | Sanger-internal `.sif`                                  | index building, candidate unitig dump                                       |
| Themisto2          | 0.0.1                       | `quay.io/sangerpathogens/themisto2:0.0.1`               | index building, ATB pseudoalignment                                         |
| pandas             | 2.2.1                       | `quay.io/sangerpathogens/pandas:2.2.1`                  | colour mapping, group-specificity filter, ATB check, marker post-processing |
| seqkit             | 2.10.0                      | `quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0`       | per-stage count checkpoints                                                 |
| primer3            | 2.6.1                       | `quay.io/biocontainers/primer3:2.6.1--pl5321h503566f_7` | primer design                                                               |

All software dependencies are containerised (Docker/Singularity).

## Issues and contributions

- **GitHub/GitLab:** Log an issue or open a merge request on this repository
- **Sanger users:** Raise an issue on the PAM Freshservice portal: https://sanger.freshservice.com/support/catalog/items/426
