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

**Run the test dataset:**

```bash
git clone --recurse-submodules git@gitlab.internal.sanger.ac.uk:sanger-pathogens/pipelines/lsmd.git
cd lsmd

nextflow run main.nf -profile test,sanger_local --outdir test_output
```

**Inspect results** in `test_output/` (11-genome _V. tarriae_ dataset). Once satisfied, clean up:

```bash
rm -rf work .nextflow*
```

## Running your own data

### Prerequisites

- Nextflow ≥ 21.04.0
- Docker or Singularity
- A manifest TSV (one row per species; see [Input](#input))
- Assembly FASTA files + a metadata CSV with Sample_ID and grouping column (e.g. Lineage, GPSC)

### Steps

1. **Prepare your manifest** (manifest.tsv):

| species   | metadata     | assemblies           | target_groups     |
| --------- | ------------ | -------------------- | ----------------- |
| v_tarriae | metadata.csv | /path/to/assemblies/ | Lineage1,Lineage2 |

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

| Column               | Required | Meaning                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      |
| -------------------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| species              | yes      | Short name for the species (no whitespace or `/`); used as the output-file prefix                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            |
| metadata             | yes      | Path to a CSV file with Sample_ID + grouping column (--group_label); see below                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| assemblies           | yes      | Directory of assembly FASTAs, or a .txt file listing one assembly path per line                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              |
| target_groups        | no       | Comma-separated `--group_label` values to discover markers for (e.g. `GPSC1,GPSC2`). Blank/absent = every group with ≥ `--candidate_min_genome_count` genomes. See the note below.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           |
| atb_target_species   | no       | Space-separated ATB colour name(s) for the ATB cross-species check (candidate markers get pseudoaligned against [AllTheBacteria](https://github.com/AllTheBacteria/AllTheBacteria) to confirm they're specific to your species). **Write it as the full binomial with an underscore, lowercase** — e.g. `vibrio_cholerae`, not `v_cholerae` — because ATB's colour names don't match the short `species` key you use elsewhere in the manifest. Look the name(s) up in ATB's `color_names.txt` if unsure (`grep -i '<species>' <atb_color_names>`). Some species span more than one ATB colour (e.g. _S. pneumoniae_ splits into lettered chunks `a`/`b`/`c`/...) — list the base name(s) here and the pipeline picks up clean lettered splits automatically. **Leave blank** if your species isn't in ATB at all: candidate markers then pass through this check unverified, with a loud warning in the log, rather than the run failing or markers being silently dropped. |
| atb_exclude_species  | no       | Space-separated ATB colour name(s) left out of the ATB check's "absent from every other species" (`atb_max_outside`) test for this species, e.g. close relatives ATB can't reliably tell apart from your target. Same naming as `atb_target_species`. **Blank** = `unknown` (ATB's catch-all bucket for unassigned/low-confidence genomes). If you fill it in, add `unknown` yourself to keep excluding that bucket.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         |
| label_map            | no       | Absolute path to a TSV of exact label fixes (`raw_label` → `group`). Applied first and final. **Blank** = no map. See [Cleaning group labels](#cleaning-group-labels).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| label_missing        | no       | `\|`-separated labels that mean "no value", e.g. `NA\|unknown\|not applicable`. Case-insensitive. **Blank** = the built-in list. Filling it in **replaces** the list.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        |
| label_multi          | no       | What to do with labels containing `;`: `keep`, `smallest` or `unclassified`. **Blank** = `keep`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             |
| unclassified_genomes | no       | `keep` or `drop` genomes that end up `unclassified`. **Blank** = `keep`. See [Unclassified genomes](#unclassified-genomes).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  |

**Example:**

| species    | metadata     | assemblies          | target_groups     | atb_target_species       | atb_exclude_species | label_map | label_missing | label_multi | unclassified_genomes |
| ---------- | ------------ | ------------------- | ----------------- | ------------------------ | ------------------- | --------- | ------------- | ----------- | -------------------- |
| s_pneu     | metadata.csv | /data/s_pneumoniae/ | GPSC1,GPSC2,GPSC3 | streptococcus_pneumoniae |                     |           |               | smallest    | keep                 |
| v_cholerae | metadata.csv | /data/v_cholerae/   | 7PET              | vibrio_cholerae          | unknown             |           |               |             | drop                 |

Each row is processed independently; one run can build indexes for multiple species.

**About `target_groups`:**

- **Specify groups explicitly** (e.g. `GPSC1,GPSC2,GPSC3`): Candidate markers are discovered only for the listed groups. Use the exact values from your `--group_label` column in the metadata, after [label cleaning](#cleaning-group-labels)
- **Leave blank**: Candidate markers are discovered for **every** group with ≥ `--candidate_min_genome_count` genomes, after label cleaning. Only the exact label `unclassified` is skipped. Useful for exploring new datasets and discovering markers for all major groups automatically
- **Catch-all labels are real groups, not targets you want.** A label like `Non-7PET_unclassified` is not `unclassified`: it becomes its own group, and would become a target if `target_groups` were blank. For _V. cholerae_, **7PET is the only target**, so list `7PET` explicitly. `Non-7PET_unclassified` still does useful work as a non-target: it counts as an outside group, so markers that also appear in non-7PET genomes are rejected. `unclassified_genomes` doesn't affect it
- For species with many groups (e.g. _S. pneumoniae_ with 765+ GPSCs), specifying groups explicitly is more efficient than discovering markers for all of them

**Note on folder naming:** `<group>` in the output paths is a placeholder for each value in your `--group_label` column. With `--group_label GPSC` you get `candidate_marker_filtering/s_pneu_GPSC1_…`, `atb_cross_species/GPSC1/`, …; with `--group_label Lineage` and species `v_cholerae` you get `atb_cross_species/7PET/`.

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

Sample 4's GPSC is a missing value, so it's labelled `unclassified` (see [Unclassified genomes](#unclassified-genomes)).

### Cleaning group labels

Group labels become colour groups. Messy labels create fake groups: `1215;5` would become its own group instead of GPSC5, and would then reject GPSC5's own markers as an "outside group". Three optional manifest columns clean labels **before** the index is built, so every later step (targets, genome-count thresholds, specificity filtering) uses the cleaned labels.

Metadata is read as plain text: `3` stays `3`, not `3.0`. Headers, sample IDs and labels are trimmed of surrounding spaces.

Each label goes through these rules **in order**. The first match decides:

| Order | Column          | What it does                                                                                                               |
| ----- | --------------- | -------------------------------------------------------------------------------------------------------------------------- |
| 1     | `label_map`     | Exact `raw_label` → `group` override. Final: no later rule touches it. Map to `unclassified` to send a label to background |
| 2     | `label_missing` | Labels that mean "no value" become `unclassified`                                                                          |
| 3     | `label_multi`   | Labels containing `;`: `keep` as written, `smallest` whole number (`1215;5` → `5`), or send to `unclassified`              |

Leave all three blank and labels are used as written, apart from the missing-value list.

**`label_map` file** (tab-separated, header required):

```tsv
raw_label	group
235;9	235_9
```

- Use an **absolute path**. Relative paths resolve against the launch directory, not the manifest
- Editing the map triggers a rebuild of that species' index, even with `-resume`
- A map entry's `group` can't be a missing value such as `NA`. Use `unclassified` instead
- Also useful for renaming labels that would break output folder names (e.g. `;`, `/` or spaces), since groups appear in output paths

**Built-in missing values** (used when `label_missing` is blank; matched case-insensitively):

`""` `NA` `N/A` `#N/A` `NaN` `null` `none` `unknown` `missing` `-` `?` `.` `not applicable` `not available` `not collected` `not provided`

If you fill in `label_missing`, it **replaces** this list. Empty labels always count as missing.

**The run stops, before any jobs start, if:**

- `label_multi` or `unclassified_genomes` has an unrecognised value
- the `label_map` path doesn't exist

**The colour-mapping step stops if:**

- the map lacks `raw_label` or `group` columns, lists a `raw_label` twice, has a blank `raw_label` or `group`, or maps to a missing value
- `smallest` meets a `;` label that isn't all whole numbers (e.g. `5;abc`). All bad labels are listed; fix them with `label_map`

**Checking the result:** every changed label is listed in `color_mapping/<species>_stats.json` under `label_changes` (raw label, new label, genome count, and which rule changed it), and printed in the colour-mapping log. Map entries that matched no label are listed under `label_map_unmatched`.

### Unclassified genomes

A genome ends up `unclassified` when:

- its label is blank or a missing value
- `label_map` or `label_multi` sent it there
- its metadata label literally reads `unclassified`
- its assembly file has no metadata row (Sample_ID taken from the filename)

`unclassified_genomes` decides what happens to them:

| Value            | Effect                                                                                                                                                                                                                                     | Use when                                                             |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------- |
| `keep` (default) | Kept in the index as one `unclassified` group. **Never a target, but it does count as an outside group** (if ≥ `--candidate_min_genome_count` genomes), so it can reject markers                                                           | They're known non-targets, e.g. GPS novel clusters                   |
| `drop`           | Removed from the index. Listed in `color_mapping/<species>_dropped_unclassified.tsv` with the reason (`label_missing`, `label_map`, `label_multi`, `labelled_unclassified` or `no_metadata_row`). **Markers are not checked against them** | Some may really be your target, so they'd wrongly reject its markers |

With `drop`, there's no evidence that markers are absent from the dropped genomes. Classify unknown genomes where possible, and use marker validation (e.g. BLAST against near-neighbours) as a backstop.

**Only genomes left out entirely, whatever the setting:** metadata rows whose assembly FASTA can't be found, and non-existent paths in a `.txt` assembly list.

### Worked examples

**_S. pneumoniae_ (GPSC).** Per the [GPS notes](https://www.pneumogen.net/gps/assigningGPSCs.html), a label like `1215;5` is a **merge history**: the canonical GPSC is the smaller number. `235;9` is the exception, a mixture rather than a merge. `NA` is a **novel cluster** with no GPSC yet.

| Column                 | Value                       | Why                                                                                                         |
| ---------------------- | --------------------------- | ----------------------------------------------------------------------------------------------------------- |
| `label_multi`          | `smallest`                  | `1215;5` → `5`, `1250;156` → `156`                                                                          |
| `label_map`            | optional: `235;9` → `235_9` | Safeguard. `235;9` isn't in the current GPS metadata, but `smallest` would fold it into GPSC9 if it appears |
| `label_missing`        | blank                       | `NA` is in the built-in list                                                                                |
| `unclassified_genomes` | `keep`                      | Novel clusters are real non-targets, so they should still reject cross-reacting markers                     |

Result: 325 targets with `target_groups` blank, and GPSC5 grows from 1107 to 1135 genomes. The 74 `NA` genomes form one pooled outside group (see [Known limitations](#known-limitations-and-follow-ups)).

**_V. cholerae_ (7PET).** 7PET is the only target.

| Column                 | Value  | Why                                                                                                                                  |
| ---------------------- | ------ | ------------------------------------------------------------------------------------------------------------------------------------ |
| `target_groups`        | `7PET` | Only 7PET gets markers. Every other lineage, including `Non-7PET_unclassified`, is an outside group that markers are checked against |
| `unclassified_genomes` | `drop` | Genomes with no metadata row are of unknown lineage and may be 7PET, so they shouldn't reject 7PET markers                           |

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
├── color_mapping/
│   ├── <species>_file_colors_input.txt         # Themisto input
│   ├── <species>_label_mapping.tsv             # Sample_ID → group, ordered by colour ID
│   ├── <species>_stats.json                    # Genome counts + label-cleaning summary
│   └── <species>_dropped_unclassified.tsv      # With unclassified_genomes = drop
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

| File                                                                   | Description                                                                                                                                                    |
| ---------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `color_mapping/<species>_stats.json`                                   | Genomes written and dropped, groups, and every label change                                                                                                    |
| `color_mapping/<species>_dropped_unclassified.tsv`                     | With `unclassified_genomes = drop`: genomes left out of the index, with the reason                                                                             |
| `candidate_marker_filtering/<species>_<group>_candidate_unitigs.fasta` | Candidate markers after group-specificity filtering, before the ATB check                                                                                      |
| `candidate_marker_filtering/<species>_<group>_specificity.tsv`         | Per-unitig within-group / max-outside-group presence (diagnostic)                                                                                              |
| `atb_cross_species/<group>/<group>_atb_check_PASS.fasta`               | **Final markers**: candidates that passed the ATB cross-species check. For a species with no `atb_target_species`, the candidates go forward unchecked instead |
| `atb_cross_species/<group>/<group>_atb_check_validation.tsv`           | Why each candidate passed or failed the ATB check                                                                                                              |
| `post_processed_markers/<species>_<group>_markers.fasta`               | `--marker_post_processing`: markers filtered by length/GC and soft-masked (a subset of the final markers)                                                      |
| `<group>/primers/<species>_<group>_primers.tsv`                        | `--primer3_design`: designed primer pairs                                                                                                                      |
| `checkpoints/pipeline_counts.tsv`                                      | `--publish_intermediate`: marker counts at each stage                                                                                                          |

### stats.json fields

`color_mapping/<species>_stats.json` summarises how metadata and assemblies were matched up, then how labels were cleaned:

| Field                                | Meaning                                                                                                                                                            |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `species`                            | Species name (the manifest `species` value)                                                                                                                        |
| `group_label_column`                 | The `--group_label` column used                                                                                                                                    |
| `assemblies_total`                   | `assemblies_written` + `assemblies_dropped_fasta_not_found` + `assemblies_dropped_unclassified`                                                                    |
| `assemblies_written`                 | Genomes included in the index (one colour per line in the colour file)                                                                                             |
| `assemblies_dropped_fasta_not_found` | Metadata rows whose expected assembly FASTA wasn't found on disk                                                                                                   |
| `assemblies_dropped_unclassified`    | Genomes removed because `unclassified_genomes = drop` (listed in `_dropped_unclassified.tsv`)                                                                      |
| `assembly_paths_missing_file`        | Entries in a `.txt` assembly path-list that don't exist on disk                                                                                                    |
| `relabelled_unclassified`            | Metadata rows whose final label is `unclassified` because it was blank, a missing value, or sent there by `label_map` / `label_multi`                              |
| `assemblies_without_metadata_row`    | Assembly files with no metadata row (labelled `unclassified`)                                                                                                      |
| `assemblies_per_group`               | Genome count per group after cleaning, including `unclassified` when kept                                                                                          |
| `label_settings`                     | The settings used: `label_map`, `label_missing` (the values), `label_missing_source` (`built-in list` or `--label-missing`), `label_multi`, `unclassified_genomes` |
| `label_changes`                      | One entry per changed label: `raw_label`, `new_label`, `genomes`, and `changed_by` (`label_map`, `label_missing` or `label_multi`)                                 |
| `label_map_unmatched`                | Map entries that matched no label (usually a typo)                                                                                                                 |

## How it works

### Pipeline stages

1. **BUILD_COLOR_INDEX**: clean group labels (see [Cleaning group labels](#cleaning-group-labels)), colour-map assemblies by group, and build a species-wide SBWT/Themisto2 index from all genomes. Changing any label column rebuilds this index, which takes a long time for large species (~42k genomes for GPSC), so set them before a run
2. **Group-specificity filtering**: for each target group (from `target_groups`, or every group with ≥ `--candidate_min_genome_count` genomes when it's left blank), keep only k-mers that are group-core (present in ≥ `--candidate_min_freq` of the group) and group-specific (present in ≤ `--specificity_max_outside` of any single other group)
3. **Rebuild candidate index**: rebuild each group's filtered k-mers into its own SBWT/Themisto2 index (a QC gate), then dump it back to candidate unitigs
4. **ATB cross-species check**: pseudoalign the candidates against the AllTheBacteria species index (`--atb_index`) and keep only markers found in ≥ `--atb_min_within` of the target species' k-mers and ≤ `--atb_max_outside` of any other ATB species → final markers. Species with no `atb_target_species` skip this step with a warning, and their candidates go forward unchecked
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

| Option                    | Type    | Default | Description                                                                          |
| ------------------------- | ------- | ------- | ------------------------------------------------------------------------------------ |
| `--color_index_kmer_size` | integer | 31      | k-mer size for GGCAT/SBWT/Themisto2                                                  |
| `--temp_dir`              | path    | —       | Scratch root for GGCAT/SBWT temp files (optional; falls back to task-local work dir) |
| `--temp_space`            | integer | 10000   | Temp storage (MB) requested for processes needing it                                 |

### Group-specificity filtering

| Option                         | Type    | Default | Description                                                                                                                            |
| ------------------------------ | ------- | ------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `--candidate_min_freq`         | string  | `core`  | Within-group presence cutoff: `core` (≥0.95), `relaxed` (≥0.5), `catchall` (>0), or a literal 0.0–1.0 (e.g. `0.8`)                     |
| `--candidate_min_genome_count` | integer | 5       | Absolute genome-count floor (on top of `--candidate_min_freq`). Also the size below which another group is ignored as an outside group |
| `--specificity_max_outside`    | string  | `0.05`  | Max presence fraction (0.0–1.0) allowed in any single other group; `null` = core-only (no specificity test)                            |

### ATB cross-species check

The defaults point at the ATB species index on the Sanger farm. Off the farm, supply your own.

| Option              | Type  | Default                        | Description                                                            |
| ------------------- | ----- | ------------------------------ | ---------------------------------------------------------------------- |
| `--atb_index`       | path  | ATB-species.thm2 (Sanger farm) | Themisto2 index of AllTheBacteria, one colour per species              |
| `--atb_color_names` | path  | color_names.txt (Sanger farm)  | Colour ID → ATB species name, matching `--atb_index`                   |
| `--atb_min_within`  | float | 0.95                           | Minimum fraction of a marker's k-mers found in the target species      |
| `--atb_max_outside` | float | `--specificity_max_outside`    | Maximum fraction of a marker's k-mers allowed in any other ATB species |

### Marker post-processing (--marker_post_processing)

| Option                     | Type    | Default                   | Description                                                                                |
| -------------------------- | ------- | ------------------------- | ------------------------------------------------------------------------------------------ |
| `--marker_post_processing` | boolean | `false`                   | Filter/mask markers for downstream assay design (primer3 and/or bait capture)              |
| `--marker_min_length`      | integer | 100                       | Minimum marker length (bp)                                                                 |
| `--marker_gc_min`          | float   | 35.0                      | Minimum global GC%                                                                         |
| `--marker_gc_max`          | float   | 60.0                      | Maximum global GC%                                                                         |
| `--marker_window_size`     | integer | `--color_index_kmer_size` | Sliding-window size (bp) for local GC check; out-of-range windows soft-masked (lowercased) |
| `--marker_write_rejected`  | boolean | `true`                    | Write rejected candidates to separate FASTA                                                |
| `--marker_plot`            | boolean | `true`                    | Generate length/GC diagnostic plot                                                         |

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

**A:** Check `label_changes` and `assemblies_per_group` in `color_mapping/<species>_stats.json`. A label may have been merged (`label_multi`), sent to `unclassified` (missing value or map), or kept as written because nothing matched. Fix it with `label_map`.

### Q: A big group is losing markers it should have

**A:** Look at the `outside_lineage` column in `<species>_<group>_specificity.tsv`. If the group rejecting markers is `unclassified` and may contain your target's genomes, set `unclassified_genomes = drop` or classify those genomes. If it's a label like `1215;5`, set `label_multi = smallest`.

### Q: The ATB check step is queued for a long time

**A:** Each group's ATB pseudoalignment loads the full ATB species index, so it asks for 350 GB (doubling on each retry); under the `standard` profile on LSF that goes to the `hugemem` queue. With many target groups, those jobs queue behind each other. List fewer groups in `target_groups` to reduce the load.

### Q: Can I reuse indexes from a previous run?

**A:** Yes, with `-resume` from the same launch directory and work directory. The species-wide index is reused as long as the metadata, assemblies and label columns are unchanged, so changing `--candidate_min_freq`, `target_groups` or the ATB settings only reruns the later steps.

## Known limitations and follow-ups

- **Group labels with `/` or spaces** end up in output paths and may break them. Until this is checked and handled, rename such labels with `label_map`
- **GPS `NA` genomes are one pooled outside group.** They're novel clusters, kept as real non-targets, but pooling different novel clusters into one group can hide a marker that's common in one of them. If popPUNK's `_clusters.csv` is available, they could be split into their real clusters
- **`unclassified_genomes = drop` means markers aren't checked against dropped genomes.** A classification step for unknown genomes is planned before the pipeline runs

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
