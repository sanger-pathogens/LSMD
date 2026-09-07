# LSMD: Lineage-Specific Marker Discovery

## What is this?

LSMD discovers k-mer markers that distinguish a target group (lineage, sublineage, serotype, etc.) from related groups within a species. It combines the accuracy of whole-genome analysis with the speed of k-mer-based approaches — discovering markers in hours rather than weeks.

The pipeline was developed and tested for **lineage-level grouping** (e.g. GPSC for *S. pneumoniae*, 7PET for *V. cholerae*), but is flexible enough to work with any grouping you define in your metadata — serotypes, clades, resistance phenotypes, or any other categorical column. You provide the genome assemblies and define the groupings; the pipeline discovers markers that distinguish your chosen groups.

Throughout this README, **group** means whatever categorical column you point `--group_label` at (GPSC, lineage, serotype, …). "Lineage" is just the most common example.

[[_TOC_]]

## Why k-mers?

K-mer analysis enables LSMD to achieve what older discovery methods cannot:

- **Alignment-free**: No need to align genomes or call variants — faster and more objective
- **Annotation-independent**: Works with non-coding sequences; doesn't require gene calls or functional annotation
- **Unbiased**: Uses all k-mers in the input genomes, not a subset (unlike sketching/sub-sampling approaches)
- **Highly scalable**: Tested on 40,000+ genomes; computational cost scales with data volume, not complexity
- **Deployment-ready**: Output k-mers are short, concrete sequences suitable for immediate PCR primer design or bait-capture panel construction

The pipeline uses two key tools: **Themisto2** for fast, colour-mapped k-mer indexing, and **SBWT** for efficient set-difference operations.

## Quick start

**Prerequisites:** Nextflow, Docker or Singularity

**Run the test dataset:**

```bash
git clone --recurse-submodules git@gitlab.internal.sanger.ac.uk:sanger-pathogens/pipelines/lsmd.git
cd lsmd

nextflow run main.nf -profile test,sanger_local --outdir test_output
```

**Inspect results** in `test_output/` (11-genome *V. tarriae* dataset). Once satisfied, clean up:

```bash
rm -rf work .nextflow*
```

## Running your own data

### Prerequisites

- Nextflow ≥ 21.04.0
- Docker or Singularity
- A manifest TSV (one row per species; see [Input](#input))
- Assembly FASTA files + metadata table with Sample_ID and grouping column (e.g. Lineage, GPSC)

### Steps

1. **Prepare your manifest** (manifest.tsv):

| species | metadata | assemblies | target_groups |
|---------|----------|-----------|----------------|
| v_tarriae | metadata.csv | /path/to/assemblies/ | Lineage1,Lineage2 |

2. **Prepare your metadata table** (metadata.csv):

| Sample_ID | Lineage | other_columns |
|-----------|---------|----------------|
| sample_1 | Lineage1 | data... |
| sample_2 | Lineage1 | data... |
| sample_3 | Lineage2 | data... |

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

By default, the pipeline stops at the markers SBWT index and FASTA. `--marker_post_processing` filters markers by length/GC% and soft-masks them for downstream assay design (primer3 PCR primers and/or bait-capture tiling); `--primer3_design` then designs primers with primer3:

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

Filtered/soft-masked markers go to `my_output/post_processed_markers/`; designed primers to `my_output/markers_<species>_<group>/primers/` (one `_primers.tsv` per group, plus `_no_primers.tsv` listing markers primer3 could not design against).

## Input

### Manifest TSV (--manifest, required)

One row per species. Columns:

| Column | Required | Meaning |
|--------|----------|---------|
| species | yes | Short name for the species (no whitespace or `/`); used as the output-file prefix |
| metadata | yes | Path to CSV/TSV file with Sample_ID + grouping column (--group_label); see below |
| assemblies | yes | Directory of assembly FASTAs, or a .txt file listing one assembly path per line |
| target_groups | no | Comma-separated `--group_label` values to discover markers for (e.g. `GPSC1,GPSC2`). Blank/absent = every group with ≥ `--candidate_min_genome_count` genomes. See the note below. |

**Example:**

| species | metadata | assemblies | target_groups |
|---------|----------|-----------|----------------|
| s_pneu | metadata.csv | /data/s_pneumoniae/ | GPSC1,GPSC2,GPSC3 |
| v_cholerae | metadata.csv | /data/v_cholerae/ | 7PET,Non-7PET |

Each row is processed independently; one run can build indexes for multiple species.

**About `target_groups`:**
- **Specify groups explicitly** (e.g. `GPSC1,GPSC2,GPSC3`): Candidate markers are discovered only for the listed groups. Use the exact values from your `--group_label` column in the metadata
- **Leave blank**: Candidate markers are discovered for **every** group in the metadata with ≥ `--candidate_min_genome_count` genomes (`unclassified` is skipped). Useful for exploring new datasets and discovering markers for all major groups automatically
- For species with many groups (e.g. *S. pneumoniae* with 765+ GPSCs), specifying groups explicitly is more efficient than discovering markers for all of them

**Note on folder naming:** `<group>` in the output paths is a placeholder for each value in your `--group_label` column. With `--group_label GPSC` you get `markers_s_pneu_GPSC1/`, `markers_s_pneu_GPSC2/`, …; with `--group_label Lineage` and species `v_cholerae` you get `markers_v_cholerae_7PET/`, etc.

### Metadata table (per-species)

CSV or TSV with one row per genome assembly. Required columns:

| Column | Meaning |
|--------|---------|
| --sample_col (default: `Sample_ID`) | Genome identifier, matched against assembly filenames after normalisation (e.g. metadata `VC_O1_8` matches file `VC_O1_8.fasta`) |
| --group_label (required) | The column whose values become the groups (e.g. `GPSC`, `Lineage`, `sublineage`) |

**Unclassified genomes:**
- Metadata rows with a blank --group_label are kept and labelled `unclassified`
- Assembly files with no metadata row are kept and labelled `unclassified` (Sample_ID taken from filename)
- **Only excluded:** metadata rows whose assembly FASTA can't be found on disk, and non-existent paths in a `.txt` assembly list

**Example:**

| Sample_ID | GPSC | Country | Resistance |
|-----------|------|---------|------------|
| sample_1 | GPSC1 | UK | sensitive |
| sample_2 | GPSC1 | USA | resistant |
| sample_3 | GPSC2 | UK | sensitive |
| sample_4 | — | UK | sensitive |

Sample 4 has no GPSC value, so it's labelled `unclassified` in the output.

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
│   └── <species>_stats.json                    # Reconciliation summary
├── themisto2/
│   ├── <species>_build/                        # Species-wide SBWT/Themisto2 index
│   ├── <species>_export/                       # Exported unitigs & colour sets
│   └── candidate_<group>_build/                # Per-group candidate index (QC gate only)
├── sbwt/
│   ├── <species>/                              # Checked species-wide SBWT
│   ├── candidate_<group>/                      # Checked candidate SBWT
│   └── bg_excl_<species>/                      # Checked background-exclusion SBWT
├── candidate_marker_filtering/
│   ├── <species>_<group>_candidate_unitigs.fasta       # Filtered candidates pre-rebuild
│   ├── <species>_<group>_specificity.tsv               # Unitig-level specificity scores
│   └── <species>_<group>_stats.txt                     # Group-level filtering stats
├── markers_<species>_<group>/
│   ├── markers_<species>_<group>.sbwt          # Final markers index
│   ├── markers_<species>_<group>_unitigs.fasta # Markers as FASTA
│   └── primers/                                # With --primer3_design: primer outputs
├── post_processed_markers/                     # With --marker_post_processing
│   ├── <species>_<group>_markers.fasta         # Soft-masked (length/GC filtered)
│   ├── <species>_<group>_rejected_markers.fasta        # Rejected candidates
│   └── <species>_<group>_marker_analysis.png           # Length/GC diagnostic plot
└── ggcat/, checkpoints/                        # With --publish_intermediate only
```

### Key output files

Paths are relative to `--outdir`.

| File | Description |
|------|-------------|
| `color_mapping/<species>_stats.json` | Reconciliation summary: genomes written, dropped, labelled unclassified |
| `candidate_marker_filtering/<species>_<group>_candidate_unitigs.fasta` | Candidate k-mers before background subtraction (group-specificity filter output) |
| `candidate_marker_filtering/<species>_<group>_specificity.tsv` | Per-unitig within-group / max-outside-group presence (diagnostic) |
| `markers_<species>_<group>/markers_<species>_<group>_unitigs.fasta` | **Final markers** (after background subtraction), as FASTA |
| `markers_<species>_<group>/markers_<species>_<group>.sbwt` | Same final markers, as an SBWT index |
| `post_processed_markers/<species>_<group>_markers.fasta` | `--marker_post_processing`: markers filtered by length/GC and soft-masked (a subset of the final markers) |
| `markers_<species>_<group>/primers/<species>_<group>_primers.tsv` | `--primer3_design`: designed primer pairs |

### stats.json fields

`color_mapping/<species>_stats.json` summarises how metadata and assemblies reconciled:

| Field | Meaning |
|-------|---------|
| `species` | Species name (the manifest `species` value) |
| `group_label_column` | The `--group_label` column used |
| `assemblies_total` | `assemblies_written` + `assemblies_dropped_fasta_not_found` |
| `assemblies_written` | Genomes included in the index (one colour per line in the colour file) |
| `assemblies_dropped_fasta_not_found` | Metadata rows whose expected assembly FASTA wasn't found on disk |
| `assembly_paths_missing_file` | Entries in a `.txt` assembly path-list that don't exist on disk |
| `relabelled_unclassified` | Metadata rows with a blank `--group_label` (kept, labelled `unclassified`) |
| `assemblies_without_metadata_row` | Assembly files with no metadata row (kept, labelled `unclassified`) |
| `assemblies_per_group` | Genome count per group, including the `unclassified` bucket |

## How it works

### Pipeline stages

1. **BUILD_COLOR_INDEX**: colour-map assemblies by group; build a species-wide SBWT/Themisto2 index from all genomes
2. **Group-specificity filtering**: for each target group (from `target_groups`, or every group with ≥ `--candidate_min_genome_count` genomes when it's left blank), keep only k-mers that are group-core (present in ≥ `--candidate_min_freq` of the group) and group-specific (present in ≤ `--specificity_max_outside` of any single sibling group)
3. **Rebuild candidate index**: rebuild the filtered k-mers into a per-group SBWT/Themisto2 index
4. **SET_DIFF_CALCULATIONS**: subtract the background (ATB) with `sbwt difference` → final markers
5. **POST_PROCESS_MARKERS** (opt-in, `--marker_post_processing`): filter on length and global GC%; soft-mask (lowercase) local windows outside the GC range
6. **DESIGN_PRIMERS** (opt-in, `--primer3_design`): run primer3 on the soft-masked markers

### Key concepts

- **species_index**: all genomes of a species, coloured by group (`--group_label`)
- **bg_excl**: background index − species_index (k-mers in the external collection but not in the target species)
- **candidate_index**: k-mers that are group-core and group-specific (per target group)
- **markers**: candidate_index − bg_excl (final group-specific k-mers that don't appear in the background)

### Why not a colour index per group?

Building a full colour index for every group doesn't scale (e.g. *S. pneumoniae* has 765+ GPSCs). Instead, LSMD filters the species-wide Themisto2 export and rebuilds only the survivors into each group's candidate index — much faster and more memory-efficient.

## Parameters

Run `nextflow run main.nf --help` for the full, always-up-to-date list.

### Core parameters

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--manifest` | path | — | **Required.** Manifest TSV (see [Input](#input)) |
| `--group_label` | string | — | **Required.** Metadata column whose values become the groups |
| `--sample_col` | string | `Sample_ID` | Metadata column matched to assembly filenames |
| `--assembly_suffix` | string | `.fasta` | Suffix appended to `--sample_col` to form the expected filename (directory input only) |
| `--outdir` | path | `./results` | Output directory |
| `--publish_intermediate` | boolean | `false` | Also publish GGCAT unitigs, raw SBWT builds and per-stage k-mer counts (`ggcat/`, `checkpoints/`) |

### Index building

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--color_index_kmer_size` | integer | 31 | k-mer size for GGCAT/SBWT/Themisto2 (must match background index) |
| `--temp_dir` | path | — | Scratch root for GGCAT/SBWT temp files (optional; falls back to task-local work dir) |
| `--temp_space` | integer | 10000 | Temp storage (MB) requested for processes needing it |

### Group-specificity filtering

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--candidate_min_freq` | string | `core` | Within-group presence cutoff: `core` (≥0.95), `relaxed` (≥0.5), `catchall` (>0), or a literal 0.0–1.0 (e.g. `0.8`) |
| `--candidate_min_genome_count` | integer | 5 | Absolute genome-count floor (on top of `--candidate_min_freq`) |
| `--specificity_max_outside` | string | `0.05` | Max presence fraction (0.0–1.0) allowed in any single other group; `null` = core-only (no specificity test) |

### Background subtraction

Off the Sanger farm there is no default background — supply `--bg_index` (an SBWT index built at the same `--color_index_kmer_size`), or `--bg_excl_index` if you already have a `bg_excl` from a previous run.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--bg_index` | path | ATB (Sanger farm only) | Background SBWT index (must be built at same `--color_index_kmer_size`) |
| `--bg_excl_index` | path | — | Pre-computed bg_excl index to reuse (skips the large background-subtraction diff) |

### Marker post-processing (--marker_post_processing)

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--marker_post_processing` | boolean | `false` | Filter/mask markers for downstream assay design (primer3 and/or bait capture) |
| `--marker_min_length` | integer | 100 | Minimum marker length (bp) |
| `--marker_gc_min` | float | 35.0 | Minimum global GC% |
| `--marker_gc_max` | float | 60.0 | Maximum global GC% |
| `--marker_window_size` | integer | `--color_index_kmer_size` | Sliding-window size (bp) for local GC check; out-of-range windows soft-masked (lowercased) |
| `--marker_write_rejected` | boolean | `true` | Write rejected candidates to separate FASTA |
| `--marker_plot` | boolean | `true` | Generate length/GC diagnostic plot |

### Primer design (--primer3_design)

Requires `--marker_post_processing` to be enabled.

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--primer3_design` | boolean | `false` | Run primer3_core on soft-masked markers |

## Troubleshooting

### Q: My group has no markers (empty FASTA)

**A:** The group has no k-mers that pass both `--candidate_min_freq` and `--specificity_max_outside`.

Try one of:
- **Relax the within-group cutoff:** `--candidate_min_freq relaxed` (≥0.5 instead of ≥0.95)
- **Increase sibling-group tolerance:** `--specificity_max_outside 0.1` (allow up to 10% presence in siblings)
- **Check the filtering stats:** open `candidate_marker_filtering/<species>_<group>_specificity.tsv` to see where k-mers were filtered out

### Q: Background subtraction is taking forever

**A:** `sbwt difference` is memory-intensive. 

On Sanger farm:
```bash
bsub.py --threads 32 64 lsmd_run "lsmd ..."  # Request 64GB
```

If you've already computed bg_excl:
```bash
nextflow run main.nf --manifest ... --bg_excl_index <path>  # Reuse it
```

### Q: How do I use a custom background index?

**A:** Build it at the same --color_index_kmer_size (default 31), then pass it:
```bash
nextflow run main.nf --manifest ... --bg_index /path/to/custom.sbwt
```

### Q: Can I reuse indexes from a previous run?

**A:** Yes. If you've already computed bg_excl, pass `--bg_excl_index <path>` to skip the large background-subtraction step. Candidate indexes are species/group-specific, so they can't be reused across runs.

## Software versions

| Software | Version | Container | Used by |
|----------|---------|-----------|---------|
| GGCAT | 2.2.0 | `quay.io/biocontainers/ggcat:2.2.0--hf1b6044_0` | index building |
| SBWT (sbwt-rs-cli) | 0.4.2 (patched, `-f93d92c`) | Sanger-internal `.sif` | index building, set-diff |
| Themisto2 | 0.0.1 | `quay.io/sangerpathogens/themisto2:0.0.1` | index building |
| pandas | 2.2.1 | `quay.io/sangerpathogens/pandas:2.2.1` | colour mapping, group-specificity filter |
| Biopython | 1.84 | `quay.io/biocontainers/biopython:1.84` | marker post-processing |
| primer3 | 2.6.1 | `quay.io/biocontainers/primer3:2.6.1--pl5321h503566f_7` | primer design |

All software dependencies are containerised (Docker/Singularity).

## Issues and contributions

- **GitHub/GitLab:** Log an issue or open a merge request on this repository
- **Sanger users:** Raise an issue on the PAM Freshservice portal: https://sanger.freshservice.com/support/catalog/items/426