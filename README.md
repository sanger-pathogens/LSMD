# lsmd

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.04.0-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

[[_TOC_]]

## Pipeline overview

lsmd (**l**ineage-**s**pecific **m**arker **d**iscovery) finds k-mer markers that identify a target lineage within a species -- specific enough that they aren't shared with sibling lineages of the same species, *and* aren't shared with a much larger external background collection (e.g. AllTheBacteria/ATB). It's species-agnostic, lineage-agnostic and background-db-agnostic: point it at any species' assemblies plus a lineage/grouping column and a background SBWT index.

The pipeline builds a chain of SBWT/Themisto2 indexes and then subtracts them from one another (`sbwt difference`) until what's left is only the k-mers unique to the target lineage. Each intermediate index has a stage name, used as its `emit:` channel name and in `--outdir` filenames/paths:

| Stage name | Meaning |
| --- | --- |
| `species_index` | Species-wide index -- all genomes of the target species, coloured by `--group_label`. |
| `lineage_index` | Per-lineage index -- the subset of `species_index` for one `--target_groups` lineage. |
| ATB / GTDB | Background index, built independently from a large external genome collection. |
| `bg_excl` | `background − species_index` -- background exclusion set (currently ATB only; GTDB not yet wired in). |
| `xlin_bg` | `species_index − lineage_index` -- k-mers found elsewhere in the species but not in this lineage. |
| `candidate_index` | Threshold-filtered candidates from `lineage_index` (`core`/`relaxed`/`catchall` presence-fraction modes). |
| `lin_cand` | `candidate_index − xlin_bg` -- candidates confirmed specific *within* the species. |
| `markers` | `lin_cand − bg_excl` -- final marker set, also confirmed specific against the outside background. |

Concretely, the pipeline runs as three stages:

1. **`BUILD_COLOR_INDEX`** ([assorted-sub-workflows/themisto2](assorted-sub-workflows/themisto2)) -- maps metadata + assemblies to a Themisto2 colour file, builds `species_index` and (if `--target_groups` is set) `lineage_index`, then filters and rebuilds `lineage_index` down to `candidate_index`.
2. **`SET_DIFF_CALCULATIONS`** ([assorted-sub-workflows/themisto2](assorted-sub-workflows/themisto2)) -- computes `bg_excl`, `xlin_bg`, `lin_cand` and `markers` by chaining `sbwt difference` (each diff is immediately re-verified with `sbwt check`, since a corrupted diff has been observed to exit `0`).
3. **`POST_PROCESS_MARKERS`** ([modules/post_processing_markers.nf](modules/post_processing_markers.nf), opt-in via `--primer_post_processing`) -- dumps `markers`' unitigs to FASTA and filters/masks them for PCR/primer-design suitability (length + global/local GC%).

> This mirrors the numbered `00`-`09` stage documentation kept alongside the pipeline's working data (outside this repo) -- see that doc set for the full biological rationale and worked examples behind each stage.

### Current development status

This pipeline is still early in development (see open TODOs in [main.nf](main.nf) and the sub-workflows):

- **Single run per invocation.** `--metadata`/`--assembly_input` only support one species/lineage-set per run; multi-species support (a samplesheet of `species_id, metadata, assembly, target_groups` rows) is planned but not yet implemented.
- **`--target_groups` is a single global list**, not per-species -- fine for one species per run, but will need to move onto that future samplesheet.
- **Step 09/10 tooling is partial.** [modules/post_processing_markers.nf](modules/post_processing_markers.nf) is wired in behind `--primer_post_processing`; [modules/primer3.nf](modules/primer3.nf) and [modules/bait_capture.nf](modules/bait_capture.nf) are still unfilled module templates for a later primer-design/bait-capture step, not yet called from `main.nf`.
- **GTDB-based background exclusion isn't implemented** -- only the ATB-based `bg_excl`/`markers` path currently runs.
- A couple of diagnostic scripts (`bin/plot_specificity.py`, `bin/validation_classify_gpsc.py`) exist but aren't yet wired into a module/subworkflow.

## Usage

### Quickstart

#### From source code

1. Clone this repository (including submodules):

   ```bash
   git clone --recurse-submodules git@gitlab.internal.sanger.ac.uk:sanger-pathogens/pipelines/lsmd.git
   cd lsmd
   ```

2. Run with `-profile sanger_local` on the Sanger farm (uses Singularity; also raises `max_cpus`/`max_memory` so per-process `cpu_*`/`mem_*` labels -- e.g. GGCAT/Themisto2's `cpu_32`/`mem_64` -- aren't silently clipped under local execution):

   ```bash
   nextflow run main.nf \
       -profile sanger_local \
       --metadata metadata.csv \
       --assembly_input assemblies/ \
       --assembly_suffix .fasta.gz \
       --group_label Lineage \
       --sample_col name \
       --outdir my_output
   ```

   `docker`/`singularity` profiles are also available (inherited from [nextflow-commons](https://github.com/sanger-pathogens/nextflow-commons)).
   :warning: If no profile is specified the pipeline runs with a Sanger HPC-specific configuration, including use of temp storage (`--temp_space`). Configure appropriately for other systems.

3. A bundled `test` profile runs a small (11-genome) *V. tarriae* dataset, exercising the `--target_groups` code path:

   ```bash
   nextflow run main.nf -profile test,sanger_local --outdir test_output
   ```

4. Once the run has finished and you've inspected the output, clean up intermediate files (keep `work/`/`.nextflow.log` until you're satisfied outputs are correct):

   ```bash
   rm -rf work .nextflow*
   ```

#### Using on the Sanger farm

First load the pipeline module:

```bash
module load lsmd
```

Then run on the command line with `lsmd <options>`:

```bash
lsmd --help
```

Submit to LSF using this team's `bsub.py` wrapper (not raw `bsub`) -- sizing here (32 threads, 64GB) matches `-profile sanger_local`'s `max_cpus`/`max_memory`, per the comment in [nextflow.config](nextflow.config#L86):

```bash
module load lsmd bsub.py

bsub.py --threads 32 64 lsmd_run \
    "lsmd -profile sanger_local \
        --metadata metadata.csv \
        --assembly_input assemblies/ \
        --assembly_suffix .fasta.gz \
        --group_label Lineage \
        --sample_col name \
        --outdir my_output"
```

Not using `-profile sanger_local`? Then the pipeline submits one LSF job per process itself (via the inherited `standard` profile), so the outer job just needs to stay alive to do that submitting -- a much smaller request (e.g. `bsub.py 4 lsmd_run "lsmd --metadata ... --outdir my_output"`) is enough.

### Input

#### Metadata (`--metadata`)

A CSV/TSV file with one row per genome assembly, including a sample-identifier column (`--sample_col`, default `Sample_ID`) matched against assembly filenames, and a grouping column (`--group_label`, required) used to colour assemblies into lineages -- e.g. a GPSC column.

#### Assemblies (`--assembly_input`, `--assembly_suffix`)

Either a directory of assembly FASTA files, or a `.txt` file listing one assembly path per line -- detected automatically. When a directory, `--assembly_suffix` (default `.contigs.fasta`) is appended to each `--sample_col` value to form the assembly filename.

#### Target lineages (`--target_groups`, optional)

Comma-separated label(s) from `--group_label`, e.g. `GPSC1,GPSC2`, to build `lineage_index` for -- required for any output past `species_index`, since `SET_DIFF_CALCULATIONS` needs a lineage index to diff against. Leave empty to only build the species-wide index.

#### Background index (`--bg_index`, `--bg_excl_index`)

`--bg_index` (default: a pre-built ATB SBWT index on the Sanger farm) is what `bg_excl`/`markers` are diffed against. If you've already computed `bg_excl` for this species in a previous run, pass it directly via `--bg_excl_index` to skip re-running that hugemem-scale diff. Both must be built at the same `--kmer_size` as the rest of the pipeline's indexes (default `31`).

### Output

Results are written to `--outdir` (default: `./results`). Under `colour_mapping/`, `ggcat/`, `sbwt/` and `themisto2/`, `<ID>` is either the species run's own ID (`species_index`) or a `--target_groups` label like `GPSC1` (`lineage_index` and its `candidate_index` rebuild):

```
results/
├── colour_mapping/<ID>/
│   ├── index_species/
│   └── index_target_group/<group>/
├── ggcat/<ID>/
├── sbwt/
│   ├── <ID>/
│   │   └── candidate/<lineage>/
│   └── <stage>/<ID>/
├── themisto2/<ID>/
│   ├── build/
│   └── export/
├── candidate_filter/<lineage>/
└── candidate_markers/markers_<species>_<lineage>/
```

| Path | Contents |
| --- | --- |
| `colour_mapping/<ID>/index_species/` | `species_file_colors_input.txt`, `species_label_mapping.tsv`, `species_stats.json` |
| `colour_mapping/<ID>/index_target_group/<group>/` | Same three files, scoped to one `--target_groups` label -- only present when `--target_groups` is set |
| `ggcat/<ID>/` | Unitigs FASTA built from the colour file |
| `sbwt/<ID>/` | SBWT index + LCS array for `species_index`/`lineage_index` |
| `sbwt/<ID>/candidate/<lineage>/` | SBWT index + LCS array for the `candidate_index` rebuild |
| `sbwt/<stage>/<ID>/` | Checked set-difference indexes -- `<stage>` is `bg_excl`, `xlin_bg` or `lin_cand` |
| `themisto2/<ID>/build/` | `index.thm2` |
| `themisto2/<ID>/export/` | `export.unitigs.fa`, `export.color_sets.txt` (optionally gzipped), `export.metadata.txt` |
| `candidate_filter/<lineage>/` | `{lineage}_{min_freq_label}_candidate_unitigs.fasta` + `_stats.txt` -- `candidate_index`, pre-rebuild |
| `candidate_markers/markers_<species>_<lineage>/` | Final `markers` outputs -- see below |

`<ID>` is either the species run's own ID (`species_index`) or a `--target_groups` label like `GPSC1` (`lineage_index` and its `candidate_index` rebuild).

**`candidate_markers/markers_<species>_<lineage>/` contents:**

| File | Written when |
| --- | --- |
| `markers_<species>_<lineage>.sbwt` | Always -- checked `markers` SBWT index |
| `markers_<species>_<lineage>_unitigs.fasta` | Always -- `markers` dumped to FASTA |
| `<lineage>_filtered_markers.fasta` | `--primer_post_processing` |
| `<lineage>_rejected_markers.fasta` | `--primer_post_processing` and `--primer_write_rejected` |
| `<lineage>_marker_analysis.png` | `--primer_post_processing` and `--primer_plot` |

See the [themisto2 sub-workflow README](assorted-sub-workflows/themisto2/README.md) for the full `stats.json` field reference and `lineage_index` directory layout.

### Parameters

Run `nextflow run main.nf --help` for the full, always-up-to-date list (rendered from [schema.json](schema.json) plus the [themisto2 sub-workflow's schema](assorted-sub-workflows/themisto2/schema.json)). Summary below:

**General options**

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--outdir` | `path` | `./results` | Directory where results are written. |
| `--monochrome_logs` | `boolean` | `false` | Output logs in plain ASCII. |
| `--temp_space` | `integer` | `10000` | Temp storage (MB) requested for processes that need it (e.g. GGCAT), via the `request_temp` process label. |

---

**Input options** (colour mapping)

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--metadata` | `path` | `null` | CSV/TSV with one row per assembly, incl. a grouping column. |
| `--sample_col` | `string` | `Sample_ID` | Metadata column matched against assembly filenames. |
| `--group_label` | `string` | `null` | Metadata column used to group assemblies into colours (required). |
| `--assembly_input` | `path` | `null` | Directory of assembly FASTAs, or a `.txt` path-list. |
| `--assembly_suffix` | `string` | `.contigs.fasta` | Suffix appended to `--sample_col` to form the assembly filename (directory input only). |
| `--target_groups` | `string` | `""` | Comma-separated lineage label(s) to also build `lineage_index` for. |

---

**Index build options** (GGCAT / SBWT / Themisto2, shared across `species_index`/`lineage_index`/`candidate_index`)

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--kmer_size` | `integer` | `31` | k-mer size used consistently across GGCAT, SBWT and Themisto2. Must match the background index. |
| `--gzip_export` | `boolean` | `false` | Gzip the Themisto2 export's `color_sets.txt`. |
| `--temp_dir` | `path` | `""` | Scratch root for GGCAT/SBWT temp/working dirs. Falls back to a task-local work dir; only set for full background-DB-scale runs. |
| `--candidate_min_freq` | `string` | `core` | Presence-fraction preset (`core` ≥0.95, `relaxed` ≥0.5, `catchall` ≥1 genome) or a literal fraction, for `candidate_index` filtering. |
| `--candidate_min_genome_count` | `integer` | `5` | Absolute genome-count floor for `candidate_index` filtering, alongside `--candidate_min_freq`. |

---

**Step08 -- set-difference options**

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--bg_index` | `path` | Sanger farm ATB index | Background SBWT index that `bg_excl`/`markers` are diffed against, unless `--bg_excl_index` is set. |
| `--bg_excl_index` | `path` | `""` | An already-computed `bg_excl` index to reuse, skipping the ATB-scale diff. |

---

**Step09 -- candidate marker post-processing** (`--primer_post_processing`)

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--primer_post_processing` | `boolean` | `false` | Filter/mask `markers` for PCR/primer-design suitability. Off by default -- the pipeline stops at the `markers` SBWT index/FASTA otherwise. |
| `--primer_min_length` | `integer` | `100` | Minimum candidate marker length in bp. |
| `--primer_gc_min` | `float` | `35.0` | Minimum global GC% a candidate's whole sequence must fall within. |
| `--primer_gc_max` | `float` | `60.0` | Maximum global GC% a candidate's whole sequence must fall within. |
| `--primer_window_size` | `integer` | `--kmer_size` | Local sliding-window size (bp) for the GC check. Out-of-range windows are soft-masked (lowercased), not rejected. |
| `--primer_write_rejected` | `boolean` | `true` | Write rejected (too-short / out-of-range) candidates to their own FASTA. |
| `--primer_plot` | `boolean` | `true` | Generate the length/GC diagnostic plot. |

### Dependencies

- Nextflow ≥ 21.04.0
- All software dependencies are containerised (Docker/Singularity images).
- Step08 requires a pre-built background SBWT index at the same `--kmer_size` (Sanger HPC default: ATB, via `--bg_index`).

## Software versions

| Software | Version | Image |
| --- | --- | --- |
| GGCAT | 2.2.0 | `quay.io/biocontainers/ggcat:2.2.0--hf1b6044_0` |
| SBWT (sbwt-rs-cli) | 0.4.2 | Sanger-internal `.sif` build (see [sbwt.nf](assorted-sub-workflows/themisto2/modules/sbwt.nf)) |
| Themisto2 | 0.0.1 | `quay.io/sangerpathogens/themisto2:0.0.1` |
| pandas / Biopython / matplotlib | 2.2.1 / 1.87 / 3.10.9 | `quay.io/sangerpathogens/pandas:2.2.1` |

See `assorted-sub-workflows/themisto2/modules/` and `modules/` for pinned container versions.

## Issues and Contributions

**GitHub/GitLab users:** if you find an issue with this pipeline, or would like to suggest an improvement, please log an issue or open a merge request on this repository.

**Sanger users:** if you need internal support, you can raise an issue on the PAM Freshservice portal: https://sanger.freshservice.com/support/catalog/items/426
