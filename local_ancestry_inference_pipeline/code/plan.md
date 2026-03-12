# Benchmarking Local Ancestry Inference in 1000 Genomes Trio Data

## Scope

This repository currently supports a two-tool benchmark:

- FLARE
- RFMix

The workflow is built around the high-coverage phased 1000 Genomes panel and a
trio-aware child-only masking design for truth-free consistency checks.

## Dataset

- cohort: 1000 Genomes high-coverage phased panel
- pedigree source: `1kGP.3202_samples.pedigree_info.txt`
- benchmark unit: complete trios with autosomal phased genotypes
- chromosome-first workflow: chr22 is the initial validation target

## Benchmark Design

### Core sample sets

`build_sample_sets.sh` generates:

- `trios.complete.tsv`
- `samples.children.tsv`
- `samples.parents.tsv`
- `samples.trio_members.tsv`
- `query.strict.tsv`
- `reference.strict.tsv`
- `query.child_only.tsv`
- `reference.child_only.tsv`

### Active masking condition

The production workflow uses child-only masking:

- query set: complete-trio children
- reference set: non-child samples
- labels: population-level labels staged separately for FLARE and RFMix

### Reference-panel contracts

`build_population_maps.R` produces:

- FLARE reference-panel files: `flare_ref_panel.*.txt`
- RFMix sample maps: `rfmix_sample_map.*.tsv`
- shared metadata maps: `sample_to_population.tsv`, `sample_to_superpopulation.tsv`, `sample_to_pop_superpop.tsv`

## Active Workflow

### Preparation

1. Download phased panel inputs with `scripts/download.sh` if they are not already present.
2. Build trio/query/reference sets with `scripts/build_sample_sets.sh`.
3. Build population and tool-map files with `scripts/build_population_maps.R`.
4. Stage common manifests with `scripts/stage_tool_maps.sh`.
5. Download/build FLARE and RFMix with `scripts/download_tools.sh`.
6. Prepare FLARE and RFMix genetic maps with `scripts/prepare_genetic_maps.sh`.

### Inference

- FLARE is run with `scripts/run_flare_test_child_only.sh`
- RFMix is run with `scripts/run_rfmix_test_child_only.sh`
- the combined multi-child workflow is orchestrated by `Snakefile`

### Integration

`code/scripts/integration/` standardizes tool outputs into a shared interval
representation:

- `flare_integrate.py`
- `rfmix_integrate.py`
- `merge_lai_inputs.py`

The merge step builds:

- unified segment tables
- task-2 Mendelian interval inputs
- task-3 cross-tool pairwise interval inputs

## Evaluation

The downstream evaluation code under
`calculate_mendelian_violation_rate/` consumes the merged task-2 interval files
to compute:

- per-child Mendelian violation burden
- cross-tool summary violation rates
- optional position-binned violation plots

## Outputs

Key workflow outputs include:

- `out/FLARE/...`
- `out/RFMix/...`
- `out/Combined/chr22/combined_tool_summary.tsv`
- `out/Combined/chr22/combined_tool_paths.txt`
- `out/integration/*.segments.tsv`
- `out/integration/*.task2_intervals.tsv`
- `out/integration/*.task3_pairwise_intervals.tsv`

## Reproducibility Notes

- keep configuration repo-relative where possible
- avoid user-specific absolute paths in committed scripts and docs
- treat `snake_samples.yaml` and `snake_conf.yaml` as the primary runtime knobs
- use the generated metadata in `code/meta/` as the contract between
  preprocessing, inference, and evaluation
