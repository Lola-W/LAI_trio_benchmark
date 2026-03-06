# 1000 Genomes LAI Benchmark (Preprocessing Stage)

This repository prepares trio-aware sample definitions for benchmarking local ancestry inference (LAI) algorithms on the 1000 Genomes high-coverage phased panel.

Current scope:
- Build reproducible trio/query/reference sample lists from `1kGP.3202_samples.pedigree_info.txt`.
- Define strict vs child-only masking conditions.
- Prepare metadata contracts for FLARE, Gnomix, and RFMix runs.

Not included yet:
- Full HPC execution pipeline for FLARE/Gnomix/RFMix.
- Final benchmarking metrics implementation.
- Raw per-chromosome inputs and generated outputs (kept local, ignored by git).

## Project Layout

- `1kGP.3202_samples.pedigree_info.txt`: downloaded pedigree source.
- `chr*/`: per-chromosome phased VCF data (local only; ignored by git).
- `code/plan.md`: benchmark planning document.
- `code/scripts/build_sample_sets.sh`: reproducible preprocessing script.
- `code/scripts/build_population_maps.R`: population/superpopulation and tool-map generator.
- `code/scripts/stage_tool_maps.sh`: copy generated maps into `data/FLARE`, `data/RFMix`, `data/Gnomix`.
- `code/scripts/download_tools.sh`: download/build FLARE, RFMix, Gnomix (no pip).
- `code/scripts/prepare_genetic_maps.sh`: download GRCh38 maps and prepare FLARE/RFMix/Gnomix map files.
- `code/scripts/run_flare_test_child_only.sh`: child-only FLARE test run.
- `code/scripts/run_rfmix_test_child_only.sh`: child-only RFMix test run.
- `code/scripts/run_gnomix_test_child_only.sh`: child-only Gnomix test run.
- `code/docs/reference_panel_and_masking.md`: benchmark condition definitions.
- `code/meta/`: generated sample lists and QC summaries (created by script).
- `data/`: staged tool-specific maps and test files.

## Setup And Run (Numbered)

```bash
conda env create -f code/envs/lai-tools.conda.yml
conda activate lai-benchmark-tools

bash code/scripts/build_sample_sets.sh
Rscript code/scripts/build_population_maps.R
bash code/scripts/stage_tool_maps.sh TRIO0001

# Download tools (no pip)
bash code/scripts/download_tools.sh
bash code/scripts/prepare_genetic_maps.sh

# Child-only test runs (population labels)
bash code/scripts/run_flare_test_child_only.sh 22
bash code/scripts/run_rfmix_test_child_only.sh 22
bash code/scripts/run_gnomix_test_child_only.sh 22
```

1. Create and activate Conda env from [code/envs/lai-tools.conda.yml](code/envs/lai-tools.conda.yml).
2. Build trio/query/reference sample sets with [code/scripts/build_sample_sets.sh](code/scripts/build_sample_sets.sh).
3. Build population/superpopulation maps with [code/scripts/build_population_maps.R](code/scripts/build_population_maps.R).
4. Stage tool-specific map files into `data/` with [code/scripts/stage_tool_maps.sh](code/scripts/stage_tool_maps.sh).
5. Download/build tools with [code/scripts/download_tools.sh](code/scripts/download_tools.sh).
6. Prepare genetic maps with [code/scripts/prepare_genetic_maps.sh](code/scripts/prepare_genetic_maps.sh).
7. Run FLARE child-only test with [code/scripts/run_flare_test_child_only.sh](code/scripts/run_flare_test_child_only.sh).
8. Run RFMix child-only test with [code/scripts/run_rfmix_test_child_only.sh](code/scripts/run_rfmix_test_child_only.sh).
9. Run Gnomix child-only test with [code/scripts/run_gnomix_test_child_only.sh](code/scripts/run_gnomix_test_child_only.sh).

Outputs are written to `code/meta/`, `data/`, `work/`, and `out/`.
