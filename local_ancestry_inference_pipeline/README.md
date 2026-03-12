# 1000 Genomes LAI Benchmark

This repository benchmarks local ancestry inference on the 1000 Genomes
high-coverage phased panel using trio-aware child-only masking.

Current scope:
- build reproducible trio/query/reference sample lists from `1kGP.3202_samples.pedigree_info.txt`
- generate FLARE and RFMix panel/sample-map inputs
- run chr-wise FLARE and RFMix test workflows directly or through Snakemake
- standardize FLARE and RFMix outputs for downstream Mendelian and cross-tool analysis

Not included in git:
- raw per-chromosome VCF inputs
- generated `tools/`, `work/`, and `out/` outputs

## Project Layout

- `1kGP.3202_samples.pedigree_info.txt`: pedigree source file
- `code/`: workflow, scripts, configs, metadata, and integration utilities
- `data/`: staged FLARE/RFMix sample maps and per-child test manifests
- `maps/`: prepared GRCh38 map files for FLARE and RFMix
- `tools/`: downloaded FLARE jar and RFMix build tree
- `work/`: per-tool working directories
- `out/`: FLARE/RFMix outputs plus combined summaries

## Quick Start

```bash
conda env create -f ./code/envs/lai-tools.conda.yml
conda activate lai-benchmark-tools

bash ./code/scripts/build_sample_sets.sh
Rscript ./code/scripts/build_population_maps.R
bash ./code/scripts/stage_tool_maps.sh TRIO0001

bash ./code/scripts/download_tools.sh
bash ./code/scripts/prepare_genetic_maps.sh

bash ./code/scripts/run_flare_test_child_only.sh 22
bash ./code/scripts/run_rfmix_test_child_only.sh 22
```

`./code/snake_samples.yaml` now defaults to `population: "ALL"` and
`superpopulation: "ALL"`, so the Snakemake workflow runs every child with a
complete trio unless you narrow that config.

## Snakemake

```bash
cd ./code
snakemake -s Snakefile --configfile snake_conf.yaml --cores 8
```

Final combined outputs:
- `out/Combined/chr22/combined_tool_summary.tsv`
- `out/Combined/chr22/combined_tool_paths.txt`

## Integration Utilities

The repo now includes FLARE/RFMix normalization and merge helpers under
`code/scripts/integration/`:
- `rfmix_integrate.py`
- `flare_integrate.py`
- `merge_lai_inputs.py`

See `code/scripts/integration/README.md` for example commands.

## More Detail

`code/README.md` documents the workflow files, helper scripts, and cluster
submission example in more detail.
