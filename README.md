# LAI Trio Benchmark

This repository contains a two-stage benchmark for local ancestry inference in
1000 Genomes trio data:

1. a FLARE/RFMix inference pipeline that builds trio-aware benchmark inputs,
   runs local ancestry inference, and standardizes outputs
2. a Mendelian-violation analysis module that scores consistency of inferred
   child ancestry labels against parental population labels

The current benchmark scope is FLARE + RFMix only.

## Repository Layout

- `local_ancestry_inference_pipeline/`: preprocessing, inference, integration,
  and combined summaries
- `calculate_mendelian_violation_rate/`: downstream violation-rate summaries,
  cross-tool comparison, and plotting scripts

## Workflow Overview

### 1. Build benchmark inputs

From the repository root:

```bash
conda env create -f ./local_ancestry_inference_pipeline/code/envs/lai-tools.conda.yml
conda activate lai-benchmark-tools

bash ./local_ancestry_inference_pipeline/code/scripts/build_sample_sets.sh
Rscript ./local_ancestry_inference_pipeline/code/scripts/build_population_maps.R
bash ./local_ancestry_inference_pipeline/code/scripts/stage_tool_maps.sh TRIO0001
```

### 2. Prepare tools and maps

```bash
bash ./local_ancestry_inference_pipeline/code/scripts/download_tools.sh
bash ./local_ancestry_inference_pipeline/code/scripts/prepare_genetic_maps.sh
```

### 3. Run inference

Direct test runs:

```bash
bash ./local_ancestry_inference_pipeline/code/scripts/run_flare_test_child_only.sh 22
bash ./local_ancestry_inference_pipeline/code/scripts/run_rfmix_test_child_only.sh 22
```

Or run the Snakemake workflow:

```bash
cd ./local_ancestry_inference_pipeline/code
snakemake -s Snakefile --configfile snake_conf.yaml --cores 8
```

### 4. Standardize and merge tool outputs

From the repository root:

```bash
./local_ancestry_inference_pipeline/code/scripts/integration/flare_integrate.py \
  --flare-dir ./local_ancestry_inference_pipeline/out/FLARE/child_only_population/chr22 \
  --threads 12 \
  --trio-tsv ./local_ancestry_inference_pipeline/data/FLARE/test/trio.test.tsv \
  --out-prefix ./local_ancestry_inference_pipeline/out/integration/flare_chr22

./local_ancestry_inference_pipeline/code/scripts/integration/rfmix_integrate.py \
  --rfmix-msp ./local_ancestry_inference_pipeline/out/RFMix/child_only_population/chr22/HG00405/test_chr22.msp.tsv \
  --trio-tsv ./local_ancestry_inference_pipeline/data/RFMix/test/HG00405/trio.test.tsv \
  --out-prefix ./local_ancestry_inference_pipeline/out/integration/rfmix_chr22_HG00405

./local_ancestry_inference_pipeline/code/scripts/integration/merge_lai_inputs.py \
  --rfmix-segments ./local_ancestry_inference_pipeline/out/integration/rfmix_chr22_HG00405.segments.tsv \
  --flare-segments ./local_ancestry_inference_pipeline/out/integration/flare_chr22.segments.tsv \
  --out-prefix ./local_ancestry_inference_pipeline/out/integration/merged_chr22
```

## Mendelian-Violation Analysis

The analysis module expects task-2 interval files produced by the merge step.

Portable entry points:

- `calculate_mendelian_violation_rate/calculate_mvr.py`
- `calculate_mendelian_violation_rate/compare_tools_mvr.py`
- `calculate_mendelian_violation_rate/plot_violations.py`

Shell wrappers:

- `run_calculate_mvr.sh`
- `run_compare_tools_mvr.sh`

Please call the Python scripts directly with explicit input/output paths.

Per-child violation burden:

```bash
python3 ./calculate_mendelian_violation_rate/calculate_mvr.py \
  ./local_ancestry_inference_pipeline/out/integration/merged_chr22.task2_intervals.tsv \
  ./local_ancestry_inference_pipeline/code/meta/sample_to_pop_superpop.tsv \
  ./local_ancestry_inference_pipeline/1kGP.3202_samples.pedigree_info.txt \
  -o ./calculate_mendelian_violation_rate/results/mendelian_violation_rates_output.tsv
```

Cross-tool comparison:

```bash
python3 ./calculate_mendelian_violation_rate/compare_tools_mvr.py \
  ./local_ancestry_inference_pipeline/out/integration/merged_chr22.task2_intervals.tsv \
  ./local_ancestry_inference_pipeline/code/meta/sample_to_pop_superpop.tsv \
  ./local_ancestry_inference_pipeline/1kGP.3202_samples.pedigree_info.txt \
  -o ./calculate_mendelian_violation_rate/results/tool_comparison.csv
```

Optional plotting:

```bash
python3 ./calculate_mendelian_violation_rate/plot_violations.py \
  -i ./calculate_mendelian_violation_rate/results/tool_comparison.csv \
  -o ./calculate_mendelian_violation_rate/results/violation_plot.png
```

## Key Outputs

- `local_ancestry_inference_pipeline/out/Combined/chr22/combined_tool_summary.tsv`
- `local_ancestry_inference_pipeline/out/Combined/chr22/combined_tool_paths.txt`
- `local_ancestry_inference_pipeline/out/integration/*.segments.tsv`
- `local_ancestry_inference_pipeline/out/integration/*.task2_intervals.tsv`
- `calculate_mendelian_violation_rate/results/*.tsv`
- `calculate_mendelian_violation_rate/results/*.csv`

## Design Principles

- clear separation between inference, integration, and evaluation
- child-only masking for trio-based truth-free validation
- reproducible metadata contracts in `local_ancestry_inference_pipeline/code/meta/`
