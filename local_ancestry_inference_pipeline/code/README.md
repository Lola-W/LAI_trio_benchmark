# Code Directory

This folder contains the FLARE/RFMix benchmark workflow, preprocessing scripts,
and integration utilities for the 1000 Genomes trio benchmark.

## Files and Folders

- `plan.md`: benchmark planning notes
- `Snakefile`: Snakemake workflow for chr-wise sample selection, masking, FLARE, RFMix, and combined summaries
- `snake_conf.yaml`: primary Snakemake config (repo paths, chromosome, threads, RFMix guardrails)
- `snake_samples.yaml`: sample-selection config; defaults to `ALL/ALL` children with complete trios
- `scripts/build_sample_sets.sh`: generates trio/query/reference sample sets from pedigree data
- `scripts/build_population_maps.R`: builds population/superpopulation metadata plus FLARE and RFMix panel-map files
- `scripts/stage_tool_maps.sh`: stages metadata and per-child test manifests into `../data/`
- `scripts/stage_rfmix_inputs_from_flare.sh`: builds an RFMix-ready bundle from refreshed FLARE inputs
- `scripts/download.sh`: downloads 1000 Genomes phased panel inputs into the repo root by default
- `scripts/download_tools.sh`: downloads/builds FLARE and RFMix
- `scripts/prepare_genetic_maps.sh`: prepares GRCh38 FLARE and RFMix genetic maps
- `scripts/run_flare_test_child_only.sh`: child-only FLARE test runner
- `scripts/run_rfmix_test_child_only.sh`: child-only RFMix test runner
- `scripts/integration/`: FLARE/RFMix standardization and merge helpers
- `docs/reference_panel_and_masking.md`: strict vs child-only masking definitions and label contracts
- `envs/lai-tools.conda.yml`: conda environment for workflow/runtime dependencies
- `meta/`: generated preprocessing outputs

## Notes

1. `download_tools.sh` auto-detects the RFMix build style and builds either the
   CMake or Autotools layout.
2. `prepare_genetic_maps.sh` writes `maps/plink/plink.chr<CHR>.GRCh38.map` for
   FLARE and `maps/plink_allchr.GRCh38.map.tsv` for RFMix.
3. `run_rfmix_test_child_only.sh` auto-detects the RFMix binary under
   `tools/rfmix/build/rfmix` or `tools/rfmix/rfmix`.
4. `run_rfmix_test_child_only.sh` keeps full-panel defaults unless you override
   `MAX_REF_PER_POP`, `REF_MIN_AF`, `RFMIX_TREES`, `RFMIX_CRF_SPACING`, or
   `RFMIX_RF_WINDOW_SIZE`.
5. `snake_conf.yaml` exposes `rfmix_max_ref_per_pop` and `rfmix_ref_min_af`
   for OOM-safe Snakemake runs.

## Recommended Order

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

## Snakemake

```bash
cd ./code
snakemake -s Snakefile --configfile snake_conf.yaml --cores 8
```

Cluster submission example:

```bash
cd ./code
mkdir -p log

snakemake -s Snakefile \
  --configfile snake_conf.yaml \
  -j 150 \
  --cluster-config cluster_config.yaml \
  --cluster "sbatch \
    -N {cluster.nodes} \
    -n {cluster.cpus_per_task} \
    --mem={cluster.mem} \
    -t {cluster.time} \
    -p {cluster.partition} \
    -q {cluster.qos} \
    -A {cluster.account} \
    --mail-user={cluster.mail_user} \
    --mail-type={cluster.mail_type} \
    --output=./log/slurm-%j.out" \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 120 \
  -n
```

Add any site-specific `--use-singularity` and `--singularity-args` options that
your cluster requires.

Final combined output:
- `out/Combined/chr22/combined_tool_summary.tsv`
- `out/Combined/chr22/combined_tool_paths.txt`

## Integration Utilities

The workflow now includes:
- `scripts/integration/rfmix_integrate.py`
- `scripts/integration/flare_integrate.py`
- `scripts/integration/merge_lai_inputs.py`

These scripts normalize FLARE and RFMix outputs into a shared segment schema and
then build downstream task-2/task-3 interval inputs. See
`scripts/integration/README.md` for command examples.

## Current Generated Outputs (`meta/`)

- `trios.complete.tsv`
- `samples.children.tsv`
- `samples.parents.tsv`
- `samples.trio_members.tsv`
- `samples.all_3202.tsv`
- `query.strict.tsv`
- `reference.strict.tsv`
- `query.child_only.tsv`
- `reference.child_only.tsv`
- `summary.counts.tsv`

Additional outputs from `build_population_maps.R`:
- `sample_to_population.tsv`
- `sample_to_superpopulation.tsv`
- `sample_to_pop_superpop.tsv`
- `kgpe_3202_metadata.tsv`
- `flare_ref_panel.*.txt`
- `rfmix_sample_map.*.tsv`
- `summary.population_maps.tsv`
