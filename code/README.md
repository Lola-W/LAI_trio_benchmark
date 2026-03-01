# Code Directory

This folder contains benchmark planning and preprocessing assets.

## Files and Folders

- `plan.md`: full benchmark design notes.
- `scripts/build_sample_sets.sh`: generates trio/query/reference sample sets from pedigree data.
- `scripts/build_population_maps.R`: builds population/superpopulation maps and tool-specific panel-map files.
- `scripts/stage_tool_maps.sh`: stages tool-map outputs into `../data/`.
- `scripts/download_tools.sh`: downloads/builds FLARE, RFMix, and Gnomix (no pip).
- `scripts/prepare_genetic_maps.sh`: downloads GRCh38 maps and prepares FLARE/RFMix/Gnomix map inputs.
- `scripts/run_flare_test_child_only.sh`: child-only FLARE test script.
- `scripts/run_rfmix_test_child_only.sh`: child-only RFMix test script.
- `scripts/run_gnomix_test_child_only.sh`: child-only Gnomix test script.
- `Snakefile`: Snakemake workflow for chr22 sample selection/masking -> FLARE/RFMix/Gnomix -> combined summary.
- `snake_conf.yaml`: primary Snakemake config (paths, chromosome, threads).
- `snake_samples.yaml`: sample-selection/masking config (population/superpopulation, n_samples, seed samples).
- `envs/lai-tools.conda.yml`: conda environment for tooling/runtime dependencies.
- `docs/reference_panel_and_masking.md`: strict vs child-only masking definitions and label contracts.
- `meta/`: generated preprocessing outputs.

## Script Index (Numbered)

1. `scripts/build_sample_sets.sh`
2. `scripts/build_population_maps.R`
3. `scripts/stage_tool_maps.sh`
4. `scripts/download_tools.sh`
5. `scripts/prepare_genetic_maps.sh`
6. `scripts/run_flare_test_child_only.sh`
7. `scripts/run_rfmix_test_child_only.sh`
8. `scripts/run_gnomix_test_child_only.sh`
9. `Snakefile` + `snake_conf.yaml` + `snake_samples.yaml`

Notes:
1. `download_tools.sh` auto-detects RFMix build style (CMake or Autotools).
2. `prepare_genetic_maps.sh` writes `maps/plink/plink.chr<CHR>.GRCh38.map`, `maps/plink_allchr.GRCh38.map.tsv`, and `maps/gnomix_map_chr<CHR>.tsv`.
3. Run scripts auto-detect map files and print recovery commands when missing.
4. `run_rfmix_test_child_only.sh` auto-detects binary path at `tools/rfmix/build/rfmix` or `tools/rfmix/rfmix`.
5. `run_rfmix_test_child_only.sh` now defaults to benchmark settings:
all child-masked references (`MAX_REF_PER_POP=0`), no extra AF thinning (`REF_MIN_AF=0`), and RFMix built-in defaults for trees/CRF/window unless explicitly overridden (`RFMIX_TREES`, `RFMIX_CRF_SPACING`, `RFMIX_RF_WINDOW_SIZE`). `RFMIX_THREADS` defaults to `THREADS` (or `SLURM_CPUS_PER_TASK`).
6. `run_gnomix_test_child_only.sh` auto-detects chromosome naming (`22` vs `chr22`), preprocesses to chromosome-restricted biallelic SNP inputs, and passes the matching chromosome label to Gnomix to avoid whole-VCF fallback and OOM.
7. `run_gnomix_test_child_only.sh` sets `model.inference: fast` in a per-run config by default (`GNOMIX_MODEL_INFERENCE` to override), and reports an explicit OOM hint on exit code 137.

## Recommended Order

```bash
conda env create -f code/envs/lai-tools.conda.yml
conda activate lai-benchmark-tools

bash code/scripts/build_sample_sets.sh
Rscript code/scripts/build_population_maps.R
bash code/scripts/stage_tool_maps.sh TRIO0001
bash code/scripts/download_tools.sh
bash code/scripts/prepare_genetic_maps.sh

bash code/scripts/run_flare_test_child_only.sh 22
bash code/scripts/run_rfmix_test_child_only.sh 22
bash code/scripts/run_gnomix_test_child_only.sh 22
```

## Snakemake (chr22 first run)

```bash
cd code
snakemake -s Snakefile --configfile snake_conf.yaml --cores 8
```

Cluster submission (generic template):

```bash
cd code
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
    --output=log/slurm-%j.out" \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 120 \
  --use-singularity \
  --singularity-args "-B /path/to/required/bind_mounts" \
  -n
```

Final combined output:
- `out/Combined/chr22/combined_tool_summary.tsv`
- `out/Combined/chr22/combined_tool_paths.txt`

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
- `gnomix_sample_map.*.tsv`
- `summary.population_maps.tsv`
