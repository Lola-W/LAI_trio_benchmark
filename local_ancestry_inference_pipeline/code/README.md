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
- `envs/gnomix-legacy.conda.yml`: dedicated Gnomix-compatible legacy environment (pins aligned with `tools/gnomix/requirements.txt`).
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
5. `run_rfmix_test_child_only.sh` standalone defaults remain full-panel (`MAX_REF_PER_POP=0`, `REF_MIN_AF=0`) with RFMix built-in tree/CRF/window defaults unless explicitly overridden (`RFMIX_TREES`, `RFMIX_CRF_SPACING`, `RFMIX_RF_WINDOW_SIZE`); `RFMIX_THREADS` defaults to `THREADS` (or `SLURM_CPUS_PER_TASK`).
6. `run_gnomix_test_child_only.sh` auto-detects chromosome naming (`22` vs `chr22`), preprocesses to chromosome-restricted biallelic SNP inputs, and passes the matching chromosome label to Gnomix to avoid whole-VCF fallback and OOM.
7. `run_gnomix_test_child_only.sh` sets `model.inference: fast` in a per-run config by default (`GNOMIX_MODEL_INFERENCE` to override), and now supports `MAX_REF_PER_POP`, `REF_MIN_AF`, and `GNOMIX_R_ADMIXED` for OOM-safe runs.
8. Snakemake wiring in `snake_conf.yaml` now exposes OOM guardrails (`rfmix_max_ref_per_pop`, `rfmix_ref_min_af`, `gnomix_max_ref_per_pop`, `gnomix_ref_min_af`, `gnomix_r_admixed`) while allowing full-panel restoration by setting them to `0` (and `gnomix_r_admixed: 1`).
9. `run_gnomix_test_child_only.sh` now supports a dedicated Python env via `GNOMIX_ENV_PREFIX` (or `GNOMIX_TOOLS_ENV_PREFIX`). If unset, it auto-detects `/tscc/nfs/home/jiweng/ps-gleesonlab5/user/jiweng/conda-envs/gnomix-legacy` when available; otherwise it falls back to `PYTHON_BIN`/`python3`.

## Recommended Order

```bash
conda env create -f code/envs/lai-tools.conda.yml
conda activate lai-benchmark-tools

# Optional but recommended for Gnomix compatibility:
# conda env create -f code/envs/gnomix-legacy.conda.yml \
#   -p /tscc/nfs/home/jiweng/ps-gleesonlab5/user/jiweng/conda-envs/gnomix-legacy
# export GNOMIX_ENV_PREFIX=/tscc/nfs/home/jiweng/ps-gleesonlab5/user/jiweng/conda-envs/gnomix-legacy

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
cd /tscc/nfs/home/jiweng/ps-gleesonlab9/user/jiweng/1000genomes/code
snakemake -s Snakefile --configfile snake_conf.yaml --cores 8
```

Cluster submission (same style as your existing workflows):

```bash
cd /tscc/nfs/home/jiweng/ps-gleesonlab9/user/jiweng/1000genomes/code
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
    --output=/tscc/nfs/home/jiweng/ps-gleesonlab9/user/jiweng/1000genomes/code/log/slurm-%j.out" \
  --rerun-incomplete \
  --keep-going \
  --latency-wait 120 \
  --use-singularity \
  --singularity-args "-B /tscc/projects/ps-gleesonlab7/,/tscc/projects/ps-gleesonlab8/,/tscc/projects/ps-gleesonlab5/,/tscc/nfs/home/xiy010/,/tscc/nfs/home/chchung/,/tscc/lustre/ddn/scratch/jiweng/,/scratch/jiweng/,/tscc/nfs/home/jiweng/" \
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
