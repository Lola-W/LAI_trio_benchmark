# Integration Scripts (RFMix + FLARE)

## Scripts

- `rfmix_integrate.py`
- `flare_integrate.py`
- `merge_lai_inputs.py`

Both parsers emit the same standardized segment schema so downstream task-2 and
task-3 code can consume FLARE and RFMix output with identical column handling.

## Standardized Segment Schema

- `Tool`
- `Chromosome`
- `Start`
- `End`
- `Sample_ID`
- `Haplotype_1_Ancestry`
- `Haplotype_2_Ancestry`
- `Haplotype_1_Code`
- `Haplotype_2_Code`
- `Window_SNP_Count`
- `Source_File`
- `Trio_ID`
- `Child_ID`
- `Father_ID`
- `Mother_ID`
- `Sample_Role`

## RFMix Parser

### Input

- RFMix `*.msp.tsv`
- optional trio metadata TSV with `trio_id child father mother`

### Usage

```bash
./code/scripts/integration/rfmix_integrate.py \
  --rfmix-msp ./out/RFMix/child_only_population/chr22/HG00405/test_chr22.msp.tsv \
  --trio-tsv ./data/RFMix/test/HG00405/trio.test.tsv \
  --out-prefix ./out/integration/rfmix_chr22_HG00405
```

## FLARE Parser

### Input

- FLARE `*.anc.vcf.gz`
- matching `*.model`
- optional trio metadata TSV

### Usage: directory scan

```bash
./code/scripts/integration/flare_integrate.py \
  --flare-dir ./out/FLARE/child_only_population/chr22 \
  --threads 12 \
  --trio-tsv ./data/FLARE/test/trio.test.tsv \
  --out-prefix ./out/integration/flare_chr22
```

### Usage: single child

```bash
./code/scripts/integration/flare_integrate.py \
  --anc-vcf ./out/FLARE/child_only_population/chr22/HG00405/test_chr22.anc.vcf.gz \
  --model ./out/FLARE/child_only_population/chr22/HG00405/test_chr22.model \
  --trio-tsv ./data/FLARE/test/HG00405/trio.test.tsv \
  --out-prefix ./out/integration/flare_chr22_HG00405
```

## Outputs

Given `--out-prefix PREFIX`, both parser scripts write:

- `PREFIX.segments.tsv`
- `PREFIX.ancestry_codes.tsv`
- `PREFIX.trio_availability.tsv` when `--trio-tsv` is provided

## Merge + Task Input Builder

`merge_lai_inputs.py` combines standardized FLARE and RFMix segment files and
then prebuilds:

- task-2 Mendelian interval inputs
- task-3 cross-tool pairwise interval inputs

### Usage

```bash
./code/scripts/integration/merge_lai_inputs.py \
  --rfmix-segments ./out/integration/rfmix_chr22_HG00405.segments.tsv \
  --flare-segments ./out/integration/flare_chr22_HG00405.segments.tsv \
  --out-prefix ./out/integration/merged_chr22_HG00405
```

### Usage with pre-standardized RFMix directory

```bash
./code/scripts/integration/merge_lai_inputs.py \
  --rfmix-dir ./out/integration/rfmix_batches \
  --rfmix-trio-tsv ./code/meta/trios.complete.tsv \
  --flare-segments ./out/integration/flare_chr22.segments.tsv \
  --out-prefix ./out/integration/merged_chr22
```

### Merge Outputs

Given `--out-prefix PREFIX`, the script writes:

- `PREFIX.unified_segments.tsv`
- `PREFIX.task2_intervals.tsv`
- `PREFIX.task2_summary.tsv`
- `PREFIX.task3_pairwise_intervals.tsv`
- `PREFIX.task3_pairwise_summary.tsv`
