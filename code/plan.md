# Benchmarking Local Ancestry Inference in 1000 Genomes Trio Data

## Background and assumptions

Local ancestry inference (LAI) assigns an ancestry label (and sometimes a probability distribution over labels) to each genomic region of a phased individual, using labeled reference haplotypes and a model of how haplotypes are copied along the genome. citeturn49search4turn31search16turn11search11  In real cohort data, “true” local ancestry is rarely observed, so benchmarking often relies on structured consistency checks and cross-method agreement rather than direct accuracy. citeturn30search2

Your proposed benchmark is a truth-free design based on trios: if we infer LAI in parents and child, then the ancestry carried by each child haplotype at a locus should match the ancestry of the corresponding transmitted parental haplotype (modulo phasing and inference uncertainty). The plan below assumes (i) autosomal analysis only (chr1–22), (ii) phased genotypes, (iii) a harmonized marker set across tools, and (iv) a fixed reference panel that excludes the evaluated child (and, for the “strict” design, excludes the evaluated trio entirely). citeturn49search0turn29view0turn25view0

image_group{"layout":"carousel","aspect_ratio":"16:9","query":["family trio pedigree diagram genetics","local ancestry inference chromosome painting plot","RFMix local ancestry output msp plot","phased haplotypes transmitted from parents to child diagram"],"num_per_query":1}

## Data acquisition and cohort construction

### Choose a primary trio-backed genotype source

The most directly “executable” and trio-aligned option on GRCh38 is the 3,202-sample high-coverage 1000 Genomes phased panel distributed via the entity["organization","International Genome Sample Resource","open human variation resource"] FTP. This resource includes 3,202 samples and *602 complete trios* (with the trio children fully represented) and was produced from high-coverage sequencing. citeturn34view0turn22view0turn21view0

Two practical reasons this dataset is a good match for your trio-based benchmark:

- The panel is already phased, and phasing used pedigree-aware correction via SHAPEIT2 duoHMM for SNVs/INDELs, which is specifically designed to infer inheritance patterns and correct switch errors in pedigrees. citeturn22view0turn30search3turn30search11  
- A simple pedigree-info file for all 3,202 samples is distributed alongside the panel, enabling automated derivation of the 602 complete trios without relying on older Phase 3 pedigree files. citeturn39view0

### Use the UCSC “trio” subset only as a pipeline smoke-test

The entity["organization","UCSC Genome Browser","genome browser project"] GRCh38 “trio” download directory contains **7 example trios** (spanning six populations) intended for visualization of phased inheritance patterns, not the full 602-trio cohort. citeturn13view0turn12view0  
This is still valuable to validate parsing, inheritance-state code, and tool wrappers on a small dataset before scaling. citeturn12view0

### Download inputs

Use the IGSR FTP directory listing (one VCF + index per chromosome) for the phased high-coverage panel. citeturn21view0turn22view0  Example (chr1 shown; repeat for chr1–chr22):

```bash
# High-coverage phased panel (GRCh38), per-chromosome VCFs (SNV/INDEL/SV)
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

# The panel README (documents filtering + phasing)
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/README_1kGP_phased_panel_110722.pdf
```

citeturn21view0turn22view0

Download the accompanying pedigree-info file for the 3,202 samples:

```bash
wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt
```

citeturn39view0

Download a GRCh38 genetic map in PLINK format (used by FLARE; also compatible with RFMix map requirements and can be converted for Gnomix training):

```bash
wget -c https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh38.map.zip
unzip -q plink.GRCh38.map.zip -d maps/
```

citeturn31search0turn24view0turn29view0turn25view0

### Build trio, reference, and query sample lists

From `1kGP.3202_samples.pedigree_info.txt`, define:

- **Children**: samples with both `fatherID != 0` and `motherID != 0` (these correspond to the 602 complete trios). citeturn39view0turn34view0  
- **Trio members**: union of children + their parents. (The expanded resource documentation notes that not all trio members are unique across pedigrees; some individuals participate in more than one trio, so “trio-member count” is not exactly 602×3.) citeturn14search13turn34view0  
- **Reference candidates**: samples not in the trio-member set, stratified by ancestry label.

Because the pedigree-info file does not include population labels, use the `kgp` R package metadata (object `kgpe`) to export a clean table containing family ID, paternal/maternal IDs, sex, population code, and superpopulation for all 3,202 samples. citeturn32search8turn45search13

Example export:

```r
# R
install.packages("kgp")
library(kgp)
data(kgpe)

# Save metadata
write.table(kgpe, file="meta/kgpe_3202.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Also export a 2-column population map convenient for RFMix/Gnomix:
# sample_id \t label
write.table(kgpe[,c("id","pop")], file="meta/sample_to_pop.tsv",
            sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
```

citeturn32search8turn45search13

## Reference-panel definition and masking strategy

### Define the ancestry label space you will benchmark

You can run two benchmark “tracks” that share the same trio consistency metric, but differ in what the LAI label means:

**Population-label track (fine labels):** infer among the population codes available in the cohort metadata (up to 26 populations in the expanded 3,202-sample resource). citeturn34view0turn45search13  

**Superpopulation-label track (coarse labels):** infer among AFR/EUR/EAS/SAS/AMR (the cohort is organized into these five superpopulation groups). citeturn37view0turn34view0  

If your scientific intent is to mirror “the family populations most associated with trios” that appear in the UCSC trio visualization (CEU, YRI, CHS, KHV, PUR, MXL), treat that as a *six-population subset* of the population-label track. citeturn12view0turn13view0

### Implement masking levels as separate benchmark conditions

To make the benchmark defensible and interpretable, implement two masking conditions:

**Strict masking (recommended default):** exclude *all trio members* (parents + child) from the reference VCF/panel. This prevents leakage where a child’s haplotype segments are essentially present in the reference via parents.  

**Child-only masking (your original design):** exclude only the child from the reference, allowing parents to remain in reference panels. This is a valid stress-test of “reference leakage sensitivity,” but it is not an independence-preserving evaluation.

Both are implementable in RFMix and FLARE because they allow excluding samples via a reference-panel/sample-map file without having to rewrite the reference VCF. citeturn29view0turn24view0

### Build tool-specific “panel map” files from exported metadata

- **FLARE** needs a `ref-panel` file with two whitespace-delimited columns: `<sample_id> <panel_name>`, and it will ignore samples not listed. It also explicitly recommends that a reference panel “should not normally contain admixed samples.” citeturn24view0  
- **RFMix** needs a tab-delimited sample map `<sample_id>\t<label>`, and it supports excluding reference samples by prefixing with `#` or `^`. citeturn29view0  
- **Gnomix** (training from scratch) takes a reference VCF plus a sample map file mapping reference samples to reference populations. citeturn25view0turn48view0  

Concretely, generate three files per labeling scheme:

- `samples.query.tsv` (all trio members you want LAI for: parents + children)
- `samples.ref.tsv` (all reference candidates not in query)
- `ref.sample_map.tsv` (two columns: sample, label)

## Running FLARE, Gnomix, and RFMix in a unified HPC pipeline

### Harmonize the variant set across tools

To avoid benchmarking “different inputs,” standardize your input marker set. A pragmatic, widely compatible choice is:

- autosomes only (chr1–22)
- `FILTER=PASS`
- biallelic SNPs only
- remove duplicate positions
- enforce **no missing genotypes** in reference+query after filtering (or impute missing genotypes)

This aligns with FLARE’s strict requirement that all genotypes be phased and have *no missing alleles*. citeturn49search0  It also avoids relying on RFMix’s missing-data support, which it explicitly states has not been well-tested. citeturn49search9

The high-coverage panel README documents that variants were filtered before phasing (including `FILTER=PASS`, genotype missingness <5%, HWE p-value threshold, Mendelian error rate threshold, and MAC≥2), but that does **not** guarantee missingness is zero, so you still need an explicit “no-missing” step for FLARE comparability. citeturn22view0turn23view0

Example preprocessing per chromosome:

```bash
# Inputs:
#   raw.vcf.gz = 1kGP_high_coverage_Illumina.chr${CHR}...vcf.gz
#   samples.ref.tsv, samples.query.tsv
# Outputs:
#   ref.bcf, query.bcf (or vcf.gz), plus a common site list if needed

CHR=1
RAW=vcf/raw/chr${CHR}.vcf.gz

# Split into reference and query sets
bcftools view -S meta/samples.ref.tsv   -Oz -o work/chr${CHR}.ref.vcf.gz   ${RAW}
bcftools view -S meta/samples.query.tsv -Oz -o work/chr${CHR}.query.vcf.gz ${RAW}

# Restrict to PASS, biallelic SNPs, and drop sites with any missing GT in either ref or query
bcftools view -f PASS -m2 -M2 -v snps work/chr${CHR}.ref.vcf.gz   -Oz -o work/chr${CHR}.ref.filt.vcf.gz
bcftools view -f PASS -m2 -M2 -v snps work/chr${CHR}.query.vcf.gz -Oz -o work/chr${CHR}.query.filt.vcf.gz

# Optional: intersect sites explicitly (RFMix does it internally; FLARE expects consistent maps/positions)
bcftools isec -n=2 -w1 -Oz -o work/chr${CHR}.ref.common.vcf.gz   work/chr${CHR}.ref.filt.vcf.gz work/chr${CHR}.query.filt.vcf.gz
bcftools isec -n=2 -w2 -Oz -o work/chr${CHR}.query.common.vcf.gz work/chr${CHR}.ref.filt.vcf.gz work/chr${CHR}.query.filt.vcf.gz

# Index
tabix -p vcf work/chr${CHR}.ref.common.vcf.gz
tabix -p vcf work/chr${CHR}.query.common.vcf.gz
```

citeturn29view0turn49search0

If you prefer to impute sporadic missing genotypes rather than drop sites, use Beagle: it is explicitly designed to phase and impute missing genotypes, and FLARE’s documentation points to Beagle as the intended solution when missing/unphased genotypes exist. citeturn49search0turn49search7

### Tool installation and pinned versions

**FLARE**: current version 0.6.0 (updated Nov 3, 2025) and requires Java 11+. citeturn24view0  
**Gnomix**: Python package installed from the repository; training-from-scratch requires a GRCh38-trained model if you use GRCh38 genotypes. The distributed “default” pre-trained models are trained on hg37 references and include multiple biogeographic labels, so they are not a clean drop-in for GRCh38. citeturn25view0turn48view0  
**RFMix v2**: compile from source; requires phased VCF/BCF reference and query, sample map, genetic map, output basename, and chromosome argument. It intersects markers between reference and query internally. citeturn29view0turn26view0

Licensing must be recorded as part of the benchmark metadata:

- **RFMix v2** is free for academic research use only; others must obtain a license via entity["organization","Stanford Office of Technology Licensing","technology licensing office"]. citeturn26view0turn29view0  
- **Gnomix** similarly includes an academic-only notice in its repository documentation (with commercial use routed via entity["company","Illumina","sequencing company"]? no—Gnomix points to Galatea Bio / Stanford licensing in its README; record the notice as written). citeturn25view0  

### Run commands in a chromosome-parallel Slurm layout

A resource-efficient layout is one Slurm job-array per tool, indexed by chromosome, because each tool naturally runs per chromosome and the IGSR distribution is per chromosome. The VCF file sizes range from hundreds of MB to a few GB per chromosome, so per-chromosome parallelization is both I/O- and memory-manageable. citeturn21view0

**FLARE** (per chromosome):
- Inputs: reference VCF, study/query VCF, reference-panel file, PLINK genetic map (`map=`), output prefix. citeturn24view0

```bash
java -Xmx32g -jar flare.jar \
  ref=work/chr${CHR}.ref.common.vcf.gz \
  ref-panel=meta/flare_ref_panel.txt \
  gt=work/chr${CHR}.query.common.vcf.gz \
  map=maps/plink.chr${CHR}.GRCh38.map \
  out=out/flare/chr${CHR}
```

citeturn24view0turn31search0

**RFMix v2** (per chromosome):
- Required options include `-f` (query), `-r` (reference), `-m` (sample map), `-g` (genetic map), `-o` (output), and `--chromosome=...`. citeturn29view0

```bash
rfmix \
  -f work/chr${CHR}.query.common.vcf.gz \
  -r work/chr${CHR}.ref.common.vcf.gz \
  -m meta/rfmix_sample_map.tsv \
  -g maps/plink_allchr.GRCh38.map.tsv \
  -o out/rfmix/chr${CHR} \
  --chromosome=chr${CHR}
```

Notes you must implement explicitly:
- RFMix expects the genetic map as tab-delimited with at least three columns: chromosome, bp position, cM position, and a genome-wide map file is acceptable. citeturn29view0turn31search0  
- It will compute the intersection of SNPs between query and reference automatically, but you should still standardize sites across tools to avoid input-driven differences. citeturn29view0  

**Gnomix** (train + infer, per chromosome):
- If training from scratch: `python3 gnomix.py <query> <out> <chr> <phase> <genetic_map> <reference> <sample_map>`; map is a 3-column TSV without headers (or headers starting with `#`). citeturn48view0turn25view0

```bash
python3 gnomix.py \
  work/chr${CHR}.query.common.vcf.gz out/gnomix/chr${CHR} ${CHR} False \
  meta/gnomix_map_chr${CHR}.tsv \
  work/chr${CHR}.ref.common.vcf.gz \
  meta/gnomix_sample_map.tsv
```

citeturn48view0turn25view0

### Standardize outputs into one canonical representation

To compute trio transmission consistency and cross-tool concordance, convert all outputs into a unified structure:

`(tool, trio_id, sample_id, hap ∈ {0,1}, chrom, start_bp, end_bp, ancestry_label, posterior_prob_optional)`

- **FLARE** outputs an ancestry-per-marker representation in a VCF (`AN1`/`AN2` fields for the most likely ancestry per haplotype; probabilities optional). citeturn24view0  
- **RFMix** outputs `.msp.tsv` (Viterbi maximum-likelihood ancestry per CRF segment) and `.fb.tsv` (forward-backward probabilities). citeturn29view0  
- **Gnomix** outputs `.msp` and `.fb` analogs and can optionally output an SNP-level `.lai` if enabled. citeturn48view1turn25view0  

This standardization step is where you ensure “marker-level” and “segment-level” concordance metrics are computed on the same coordinate system.

## Trio-based evaluation metrics

### Infer parental origin along each child haplotype

You need an inheritance-state track to decide, at each locus, which *parental haplotype* contributed to each child haplotype. A practical, robust approach is to run a simple inheritance HMM per trio per chromosome:

- Hidden state: `(f ∈ {0,1}, m ∈ {0,1})`, indicating which of the two paternal haplotypes and which of the two maternal haplotypes are transmitted at that marker.
- Emission: agreement between child hap alleles and the state-implied parental hap alleles (allow a small “error” rate to accommodate residual genotype/phasing errors).
- Transition: recombination-driven switch probability informed by genetic distance from the GRCh38 map.

This approach matches the conceptual purpose of duoHMM (combine haplotypes with pedigree structure to infer inheritance patterns, detect recombination, and correct errors), and it is well-supported by the literature on pedigree-aware phasing. citeturn30search3turn30search11turn30search7

A critical implementation detail: duoHMM-corrected haplotypes can encode a consistent paternal/maternal ordering for children in some outputs, but you should **not** assume that the VCF haplotype ordering is paternal-first without checking; instead, validate ordering on a small set of informative loci and fall back to the inheritance HMM consistently for all trios. citeturn30search7turn22view0

### Mendelian ancestry transmission violation rate

For each trio, chromosome, and tool:

1. Use the inheritance-state track to label each marker (or segment) as belonging to a transmitted parental haplotype: paternal hap `f` and maternal hap `m`.  
2. Compare the tool’s inferred ancestry for:
   - child haplotype segment inherited from father vs father haplotype `f` at that segment,
   - child haplotype segment inherited from mother vs mother haplotype `m` at that segment.

Define the *ancestry transmission violation rate* as:

- **Marker-weighted**: fraction of markers where child transmitted-hap ancestry ≠ parent transmitted-hap ancestry.
- **Basepair-weighted**: fraction of genomic length (bp) where mismatch occurs (preferable if tools output segments).

Because the high-coverage panel has already been filtered using an explicit Mendelian error-rate threshold at the variant level, gross genotype-level trio inconsistencies should already be reduced; remaining violations are therefore more likely to reflect LAI disagreement, residual phasing quirks, or reference mismatch rather than obvious genotype errors. citeturn22view0turn23view0

Report stratifications that make the benchmark interpretable:

- by trio population / superpopulation label
- by chromosome
- by parental heterozygosity density (proxies for informativeness)
- by inferred ancestry posterior probability (where available)

### Cross-tool concordance

Because the true ancestry is unknown, compute *pairwise* and *three-way* agreement across FLARE, Gnomix, and RFMix on the standardized output representation.

Recommended concordance outputs (all computed per haplotype and then aggregated):

- **Basepair agreement** between tools A and B: proportion of bp where labels match.
- **Cohen’s κ** on ancestry labels computed on a fixed grid (e.g., per-marker or per 10kb bins).
- **Boundary concordance**: Jaccard similarity on sets of ancestry-labeled intervals (especially good for detecting over-fragmentation differences).

Also report per-tool “switch rate” (number of ancestry transitions per Morgan or per Mb), since HMM/CRF smoothing differences can produce similar global ancestry but very different tract fragmentation.

## Reproducibility, compute strategy, and deliverables

### Record data provenance and filtering parameters

The phased panel README documents the panel composition (3,202 samples), the full set of variant-type counts, and the phasing approach (SHAPEIT2 duoHMM for SNVs/INDELs; SHAPEIT4 scaffolding for SVs; chrX handled separately), as well as the filtering criteria applied before phasing. Save this README verbatim in your benchmark repository and cite it in outputs. citeturn22view0turn23view0

The cohort publication also documents that the expanded cohort consists of 2,504 original + 698 added related samples completing 602 trios, and that the cohort spans 26 populations and five superpopulation labels. citeturn34view0

### Pin software and environment

For each tool, persist:

- exact version/hash (e.g., FLARE “version 0.6.0” and update date) citeturn24view0  
- command line used (including seeds/threads where applicable) citeturn24view0turn29view0turn25view0  
- genetic map provenance (Beagle GRCh38 PLINK map zip name and timestamp as available in the public index) citeturn31search0  
- licensing notices for RFMix and Gnomix included in your benchmark documentation citeturn26view0turn25view0  

### Final deliverables

A complete benchmark should emit, at minimum:

- A per-trio summary table: violation rates (paternal, maternal, combined), per tool, per chromosome, plus QC flags.
- Pairwise cross-tool concordance summaries (bp agreement, κ) per tool-pair, stratified by population label.
- A reproducible “run manifest” capturing all file checksums from the IGSR FTP indices (VCFs + `.tbi` + pedigree-info), enabling exact re-runs. citeturn21view0turn39view0