# 🧬 GEL CloudOS GWAS and rare-variant analysis with REGENIE

This repository provides worked examples for analysing Genomics England (GEL) AggV3 whole-genome sequencing data on CloudOS. It uses PLINK2 to prepare a common-variant genotype backbone, REGENIE Step 1 to generate leave-one-chromosome-out (LOCO) predictions, and REGENIE Step 2 for single-variant GWAS or rare-variant set tests.

The example phenotype is **cutaneous melanoma (`CM`)**, coded as a binary trait. The scripts contain analysis-specific phenotype names, covariates, paths, and scientific settings. Treat this repository as an example to adapt and validate for your study.

**The repository is a collection of separately submitted workflows and interactive preparation scripts. Running the root `main.nf` runs rare ncRNA/pseudogene exonic association testing only. It does not run QC, Step 1, or GWAS automatically.**

## 📚 Contents

- [🗺️ Repository map and analysis order](#repository-map-and-analysis-order)
- [🧰 Prerequisites and shared inputs](#prerequisites-and-shared-inputs)
- [☁️ Submitting workflows on CloudOS](#submitting-workflows-on-cloudos)
- [🧬 Common-variant preparation and GWAS](#common-variant-preparation-and-gwas)
- [🎭 Rare-variant annotation and mask preparation](#rare-variant-annotation-and-mask-preparation)
- [🔬 Rare-variant association testing](#rare-variant-association-testing)
- [🔎 Validation and troubleshooting](#validation-and-troubleshooting)
- [🌱 Adapting the example and recording provenance](#adapting-the-example-and-recording-provenance)
- [💌 Contact](#contact)

<a id="repository-map-and-analysis-order"></a>

## 🗺️ Repository map and analysis order

| Stage | Workflow entry script | Companion configuration |
| --- | --- | --- |
| Shared genotype QC | [`01_shared_QC/main.nf`](01_shared_QC/main.nf) | `01_shared_QC/nextflow.config` |
| LD pruning | [`02_LD_pruning/main.nf`](02_LD_pruning/main.nf) | `02_LD_pruning/nextflow.config` |
| REGENIE Step 1 | [`Step.1.REGENIE/main.nf`](Step.1.REGENIE/main.nf) | `Step.1.REGENIE/nextflow.config` |
| Single-variant GWAS | [`Step.2.GWAS/main.nf`](Step.2.GWAS/main.nf) | `Step.2.GWAS/nextflow.config` |
| VEP annotation extraction | [`WGS-RV-pipeline/RV-protein-coding/Extract.VEP.Annotations/scripts/main.nf`](WGS-RV-pipeline/RV-protein-coding/Extract.VEP.Annotations/scripts/main.nf) | `nextflow.config` in that same directory |
| Rare protein-coding association | [`WGS-RV-pipeline/RV-protein-coding/RV-protein-coding.main.nf`](WGS-RV-pipeline/RV-protein-coding/RV-protein-coding.main.nf) | `RV-protein-coding-nextflow.config` in that same directory |
| Rare melanocyte cCRE association | [`WGS-RV-pipeline/RV-cCRE/RV-cCRE.main.nf`](WGS-RV-pipeline/RV-cCRE/RV-cCRE.main.nf) | `RV-cCRE-nextflow.config` in that same directory |
| Rare ncRNA/pseudogene exonic association | [`main.nf`](main.nf) at repository root | [`nextflow.config`](nextflow.config) at repository root |

Supporting files:

- `02_LD_pruning/01_extract_LD_pruned_variants.sh` and `02_merge_HQ_pruned_variants.sh`: interactive extraction and merging of the Step 1 genotypes.
- `Step.2.GWAS/run_REGENIE_step2_chr.sh` and `generate_chr_scripts.sh`: an alternative interactive GWAS route.
- `WGS-RV-pipeline/RV-protein-coding/Extract.VEP.Annotations/scripts/`: VEP extraction, chromosome merging, and coding-mask builders.
- `WGS-RV-pipeline/RV-cCRE/masks/`: cCRE score extraction and mask construction.
- `WGS-RV-pipeline/RV-ncRNA/masks/`: ncRNA GERP/CADD extraction and mask construction. The ncRNA JARVIS extractor is one directory above `masks/`.
- `containers/`: Dockerfiles for PLINK2, REGENIE, and the VEP/annotation utilities.
- `WGS-RV-pipeline/Inputs.md` and `Step.2.GWAS/REGENIE.Step.2.interactive.session` are currently empty; they are not input specifications or executable instructions.

```text
AggV3 chromosome PGENs + site-QC lists + analysis sample list
  -> shared QC -> HQ common-SNP PGENs
  -> LD pruning -> chromosome-specific retained marker lists
  -> interactive extraction -> chromosome-specific pruned PGENs
  -> interactive merge -> genome-wide Step 1 PGEN/PVAR/PSAM
  -> REGENIE Step 1 + phenotype/covariates -> LOCO predictions

Original AggV3 PGENs + site-QC lists + phenotype/covariates + LOCO
  -> REGENIE Step 2 GWAS
  OR
  + prepared masks + gnomAD exclusions -> rare-variant Step 2
       |- protein-coding genes
       |- melanocyte cis-regulatory elements (cCREs)
       `- ncRNA/pseudogene exons
```

The common-SNP/MAF filter belongs to the **Step 1 backbone**. Do not substitute the HQ-common or LD-pruned genotypes/lists for the original genotypes and site-QC lists in rare-variant testing.

<a id="prerequisites-and-shared-inputs"></a>

## 🧰 Prerequisites and shared inputs

### Access and software

You need a GEL CloudOS account, access to the relevant AggV3 data and your study inputs, and a workspace capable of running Nextflow jobs. GEL describes access and platform support in its [CloudOS guide](https://re-docs.genomicsengland.co.uk/cloudos/); detailed CloudOS application documentation is available through the support route linked there.

The workflows use Nextflow DSL2. Existing QC/pruning documentation reports runs with Nextflow `24.04.4`; the rare-variant configurations declare `>=22.10.0`. Record the actual version used by your CloudOS submission rather than assuming every release has been tested.

| Purpose | Configured container image | Interactive requirements |
| --- | --- | --- |
| QC, pruning, extraction, merge | `ghcr.io/shiyuzhang0522/gel-gwas-regenie-plink2:alpha7-20260808` | PLINK2; Linux/Bash utilities |
| REGENIE Steps 1 and 2 | `ghcr.io/shiyuzhang0522/regenie:4.1.2` | REGENIE 4.1.2 for the interactive alternative |
| VEP extraction | `ghcr.io/shiyuzhang0522/gel-vep-extractor:v1.1.0` | Python, bcftools; annotation helpers also use bedtools and `bigWigToBedGraph` |
| Coding masks | No dedicated mask-builder container configured | Python >=3.10, numpy, pandas, **pyarrow** |
| cCRE/ncRNA masks | No dedicated mask-builder container configured | R with `data.table` and `optparse` |

The VEP Dockerfile includes several annotation utilities but does **not** install pyarrow or the R mask-builder dependencies. Check the active interactive environment before starting. The VEP chromosome merger requires Bash >=4.3 and Linux/GNU utilities; these helpers are intended for a CloudOS Linux session, not an unmodified macOS shell.

### Inputs prepared outside this repository

Supply the analysis cohort/keep list, phenotype and covariates, per-chromosome PASS-or-LowMLSQ variant lists, and access to the chromosome genotypes. Cohort definition, ancestry/relatedness decisions, phenotype derivation, and creation of those site-QC lists are not automated here.

Rare analyses additionally require gnomAD exclusion lists and reference resources. The repository does not generate the gnomAD lists, Ensembl gene/exon reference tables, melanocyte cCRE BED, or GERP/JARVIS bigWigs. VEP extraction requires a pre-existing functional-annotation shard manifest and VCF/index files.

### Genotypes, lists, phenotype, and covariates

`pgen_root` must point to the directory immediately above `chrom-1` through `chrom-22`:

```text
<pgen_root>/chrom-N/postproc-pgen/dragen.pgen
<pgen_root>/chrom-N/postproc-pgen/dragen.pvar
<pgen_root>/chrom-N/postproc-pgen/dragen.psam
<variant_list_dir>/chrN.PASS_or_LowMLSQ.variant_ids.txt
<gnomad_exclude_dir>/chrN.variant.list               # rare analyses only
```

Replace `N` with the chromosome number. PGEN/PVAR/PSAM must describe the same dataset. The workflows expect `.pvar`, not `.pvar.zst`.

Site-QC and exclusion lists have **one exact PVAR variant ID per line, without a header**. The site-QC list represents `FILTER=PASS` or `LowMLSQ`. For rare analyses, do not pre-filter it to common variants. The supplied rare-analysis design excludes variants with gnomAD joint NFE AF >0.01, using the same exclusion lists at every tested AAF threshold. An exclusion file may be empty, but it must exist for each selected chromosome.

Example IDs below are illustrative; actual IDs must match your PVAR, including chromosome prefix and alleles. Do not add a `DRAGEN:` prefix or rewrite IDs independently in masks and lists.

```text
chr21:10000001:A:G
chr21:10000002:C:T
```

The PLINK2 keep file contains FID/IID pairs, for example:

```text
sample001 sample001
sample002 sample002
```

The phenotype file has exactly these columns for the implemented example:

```text
FID IID CM
sample001 sample001 0
sample002 sample002 1
sample003 sample003 NA
```

Use whitespace-delimited files (TSV is recommended), with `0=control`, `1=case`, and `NA=missing`. Identifiers here are synthetic. Match FID/IID to the genotype sample identifiers and use a consistent cohort across Step 1 and Step 2.

The covariate header contains **25 columns**: `FID`, `IID`, `genetic_sex`, `study_source`, `year_of_birth`, and `PC1` through `PC20`. Supply one corresponding data row per sample. `genetic_sex` and `study_source` are categorical (`--maxCatLevels 30`); year of birth and PCs are quantitative. All covariate columns are used. The rare workflows validate this exact header composition; they do not perform comprehensive covariate/cohort QC for you.

Step 2 requires the actual **`GEL_CM_REGENIE_step1_1.loco`** file from Step 1. The Nextflow Step 2 workflows create this prediction-list entry inside each task:

```text
CM GEL_CM_REGENIE_step1_1.loco
```

Do not supply the old `*_pred.list` in place of the LOCO file: its paths may refer to a previous task's temporary directory.

<a id="submitting-workflows-on-cloudos"></a>

## ☁️ Submitting workflows on CloudOS

1. Choose a stage from the repository map and record the source commit/tag. Make the required input files and directories available through your CloudOS workspace.
2. Import your chosen GitHub workflow revision using your workspace's supported import procedure. Confirm **both** the entry script and its matching configuration; importing the unchanged repository root selects ncRNA association.
3. The existing LD-pruning guide documents an import interface that expects root-level `main.nf` and `nextflow.config`. For that interface, prepare a **separate deployment copy or branch in your own fork**, placing the chosen stage's script/config at those root names and preserving supporting directories. For coding/cCRE workflows, this also means renaming their differently named entry/config files in that deployment copy. Keep the canonical source files in their original locations. If your workspace supports explicit script/config selection, use the pair in the table instead.
4. Enter the stage's parameter names and values below. Choose file/directory inputs accessible to batch tasks; an interactive session's `/home/vscode/session_data/...` path is not automatically accessible in AWS Batch. Use workspace-resolved locations rather than inventing bucket paths.
5. Confirm the execution backend, queue, work directory, container access, and resource allocation provided by CloudOS. The repository does not supply a complete standalone AWS account/queue setup. The VEP config explicitly sets `awsbatch`, `eu-west-2`, and `workDir='work'`; ensure CloudOS applies the correct workspace work directory and settings.
6. Submit, inspect task logs and published outputs, and use a distinct output directory for each stage/run. Transfer interactive preparation outputs into CloudOS-managed storage before using them as batch inputs.

The blocks below show **workflow parameters**, corresponding to Nextflow `--parameter value` arguments. In a form with separate name/value fields, enter the name without `--`. Replace every `/PATH/...` placeholder. They are not standalone shell commands and do not configure AWS execution.

For a separately configured command-line Nextflow environment, the equivalent invocation is `nextflow run <entry-script> -c <matching-config> --parameter value ...`. Avoid accidentally loading this repository's root ncRNA config into another stage: run from an isolated stage deployment directory and inspect the resolved configuration. Do not use `-C` to discard CloudOS-generated configuration. See [Nextflow configuration](https://www.nextflow.io/docs/latest/config.html) and [command-line options](https://www.nextflow.io/docs/latest/cli.html).

Most workflows expose a `docker` profile for a separately provisioned local Docker environment. There is no general `cloudos` profile defined here; CloudOS supplies its execution settings.

<a id="common-variant-preparation-and-gwas"></a>

## 🧬 Common-variant preparation and GWAS

### 1. Shared QC

Submit `01_shared_QC/main.nf` with its companion config:

```text
--pgen_root /PATH/chrom-msvcf
--variant_list_dir /PATH/site_QC_lists
--keep_file /PATH/GEL_CM_REGENIE.keep.txt
--outdir /PATH/run01/shared_QC
```

| Parameter | Required/default | Meaning |
| --- | --- | --- |
| `pgen_root`, `variant_list_dir`, `keep_file` | Required | Inputs described above |
| `maf` | `0.01` | Minimum minor allele frequency |
| `geno` | `0.01` | Maximum variant missingness |
| `hwe` | `1e-15` | HWE threshold; command uses `--hwe <value> 0` |
| `outdir` | `results` | Published output root |

The workflow always processes autosomes 1–22. It applies the sample keep list and site-QC list, retains A/C/G/T SNPs (`--snps-only just-acgt`), applies the thresholds above, and writes PGEN files plus a SNP list. Configured resources are 4 CPUs and 16 GB per chromosome.

Outputs under `<outdir>/chrN/` are `chrN.HQ_common.{pgen,pvar,psam,snplist,summary.tsv,log}`. Check sample counts, retained variant counts, and logs before pruning. See the [shared-QC guide](01_shared_QC/README_01_shared_QC.md).

### 2. LD pruning

Submit `02_LD_pruning/main.nf` with its companion config:

```text
--hq_pgen_root /PATH/run01/shared_QC
--window 500kb
--r2 0.2
--outdir /PATH/run01/LD_pruning
```

`hq_pgen_root` is required and must directly contain `chr1/` through `chr22/`. `window=500kb`, `r2=0.2`, and `outdir=results` are defaults. PLINK2 runs `--indep-pairwise 500kb 0.2`. Each of the 22 chromosome tasks requests 4 CPUs and 16 GB; the combine task requests 1 CPU and 2 GB.

Outputs under `<outdir>/chrN/` are `chrN.LD_pruned.prune.in`, `.prune.out`, `.summary.tsv`, and `.log`. The intended combined outputs are `REGENIE_step1.LD_pruned.variant_ids.txt`, `LD_pruning.all_chr.summary.tsv`, and `LD_pruning.genomewide.summary.tsv`.

**Known combine-step issue:** `COMBINE_PRUNED_LISTS` is invoked unconditionally, and its `find ... -type f` check can miss Nextflow-staged symlinks. The existing [LD-pruning guide](02_LD_pruning/README_02_LD_pruning.md) reports this failure after all 22 chromosome jobs succeeded. There is no runtime switch to disable the combine process. Verify all chromosome tasks and their outputs before recovering manually.

For example, in Bash in an interactive session with the published results mounted:

```bash
set -euo pipefail
prune_root=/PATH/mounted/LD_pruning
combined=/PATH/writable/REGENIE_step1.LD_pruned.variant_ids.txt
for chr in {1..22}; do
    test -s "$prune_root/chr${chr}/chr${chr}.LD_pruned.prune.in"
    test -s "$prune_root/chr${chr}/chr${chr}.LD_pruned.summary.tsv"
done
for chr in {1..22}; do
    cat "$prune_root/chr${chr}/chr${chr}.LD_pruned.prune.in"
done > "$combined"
wc -l "$combined"
sort "$combined" | uniq -d | wc -l
```

The duplicate count should be zero. Check each summary's retained + removed count against its input count. This recovery produces the combined marker list, not replacement combined summary tables. The extraction script below uses the individual chromosome lists.

### 3. Extract and merge the Step 1 genotypes

Pruning produces **lists**, not the merged genotype dataset required by Step 1. Run the two PLINK2 helper scripts in order inside a Linux interactive session, after arranging their inputs:

```bash
bash 02_LD_pruning/01_extract_LD_pruned_variants.sh
bash 02_LD_pruning/02_merge_HQ_pruned_variants.sh
```

These commands assume the repository is your current directory. Both scripts contain fixed paths; they do not accept path arguments or honour arbitrary environment overrides. Either mount inputs at the expected locations or adapt paths in your own working copies before running them.

| Helper | Existing input location | Existing output location |
| --- | --- | --- |
| Extract | `/home/vscode/session_data/filesystems/results/chrN/chrN.HQ_common.*` and `/home/vscode/session_data/filesystems/chrN/chrN.LD_pruned.prune.in` | `/home/vscode/session_data/REGENIE.GWAS/GWAS.pipeline/Extracted.HQ.Pruned.Var.GT.Step.1/chrN/chrN.HQ_pruned.*` |
| Merge | The extraction output above | `/home/vscode/session_data/REGENIE.GWAS/GWAS.pipeline/REGENIE.Step1.Genotype/REGENIE.Step1.HQ_pruned.*` |

Extraction uses `--extract`, `--make-pgen`, and 4 threads. Merging uses chromosome 1 as the base plus a `--pmerge-list` for chromosomes 2–22, with 8 threads. Verify variant/sample counts and matching genotype file prefixes. Publish the merged `.pgen`, `.pvar`, and `.psam` to CloudOS storage for Step 1. A merged list of IDs alone is insufficient.

### 4. REGENIE Step 1

Submit `Step.1.REGENIE/main.nf` with its companion config:

```text
--pgen /PATH/REGENIE.Step1.HQ_pruned.pgen
--pvar /PATH/REGENIE.Step1.HQ_pruned.pvar
--psam /PATH/REGENIE.Step1.HQ_pruned.psam
--phenotype /PATH/GEL_CM_REGENIE.phenotype.tsv
--covariates /PATH/GEL_CM_REGENIE.covariates.tsv
--outdir /PATH/run01/REGENIE.Step1.results
--bsize 1000
```

All five input parameters are required; the three genotype basenames must match. Note the names **`phenotype`/`covariates` here**, versus **`pheno_file`/`covar_file` in Step 2**. Defaults are `bsize=1000` and `outdir=REGENIE.Step1.results`.

The command uses `--step 1 --bt --lowmem`, phenotype `CM`, and the covariates described above. It publishes files beginning `GEL_CM_REGENIE_step1`, including the LOCO prediction file and prediction list. Carry `GEL_CM_REGENIE_step1_1.loco` forward.

The process body declares 100 GB, but the companion config's `withName: REGENIE_STEP1` selector specifies **30 GB**, with 16 CPUs and 7 days. With that config applied, the selector overrides the process-body resource declaration; check the final CloudOS allocation. These settings are starting allocations, not a guarantee that every dataset fits.

### 5. Single-variant GWAS

Submit `Step.2.GWAS/main.nf` with its companion config:

```text
--pgen_root /PATH/chrom-msvcf
--variant_list_dir /PATH/site_QC_lists
--pheno_file /PATH/GEL_CM_REGENIE.phenotype.tsv
--covar_file /PATH/GEL_CM_REGENIE.covariates.tsv
--step1_loco /PATH/run01/REGENIE.Step1.results/GEL_CM_REGENIE_step1_1.loco
--outdir /PATH/run01/regenie.step2.GWAS
--bsize 400
```

The five inputs are required. Defaults are `bsize=400` and `outdir=regenie.step2.GWAS`. This workflow always runs all 22 autosomes; it has no `chromosomes` parameter. Each task requests 16 CPUs, 40 GB, and 7 days.

Scientific settings are `--bt --minMAC 20 --firth --approx --firth-se --pThresh 0.01`, with LOCO adjustment and site-QC extraction. Results are compressed:

```text
<outdir>/chrN/chrN.GEL_CM_REGENIE_step2_CM.regenie.gz
<outdir>/chrN/chrN.GEL_CM_REGENIE_step2.log
```

The config retries failures up to three times; this overrides the narrower retry rule in the process body. Its `executor.exitReadTimeout` is 10 minutes. Repeated failures still require investigation rather than additional blind retries.

**Interactive alternative:** after adapting its hard-coded paths and providing a prediction list with a LOCO path accessible from the result directory, run `bash Step.2.GWAS/run_REGENIE_step2_chr.sh 21` for chromosome 21. It uses 16 threads. `generate_chr_scripts.sh` writes launchers for all 22 chromosomes under `commands/`; run it from a working directory containing `run_REGENIE_step2_chr.sh`. Generating launchers does not submit or execute them. Avoid launching all chromosomes concurrently unless your session has sufficient resources.

<a id="rare-variant-annotation-and-mask-preparation"></a>

## 🎭 Rare-variant annotation and mask preparation

These preparation steps run separately from association. You may start with already validated masks in the required layouts. CADD-based cCRE/ncRNA preparation uses the chromosome-level VEP tables, so VEP extraction is shared preparation when those tables are not already available.

### 6. Extract VEP annotations and merge by chromosome

Use the VEP workflow/config pair in the repository map. Supply a CSV `manifest` with these exact column names:

```csv
chr,start,end,region,shard,subshard,func_anno_vcf,func_anno_vcf_index
21,10000001,11000000,chr21_10000001_11000000,1,1,/PATH/subshard.vcf.gz,/PATH/subshard.vcf.gz.tbi
```

The row is a schema illustration, not a real GEL shard. Use the actual GEL manifest, including its region naming and complete shard coverage; downstream merging relies on the extracted filenames. VCFs must contain the expected VEP CSQ annotations.

```text
--manifest /PATH/functional_annotation_shards.csv
--outdir /PATH/run01/VEP_extraction
```

`manifest` is required; `outdir` defaults to `results`. One task per subshard requests 4 CPUs, 8 GB, and 48 hours. The summary task requests 1 CPU, 2 GB, and 1 hour. Outputs are:

```text
<outdir>/VEP_annotation/GEL.VEP.shard<shard>.subshard<subshard>.<region>.tsv
<outdir>/VEP_annotation/VEP_annotation_summary.tsv
```

The workflow stages `extract_v8.py` using a default path relative to `projectDir` that includes the full repository subdirectory. Preserve that tree when deploying the entrypoint at root. If running the entrypoint directly from its nested directory, explicitly set `--extract_script` to the accessible `extract_v8.py` file; otherwise the default may resolve to a duplicated, nonexistent subdirectory.

Mount the extraction output and run the chromosome merger in an interactive session:

```bash
bash WGS-RV-pipeline/RV-protein-coding/Extract.VEP.Annotations/scripts/merge_VEP_subshards_by_chromosome.sh
```

Its fixed input is `/home/vscode/session_data/filesystems/VEP_annotation`, and its output is `/home/vscode/session_data/chromosome_level_VEP`. It creates `GEL.VEP.chrN.tsv`, `chromosome_level_VEP.merge_summary.tsv`, and logs. It expects all 22 chromosomes, checks shard/header consistency, and refuses to overwrite existing chromosome outputs. `MAX_PARALLEL` is an environment override (default 4); input/output directory variables are fixed in the script. Inspect logs and summary before using the merged tables.

### 7. Protein-coding masks

Supply chromosome-level VEP TSVs and the Ensembl 114 GRCh38 gene-coordinate TSV used by this example. The coordinate file requires `ensembl_gene_id`, `hgnc_symbol`, `chromosome_name`, `start_position`, and `end_position`.

Use the `chr16fix` builder for this repository's documented analysis:

```bash
python WGS-RV-pipeline/RV-protein-coding/Extract.VEP.Annotations/scripts/02_build_GEL_regenie_rare_coding_masks_chr16fix.py \
  --chrom 21 \
  --vep /PATH/chromosome_level_VEP/GEL.VEP.chr21.tsv \
  --gene-coordinates /PATH/Ensembl114_GRCh38_gene_coordinates.tsv \
  --outdir /PATH/Rare-Coding-Masks
```

Repeat for chromosomes 1–22, updating `--chrom` and `--vep`. Optional `--chunksize` defaults to 1,000,000; retained rows are still accumulated, so chunking does not cap total memory. The original builder remains available, but the `chr16fix` variant additionally excludes **ENSG00000310590 on chromosome 16**, because that gene is absent from the shared coordinate reference. This is an analysis-specific exclusion, not a general requirement for all studies.

The builder selects protein-coding transcript rows using MANE_SELECT or the documented row-specific CANONICAL fallback and creates:

| Mask | Definition in this example |
| --- | --- |
| `pLoF_only` | High-confidence LOFTEE (`LoF=HC`) class |
| `pLoF_Dmis` | pLoF plus damaging class: eligible missense consequences with REVEL >=0.773 or CADD_PHRED >=28.1, or SpliceAI DSmax >=0.20, or `LoF=LC` |
| `Missense` | Broad retained missense/protein-altering consequences |
| `Synonymous` | Synonymous final hierarchical annotation class |

Masks can overlap. The builder preserves original variant IDs and applies no cohort AF, gnomAD AF, or VCF FILTER restriction. Consult its source for transcript and annotation-priority rules. Outputs include an audit Parquet, summary TSV, and these three files per mask, directly under `--outdir`:

```text
chr21.pLoF_only.annotation.txt
chr21.pLoF_only.setlist.txt
chr21.pLoF_only.maskdef.txt
```

### 8. Melanocyte cCRE scores and masks

Required references are `ENCODE4.Epidermal.Melanocyte.cCREs.chr1-22.bed`, the GRCh38 GERP bigWig, `jarvis.bw`, chromosome PVARs, and `GEL.VEP.chrN.tsv`. Arrange the mounted inputs according to each script's path block before running. Commands from repository root for chromosome 21 are:

```bash
bash WGS-RV-pipeline/RV-cCRE/masks/01_retrieve_melanocyte_cCRE_variant_GERP.GEL.sh 21
bash WGS-RV-pipeline/RV-cCRE/masks/02_retrieve_melanocyte_cCRE_variant_CADD.GEL.sh 21
bash WGS-RV-pipeline/RV-cCRE/masks/03_retrieve_melanocyte_cCRE_variant_JARVIS.GEL.sh 21
```

The extractors write under `/home/vscode/session_data/{GERP,CADD,JARVIS}_raw/Melanocyte.cCRE.{GERP,CADD,JARVIS}.scores`, respectively. The R builder expects those directory trees **under `/home/vscode/session_data/filesystems/`**. Persist/remount the score outputs at those expected input paths, or adapt the fixed paths in a separate working copy. Then run:

```bash
Rscript WGS-RV-pipeline/RV-cCRE/masks/04_build_REGENIE_cCRE_masks.GEL.R --chr 21
```

This builder exposes `--chr` only and writes `/home/vscode/session_data/REGENIE_cCRE_masks/chr21/<MASK>/`. Repeat extraction and construction for each chromosome needed by association. Check that all three score tables represent the same variant–cCRE pairs, including missing-score rows.

### 9. ncRNA/pseudogene exonic scores and masks

Use `Ensembl116.GRCh38.autosomal.ncRNA_pseudogene.gene_union_exons.merged.tsv`, chromosome PVARs, chromosome VEP tables, and the GERP/JARVIS bigWigs. Preserve the exon-table schema and coordinate conventions expected by the extractors; this reference is distinct from the Ensembl 114 coding-gene table.

```bash
bash WGS-RV-pipeline/RV-ncRNA/masks/01_retrieve_ncRNA_pseudogene_exonic_variant_GERP.GEL.sh 21
bash WGS-RV-pipeline/RV-ncRNA/masks/02_retrieve_ncRNA_pseudogene_exonic_variant_CADD.GEL.sh 21
bash WGS-RV-pipeline/RV-ncRNA/03_retrieve_ncRNA_pseudogene_exonic_variant_JARVIS.GEL.sh 21
```

Default outputs are `/home/vscode/session_data/GERP_ncRNA`, `CADD_ncRNA`, and `JARVIS_ncRNA`, with filenames `chrN.ncRNA_pseudogene.exonic_variants.<SCORE>.tsv`. The JARVIS script supports environment overrides `PVAR`, `NCRNA_EXONS`, `JARVIS_BW`, and `OUT_DIR`; do not assume the other extractors expose those overrides.

Unlike the cCRE builder, the ncRNA builder accepts explicit input/output roots. For freshly produced scores in the current session:

```bash
Rscript WGS-RV-pipeline/RV-ncRNA/masks/05_build_REGENIE_ncRNA_pseudogene_exonic_masks.GEL.R \
  --chr 21 \
  --input-root /home/vscode/session_data \
  --exon-file /PATH/Ensembl116.GRCh38.autosomal.ncRNA_pseudogene.gene_union_exons.merged.tsv \
  --out-root /PATH/REGENIE_ncRNA_masks \
  --threads 2
```

`--input-root` must directly contain `CADD_ncRNA/`, `GERP_ncRNA/`, and `JARVIS_ncRNA/`. Its default is `/home/vscode/session_data/filesystems`; `--out-root` defaults to `/home/vscode/session_data/REGENIE_ncRNA_masks`. The builder requires concordant variant–gene pairs and consistent metadata across score tables, rejects empty masks, and preserves existing chromosome outputs by failing if they already exist. Use a new output root for reruns.

Both cCRE and ncRNA builders create `CADD` (PHRED >=20), `GERP` (>=2), `JARVIS` (>=0.99), `FUNC_ALL` (union of those score masks), and `ALL` (all eligible pairs irrespective of scores). Missing scores do not qualify for their score-specific masks. Mask construction does not apply the association AAF thresholds.

### Mask formats and directory handoff

Files are headerless and whitespace-delimited. Synthetic ncRNA example for mask `CADD`:

```text
# annotation.txt: variant_ID set_ID annotation_label
chr21:10000001:A:G ENSG00000999999 CADD

# setlist.txt: set_ID chromosome position comma_separated_variant_IDs
ENSG00000999999 21 10000001 chr21:10000001:A:G

# mask definition: mask_name annotation_label
CADD CADD
```

The explanatory `#` lines above are **not** part of the files. cCRE masks use cCRE identifiers instead of gene IDs; coding files use their coding-mask label. These runners expect each selected mask definition to contain one row with matching mask/annotation labels. Retain the builders' set coordinates and variant-to-set mappings.

| Association branch | Set `mask_dir` to | Expected paths below it |
| --- | --- | --- |
| Coding | The coding builder's output directory | `chrN.MASK.annotation.txt`, `chrN.MASK.setlist.txt`, `chrN.MASK.maskdef.txt` |
| cCRE | `REGENIE_cCRE_masks` | `chrN/MASK/chrN.MASK.annotation.txt`, `chrN/MASK/chrN.MASK.setlist.txt`, `chrN/MASK/MASK.mask.def` |
| ncRNA | **`REGENIE_ncRNA_masks/REGENIE_inputs`** | Same nested filenames as cCRE |

Persist these outputs to a CloudOS-accessible location before submitting association. In particular, pointing ncRNA `mask_dir` one level above `REGENIE_inputs` will fail.

<a id="rare-variant-association-testing"></a>

## 🔬 Rare-variant association testing

### Required parameters and scientific settings

All three rare workflows require seven input locations:

| Parameter | Meaning |
| --- | --- |
| `pgen_root` | Original AggV3 chromosome PGEN root |
| `variant_list_dir` | PASS-or-LowMLSQ lists; no common-variant filter |
| `gnomad_exclude_dir` | `chrN.variant.list` files, including empty files when appropriate |
| `pheno_file` | `FID IID CM` phenotype |
| `covar_file` | FID/IID plus the 23 agreed covariates |
| `step1_loco` | Actual Step 1 LOCO file |
| `mask_dir` | Branch-specific root from the table above |

| Optional parameter | Default/allowed values |
| --- | --- |
| `chromosomes` | `1-22`; accepts ranges or comma-separated autosomes, without duplicates/overlaps |
| `masks` | Coding: `pLoF_only,pLoF_Dmis,Missense,Synonymous`; cCRE/ncRNA: `CADD,GERP,JARVIS,FUNC_ALL,ALL` |
| `aaf_thresholds` | `0.01,0.001,0.0001`; select any unique subset |
| `task_cpus`, `task_memory`, `task_time` | `16`, `40 GB`, `7d` per task |
| `regenie_image` | `ghcr.io/shiyuzhang0522/regenie:4.1.2` |
| `outdir` | Coding: `Rare-Coding/Results`; cCRE: `Rare-cCRE/Results`; ncRNA: `Rare-ncRNA/Results` |
| `help` | `false`; `--help` prints a usage summary |

Each task analyses one whole chromosome, one mask, and one AAF threshold. AAF values are proportions: 0.01=1%, 0.001=0.1%, and 0.0001=0.01%. The scripts set both `--aaf-bins` and `--vc-maxAAF` to that task's threshold and use `--vc-tests skato,acato-full --vc-MACthr 10`. Other fixed settings are `--bt --firth --firth-se --approx --pThresh 0.05 --bsize 400 --write-mask-snplist --check-burden-files`.

`--minMAC` and `--build-mask` are left at REGENIE defaults, and standard singleton-mask behaviour is retained. GWAS's explicit `--minMAC 20` and `--pThresh 0.01` are not the rare-analysis settings. Refer to the [REGENIE options documentation](https://rgcgithub.github.io/regenie/options/) when interpreting these tests.

### Pilot, inspect, then expand

For a **coding pilot**, select the coding workflow/config and supply:

```text
--pgen_root /PATH/chrom-msvcf
--variant_list_dir /PATH/site_QC_lists
--gnomad_exclude_dir /PATH/gnomAD_exclusions
--pheno_file /PATH/GEL_CM_REGENIE.phenotype.tsv
--covar_file /PATH/GEL_CM_REGENIE.covariates.tsv
--step1_loco /PATH/GEL_CM_REGENIE_step1_1.loco
--mask_dir /PATH/Rare-Coding-Masks
--chromosomes 21
--masks pLoF_only
--aaf_thresholds 0.01
--outdir /PATH/run01/pilot_coding
```

For a **cCRE pilot**, select the cCRE workflow/config, retain the same six non-mask input parameters, and use:

```text
--mask_dir /PATH/REGENIE_cCRE_masks
--chromosomes 21
--masks CADD
--aaf_thresholds 0.01
--outdir /PATH/run01/pilot_cCRE
```

For an **ncRNA pilot**, select the repository-root workflow/config, retain those six inputs, and use:

```text
--mask_dir /PATH/REGENIE_ncRNA_masks/REGENIE_inputs
--chromosomes 21
--masks CADD
--aaf_thresholds 0.01
--outdir /PATH/run01/pilot_ncRNA
```

Each pilot creates one association task and requires nonempty mask inputs for that selection. After reviewing outputs, start a separate full run with `chromosomes=1-22`, all branch-specific masks, `aaf_thresholds=0.01,0.001,0.0001`, and a new output directory. This gives **264 coding tasks** or **330 tasks for each of cCRE and ncRNA**. There is no workflow-level concurrency cap; CloudOS/AWS Batch scheduling controls how many run together. Check capacity before expanding.

### Published rare-analysis outputs

All branches publish to `<outdir>/<MASK>/AAF_<AAF>/`. Example prefixes are:

```text
Coding: chr21.pLoF_only.AAF_0.01
cCRE:   chr21.Melanocyte_cCRE.CADD.AAF_0.01
ncRNA:  chr21.ncRNA_pseudogene_exonic.CADD.AAF_0.01
```

For each prefix, that directory contains `<prefix>_CM.regenie` and `<prefix>_masks.snplist`, both **uncompressed**. Its `logs/` subdirectory holds the REGENIE `.log`, `.console.log`, resolved `.command.sh`, local `.pred.list`, `.qc.tsv`, and optional `_masks_report.txt`. A burden report may be absent when REGENIE has nothing to report.

`<outdir>/pipeline_info/` contains timestamped Nextflow trace TSV, report HTML, and timeline HTML. Rare-task publishing uses overwrite mode: choose separate output roots to preserve earlier runs even though workflow report filenames are timestamped.

<a id="validation-and-troubleshooting"></a>

## 🔎 Validation and troubleshooting

Validate every handoff before starting the next stage. A file-presence check or completed Nextflow task does not certify that every association converged or produced a usable P value.

| Check or symptom | What to inspect/do |
| --- | --- |
| Wrong analysis starts | Confirm the chosen entrypoint/config pair. Root is ncRNA, not a GWAS driver. |
| Missing input files | Verify directory level, chromosome naming, exact suffixes, and batch-accessible locations. Interactive paths do not automatically carry into batch tasks. |
| QC/pruning completion | Check all 22 logs and summaries, sample counts, retained markers, and retained + removed = input variants. Do not use historical sample counts as expectations for a new cohort. |
| LD combine fails but chromosome jobs succeeded | Inspect the known `find -type f`/symlink issue, validate all chromosome outputs, then combine lists manually as above. |
| Step 1 genotype-prefix error | Supply `.pgen`, `.pvar`, and `.psam` with the same basename from the completed merge. |
| LOCO/prediction-file error | Pass the LOCO file to Nextflow Step 2. For interactive GWAS, ensure the manually prepared prediction-list path is valid from the task's result directory. |
| VEP extractor not found | Check `projectDir` and `extract_script`; preserve the supporting repository tree or explicitly supply the script path. |
| Mask-builder missing input | Check score-table output-to-mount handoff, reference version/schema, and exact directory names. Fixed script variables are not CLI parameters. |
| Coding chr16 coordinate error | Check whether the documented `chr16fix` builder and matching reference were used; do not silently discard additional genes. |
| Rare input-validation failure | Check phenotype coding, 25-column covariate header, mask label, four-column set list, chromosome, and nonempty inputs. Ensure IDs match PVAR exactly. |
| Exit 137/143 in rare workflows | These statuses trigger up to two retries. Investigate memory limits or interruption; retries do not automatically increase memory. Other errors terminate immediately. |
| Missing/empty association output | Read REGENIE and task logs for input, sample-matching, convergence, or filtering problems. GWAS results are `.regenie.gz`; rare results are `.regenie`. |
| Rare QC TSV shows missing P values | Inspect `missing_log10p` by test together with burden reports and REGENIE logs. Row/set counts alone are insufficient. |

Rare QC summaries count result rows, represented genes/cCREs, and missing `LOG10P` values by test. Review masks actually used, exclusions, sample counts, convergence messages, and unexpected empty/small sets. Downstream multiple-testing decisions, plots, and interpretation are separate from these workflows.

Use Nextflow `-resume` only when the relevant cache and work directory remain available, or use the corresponding supported CloudOS resume mechanism. Resubmitting with an empty/new work directory does not recover cached tasks. Preserve failed-task diagnostics before cleanup. See [Nextflow caching and resuming](https://www.nextflow.io/docs/latest/cache-and-resume.html).

<a id="adapting-the-example-and-recording-provenance"></a>

## 🌱 Adapting the example and recording provenance

To analyse a different phenotype/cohort, review cohort construction, binary versus quantitative trait handling, phenotype column names, categorical covariates, and sample matching. `CM`, output names, validation checks, and several scientific flags are hard-coded: changing a filename or inventing a `--phenoCol` workflow parameter will not update them. Adapt the relevant scripts in your own fork, then validate the resulting workflow. The current workflows cover autosomes only.

For different regulatory regions or gene annotations, reconsider the cCRE/exon reference, coordinate conventions, transcript-selection rules, score thresholds, and chromosome 16 exclusion. Do not assume that an analysis-specific reference exception transfers to a new study.

Record the source/deployment commits, Nextflow version, container tags (and resolved digests where available), input dataset/reference versions, submitted parameters, final resource configuration, cohort definition, mask-builder logs/audits, and Step 1 model used for each Step 2 run. Keep controlled-data identifiers in the authorised environment; task logs and mask lists may contain them.

This README describes the checked-in implementation and existing stage documentation. It does not assert that a new CloudOS run has been executed or validated. Software is distributed under the repository's [MIT license](LICENSE); access to GEL datasets is separate.

<a id="contact"></a>

## 💌 Contact

Questions, feedback, or ideas? Contact **Shelley (Shiyu Zhang)** at [shiyuzhang0522@gmail.com](mailto:shiyuzhang0522@gmail.com), or [open a GitHub issue](https://github.com/shiyuzhang0522/Genomics.England.CloudOS.GWAS.REGENIE.Pipeline/issues/new).
