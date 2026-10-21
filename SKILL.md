---
name: tassel-cli
description: Construct, validate, execute, and troubleshoot TASSEL 5 command-line pipelines (run_pipeline.pl / run_pipeline.bat) for bioinformatics workflows including genotype format conversion, site/taxa filtering, GWAS (GLM, MLM, FastAssociation), kinship, PCA, MDS, linkage disequilibrium, imputation, and the GBSv2 SNP-calling pipeline. Use when a user or bioinformaticist wants to analyze genetic or phenotypic data with TASSEL from the command line.
---

# TASSEL 5 Command-Line Bioinformatics Skill

TASSEL 5 (Trait Analysis by aSSociation, Evolution and Linkage) provides a command-line pipeline interface (`run_pipeline.pl` on macOS/Linux, `run_pipeline.bat` on Windows) for reproducible genomics workflows.

- **Getting Started**: https://tassel.maizegenetics.net/getting_started/
- **Pipeline Tutorial**: https://tassel.maizegenetics.net/pipelines/pipeline-tutorial/
- **Pipeline CLI Reference**: https://tassel.maizegenetics.net/pipelines/tassel5-pipeline-cli/
- **GBSv2 Pipeline Guide**: https://tassel.maizegenetics.net/gbsv2_pipeline/

## Agent Workflow

When a user requests a TASSEL analysis:
1. Understand the objective: Identify the requested analysis, input formats, expected outputs, and computational constraints.
2. Inspect the environment: Locate the TASSEL installation, determine the operating system, and verify the Java version.
3. Inspect the input data: Check file existence, headers, taxa identifiers, chromosome naming, and other relevant metadata.
4. Select the appropriate plugins: Prefer TASSEL 5 self-describing plugins over legacy options where suitable.
5. Verify plugin parameters: Use -<PluginName> -help rather than assuming parameter names or defaults.
6. Construct the pipeline: Use correct plugin ordering, forks, joins, and export flags.
7. Validate before execution: Check all input/output paths, parameters, and expected intermediate datasets.
8. Execute when authorized: Run the analysis only when the user has requested execution and sufficient resources are available.
9. Check the results: Examine execution logs, output files, and relevant quality metrics.
10. Explain the analysis: Provide the exact command, parameter choices, assumptions, and output descriptions.

Never claim a pipeline has executed successfully unless execution and its outputs have been verified.

## Preflight Validation and Execution Safeguards

Always distinguish between **generating a command** for the user and **executing an analysis** on their system.

### Preflight Data Checks (Before Running)

- **Taxa identifier matching**: Verify that sample/taxa IDs match exactly (including case and delimiters) across genotype, phenotype, covariate (`Q`), and kinship (`K`) files.
- **Duplicate sample identifiers**: Check for duplicate taxa in phenotype/covariate tables or unintended replicate merging in GBSv2 key files (`FullSampleName`).
- **Chromosome and position conventions**: Confirm that chromosome names (for example `1` vs `chr1`) and physical coordinate assemblies match across genotypes, BED files, and reference FASTAs, and that VCF/HapMap positions are sorted within each chromosome.
- **Output path safety**: Confirm that `-export` or `-o` paths will not unintentionally overwrite existing raw data or prior results.

### Execution Safeguards

- Estimate memory and runtime requirements where possible, without presenting guesses as guarantees.
- For large GBSv2 or GWAS analyses, explain the resource implications before starting.
- Preserve original genotype and phenotype files.
- Avoid silently changing analysis thresholds or dropping taxa.
- Capture command-line output and TASSEL logs for reproducibility.

## 1. Setup and Execution Basics

- **Prerequisites**:
  - **Java 21** (64-bit JDK or JRE on `PATH`). Verify before running:
    ```bash
    java -version
    ```
  - **Perl**: Required on macOS/Linux to run the `./run_pipeline.pl` and `./start_tassel.pl` wrapper scripts (`perl -v`). On Windows, use `run_pipeline.bat`.
  - **Optional native tools**: OpenBLAS for accelerated matrix operations; BWA or Bowtie2 for GBSv2 read alignment.
- **Memory (Heap Size)**: Always place `-Xms` (initial) and `-Xmx` (maximum) flags **first**, before any pipeline flags:
  ```bash
  ./run_pipeline.pl -Xms2g -Xmx16g -importGuess genotypes.vcf ...
  ```
- **Logging and Memory Diagnostics**:
  - `-debug [logfile.txt]` or `-log [logfile.txt]` writes detailed execution logs.
  - `-printMemoryUsage` can be placed anywhere in a pipeline to print current JVM memory consumption.

## 2. Discovering Any Plugin and Its Parameters

Never guess plugin parameter names. Query TASSEL directly from the CLI:

```bash
# Enumerate all self-describing plugins
./run_pipeline.pl -ListPlugins

# Enumerate self-describing plugins with full parameter descriptions
./run_pipeline.pl -ListPlugins -usage true

# Print parameter help for a specific self-describing plugin
./run_pipeline.pl -KinshipPlugin -help
./run_pipeline.pl -FilterSiteBuilderPlugin -help
./run_pipeline.pl -FixedEffectLMPlugin -help
```

Note: `-ListPlugins` enumerates **self-describing plugins** only. Legacy shorthand flags (such as `-importGuess`, `-vcf`, `-h`, `-r`, `-q`, `-k`, `-intersect`, `-mlm`, `-ld`, and `-export`) are documented separately in the [Pipeline CLI Reference](https://tassel.maizegenetics.net/pipelines/tassel5-pipeline-cli/).

## 3. Core Pipeline Syntax: Plugins, Forks, and Joins

TASSEL pipelines chain steps left-to-right, where each step's output feeds the next step's input:

1. **Self-describing plugins (preferred)**: Start with `-<PluginName>`, pass parameter flags, and close with `-endPlugin`.
   ```bash
   -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05 -siteMinCount 50 -endPlugin
   ```
2. **Legacy convenience flags**: Do not use `-endPlugin` (for example `-importGuess`, `-vcf`, `-h`, `-r`, `-intersect`, `-union`, `-export`, `-exportType`).

### Multi-Branch Pipelines (`-fork`, `-combine`, `-input`)

When an analysis requires multiple inputs (such as genotypes + phenotypes + kinship), load each in its own `-fork<id>` (no space between `-fork` and `<id>`) and merge them with `-combine<id>`:

- `-fork1 ...`: Starts branch 1.
- `-combine3 -input1 -input2 -intersect`: Starts branch 3, waits for branches 1 and 2, and joins them by common taxa (`-intersect`) or all taxa (`-union`).
- `-input<id>`: Passes the output of branch `<id>` to the immediately preceding plugin or `-combine` flag.
- Do **not** append `-runfork1` at the end of commands; it is deprecated and TASSEL runs active forks automatically.

## 4. Input File Formats and Structure (`mdp_*` Examples)

Inspect input headers before building a pipeline to verify compatibility with TASSEL.

### Genotypes (HapMap, VCF, HDF5, PLINK)

- **HapMap (`.hmp.txt` / `.hmp.txt.gz`)**: Tab-delimited; requires 11 fixed metadata columns (`rs#`, `alleles`, `chrom`, `pos`, `strand`, `assembly#`, `center`, `protLSID`, `assayLSID`, `panel`, `QCcode`) followed by taxa columns. Positions must be sorted within each chromosome. Missing calls are `N` or `NN`.
  ```text
  rs#	alleles	chrom	pos	strand	assembly#	center	protLSID	assayLSID	panel	QCcode	33-16	38-11	4226
  PZB00859.1	A/C	1	157104	+	AGPv1	Panzea	NA	NA	maize282	NA	CC	CC	CC
  PZA01271.1	C/G	1	1947984	+	AGPv1	Panzea	NA	NA	maize282	NA	CC	GG	CC
  ```
- **VCF (`.vcf` / `.vcf.gz`)**: Standard VCF v4.x (`#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <Sample1> ...`). Must be sorted by chromosome and position (`-sortPositions` on import). Missing genotypes are `./.`.
- **PLINK (`.ped` + `.map`)**: Tab-delimited 6-column PED header (`FamilyID IndividualID PaternalID MaternalID Sex Phenotype` + genotypes) and 4-column MAP (`Chrom rs# GeneticDist Pos`). TASSEL converts PLINK missing `0` to `N`.
- **HDF5 (`.hmp.h5` / `.h5`)**: TASSEL's high-speed binary genotype format. Depths are stored as log-approximated signed bytes (`-128` to `127`).

### Phenotypes and Covariates (`-r` / `-t` / `-q`)

- **TASSEL 5 `<Phenotype>` format (preferred)**: 3 header lines specifying `<Phenotype>`, column attribute types (`taxa`, `data`, `factor`, `covariate`), and column names. Missing values use `NA`, `NaN`, or `-999`.
  ```text
  <Phenotype>
  taxa	factor	data	data	covariate	covariate
  Taxa	rep	EarHT	dpoll	Q1	Q2
  33-16	1	64.75	64.5	0.014	0.972
  38-11	1	92.25	68.5	0.003	0.993
  ```
- **Legacy `<Trait>` and `<Covariate>` formats** (used in `mdp_traits.txt` and `mdp_population_structure.txt`):
  ```text
  <Trait>	EarHT	dpoll	EarDia
  33-16	64.75	64.5	-999
  38-11	92.25	68.5	37.897
  ```
  ```text
  <Covariate>
  <Trait>	Q1	Q2	Q3
  33-16	0.014	0.972	0.014
  38-11	0.003	0.993	0.004
  ```

### Kinship / Square Matrices (`-k`)

Optional `##` metadata lines, followed by a line with the integer taxon count `n`, followed by `n` tab-delimited rows of `TaxonName\tval_1\t...\tval_n`:
```text
##Matrix_Type=Centered_IBS
281
33-16	1.7978	0.0396	0.0614
38-11	0.0396	1.9099	0.0198
4226	0.0614	0.0198	1.9277
```

### GBSv2 Key File and FASTQ Naming

- **Key file**: Tab-delimited with 4 required case-sensitive headers (`Flowcell`, `Lane`, `Barcode`, `FullSampleName`). Identical `FullSampleName` rows across lanes/flowcells are merged into one taxon:
  ```text
  Flowcell	Lane	Barcode	DNASample	Enzyme	FullSampleName
  C05F2ACXX	5	CTCC	B73	ApeKI	B73:C05F2ACXX:5:250020986
  C05F2ACXX	5	TTCTC	Mo17	ApeKI	Mo17:C05F2ACXX:5:250020994
  ```
- **FASTQ filenames**: Must contain 3, 4, or 5 underscore-delimited tokens matching `Flowcell` and `Lane` in the key file: `flowcell_lane_fastq.txt.gz`, `flowcell_s_lane_fastq.txt.gz`, or `code_flowcell_s_lane_fastq.txt.gz` (accepted extensions include `.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`, `.fastq.txt.gz`, `_sequence.txt.gz`).

### Import and Export Flags

- **Import**: `-importGuess <file>`, `-vcf <file>`, `-h <file>`, `-h5 <file>`, `-plink -ped <f.ped> -map <f.map>`, `-r <file>`, `-q <file>`, `-k <file>`, `-fasta <file>`, `-table <file>`, `-sortPositions`.
- **Export**: `-export <filename> -exportType <Type>` where `<Type>` is `Hapmap`, `HapmapDiploid`, `VCF`, `HDF5`, `Plink`, `ReferenceProbability`, `Depth`, `Fasta`, `Phylip_Seq`, `Phylip_Inter`, `SqrMatrix`, `Phenotype`, `Table`, or `Text`.

## 5. Filtering Genotypes and Taxa

Filtering thresholds (MAF, missingness, heterozygosity) are **study- and biology-dependent**, not universal rules. For example, inbred diversity panels tolerate strict max-heterozygosity filters (such as `0.05`–`0.10`), whereas F1/outcrossing populations require permissive heterozygosity bounds; rare-variant studies use lower MAF cutoffs than small-N GWAS panels. Always state the chosen thresholds and rationale to the user.

```bash
./run_pipeline.pl -Xmx8g -importGuess input.vcf \
  -FilterSiteBuilderPlugin \
    -siteMinCount 80 \
    -siteMinAlleleFreq 0.05 \
    -siteMaxAlleleFreq 0.95 \
    -minHeterozygous 0.0 \
    -maxHeterozygous 0.1 \
    -removeMinorSNPStates true \
  -endPlugin \
  -FilterTaxaBuilderPlugin \
    -minNotMissing 0.8 \
    -maxHeterozygous 0.1 \
  -endPlugin \
  -export filtered.vcf -exportType VCF
```

- Region filtering on `FilterSiteBuilderPlugin`: `-bedFile regions.bed` or `-siteRangeFilterType POSITIONS -startChr 1 -startPos 100000 -endChr 1 -endPos 5000000`.
- Quick subsetting flags: `-includeTaxaInFile keep.txt`, `-excludeTaxaInFile drop.txt`, `-includeSiteNamesInFile sites.txt`, `-subsetSites 10000`, `-subsetSites 0.1 -step`, `-separate 1,2`, `-homozygous`.

## 6. Population Structure, Kinship, and LD

- **Kinship (`-KinshipPlugin`)**: Methods are `Centered_IBS` (default, Endelman-Jannink), `Normalized_IBS` (GCTA), `Dominance_Centered_IBS`, `Dominance_Normalized_IBS`.
- **PCA (`-PrincipalComponentsPlugin`)**: Converts genotypes to numeric (`0`, `0.5`, `1`) and mean-imputes missing values. Exports 3 files: PCs, eigenvalues (variance explained), and eigenvectors.
- **MDS (`-MDSPlugin`)**: Runs classical multidimensional scaling on a `-DistanceMatrixPlugin` output. Preferred over PCA when missingness is high because pairwise distances use only shared non-missing sites.
- **LD (`-ld`)**: `-ld -ldType SlidingWindow -ldWinSize 50 -ldHetTreatment Homozygous -td tab ld_out.txt -ldd png -ldplotsize 1000 -o ld_plot.png`.

## 7. Association Mapping (GWAS: GLM, MLM, FastAssociation)

### Key Scientific and Pipeline Rules for GWAS

- **Why `-intersect` vs `-union` matters**: `-intersect` keeps only taxa present in **both** the genotype and phenotype/covariate datasets (and aligns them for modeling). If taxon IDs are formatted differently between files, `-intersect` silently drops non-matching samples or produces an empty dataset. Always check the retained sample count (`Obs` / `df`) in the output.
- **Why `-excludeLastTrait` is used on `Q` matrices**: Admixture/STRUCTURE subpopulation membership proportions (`Q1 + Q2 + ... + Qk = 1.0`) are linearly dependent with the model intercept. Passing `-excludeLastTrait` after `-q pop_structure.txt` drops the final `Q` column to prevent design-matrix collinearity in GLM and MLM.
- **FastAssociation (`-EqtlAssociationPlugin`)**: Requires **no missing phenotype values** and an `-intersect` join.

## 8. Genotype Imputation: Nucleotide vs. Numerical

TASSEL provides two distinct categories of imputation with different output types and downstream compatibility:

1. **Nucleotide Genotype Imputation** (`-LDKNNiImputationPlugin`, `-FILLINImputationPlugin`, `-FSFHapImputationPlugin`):
   - **Output**: Standard `GenotypeTable` with nucleotide calls (`A`, `C`, `G`, `T`, IUPAC codes), exportable to `VCF`, `Hapmap`, or `HDF5`.
   - **Downstream use**: Compatible with **all** TASSEL plugins (Kinship, LD, PCA, GLM, MLM) and external tools.
   ```bash
   ./run_pipeline.pl -Xmx12g -importGuess masked.hmp.txt.gz \
     -LDKNNiImputationPlugin -highLDSSites 30 -knnTaxa 10 -maxLDDistance 10000000 -endPlugin \
     -export imputed.hmp.txt.gz -exportType Hapmap
   ```
2. **Numerical Genotype Imputation** (`-NumericalGenotypePlugin` followed by `-ImputationPlugin -ByMean true` or kNN):
   - **Output**: `ReferenceProbability` numeric matrix (values in `[0.0, 1.0]` representing major allele dosage), exportable via `-exportType ReferenceProbability`.
   - **Downstream use**: Can be used in PCA, Kinship, GLM (`FixedEffectLMPlugin`), MLM, and `EqtlAssociationPlugin` (with `-genotypeComponent ReferenceProbability`), but **cannot** be exported back to nucleotide VCF/HapMap alleles.

## 9. GBSv2 Discovery and Production Pipeline

GBSv2 stores tags, alignments, and SNPs in a SQLite database (`.db`). Do **not** include a trailing slash on the `-i` FASTQ directory path. Note that `TagExportToFastqPlugin` always writes **FASTQ** format (even if a legacy example names the output `.fa.gz`, prefer `.fq.gz`).

```bash
# 1. Extract barcoded tags from FASTQ directory into SQLite DB
./run_pipeline.pl -Xms16g -Xmx32g -GBSSeqToTagDBPlugin \
  -e ApeKI -i fastq_dir -db gbsv2.db -k keyfile.txt \
  -kmerLength 64 -minKmerL 20 -mnQS 20 -c 10 -mxKmerNum 50000000 -endPlugin

# 2. Export distinct tags in FASTQ format for external alignment
./run_pipeline.pl -Xmx8g -TagExportToFastqPlugin \
  -db gbsv2.db -o tagsForAlign.fq.gz -c 1 -endPlugin

# 3. External aligner (BWA or Bowtie2) -> tags.sam
# bowtie2 -M 4 --very-sensitive-local -x ref_index -U tagsForAlign.fq.gz -S tags.sam

# 4. Load SAM alignments into SQLite DB
./run_pipeline.pl -Xmx8g -SAMToGBSdbPlugin \
  -i tags.sam -db gbsv2.db -aProp 0.0 -aLen 0 -minMAPQ 0 -endPlugin

# 5. Discover SNPs across aligned tags in DB
./run_pipeline.pl -Xmx16g -DiscoverySNPCallerPluginV2 \
  -db gbsv2.db -sC 1 -eC 10 -mnLCov 0.1 -mnMAF 0.01 -ref ref.fa -endPlugin

# 6. Profile SNP quality statistics
./run_pipeline.pl -Xmx8g -SNPQualityProfilerPlugin \
  -db gbsv2.db -statFile snp_quality_stats.txt -endPlugin

# 7. Call genotypes directly to VCF (or .h5) in Production
./run_pipeline.pl -Xms16g -Xmx32g -ProductionSNPCallerPluginV2 \
  -db gbsv2.db -e ApeKI -i fastq_dir -k keyfile.txt \
  -kmerLength 64 -o production_snps.vcf -endPlugin
```

## 10. End-to-End Reproducible Workflows (`mdp_*` Tutorial Data)

### Workflow 1: VCF Import -> Filter Sites and Taxa -> Export Filtered VCF

- **Inputs**: `mdp_genotype_w_depth.vcf` (281 taxa, 3,093 sites)
- **Command**:
  ```bash
  ./run_pipeline.pl -Xmx8g \
    -vcf dataFiles/Tutorial/mdp_genotype_w_depth.vcf -sortPositions \
    -FilterSiteBuilderPlugin -siteMinCount 200 -siteMinAlleleFreq 0.05 -maxHeterozygous 0.10 -endPlugin \
    -FilterTaxaBuilderPlugin -minNotMissing 0.80 -endPlugin \
    -genotypeSummary overall,site,taxa \
    -export filtered_mdp.vcf -exportType VCF
  ```
- **Outputs & Verification**: `filtered_mdp.vcf` and summary tables. Verify in the log or summary output that the remaining site count and taxon count match expectations and that `filtered_mdp.vcf` is non-empty.

### Workflow 2: HapMap -> PCA -> Export PCs and Variance Explained

- **Inputs**: `mdp_genotype.hmp.txt`
- **Command**:
  ```bash
  ./run_pipeline.pl -Xmx8g \
    -h dataFiles/Tutorial/mdp_genotype.hmp.txt \
    -FilterSiteBuilderPlugin -siteMinAlleleFreq 0.05 -endPlugin \
    -PrincipalComponentsPlugin -covariance true -ncomponents 5 -endPlugin \
    -export mdp_pca
  ```
- **Outputs & Verification**: Produces `mdp_pca1.txt` (PC scores per taxon), `mdp_pca2.txt` (eigenvalues / proportion of variance explained), and `mdp_pca3.txt` (eigenvectors / site loadings). Confirm taxon count in `mdp_pca1.txt` and inspect cumulative variance in `mdp_pca2.txt`.

### Workflow 3: Genotype + Phenotype + Q -> Kinship -> MLM -> Export Association Results

- **Inputs**: `mdp_genotype.hmp.txt`, `mdp_traits.txt`, `mdp_population_structure.txt`
- **Command**:
  ```bash
  ./run_pipeline.pl -Xmx10g \
    -fork1 -h dataFiles/Tutorial/mdp_genotype.hmp.txt \
      -FilterSiteBuilderPlugin -siteMinCount 150 -siteMinAlleleFreq 0.05 -endPlugin \
    -fork2 -r dataFiles/Tutorial/mdp_traits.txt \
    -fork3 -q dataFiles/Tutorial/mdp_population_structure.txt -excludeLastTrait \
    -fork4 -KinshipPlugin -method Centered_IBS -endPlugin -input1 \
    -combine5 -input1 -input2 -input3 -intersect \
    -mlm -mlmVarCompEst P3D -mlmCompressionLevel Optimum \
      -input4 -mlmMaxP 1.0 -mlmOutputFile mdp_mlm_out
  ```
- **Outputs & Verification**: Produces residuals (`mdp_mlm_out1.txt`..), MLM statistics (`mdp_mlm_out*_stats.txt`), allele effects, and compression table. Verify that all 3 traits (`EarHT`, `dpoll`, `EarDia`) were modeled, check sample size `df`/`Obs`, and inspect p-value distribution for inflation.

### Workflow 4: GBSv2 Raw Reads -> Tag DB -> Alignment -> SNP Discovery -> VCF

- **Inputs**: `dataFiles/GBS/Chr9_10-20000000` (FASTQ directory), `Pipeline_Testing_key.txt`, `ZmB73_RefGen_v2_chr9_10_1st20MB.fasta`
- **Intermediate Data**: `gbsv2.db` (SQLite database), `tagsForAlign.fq.gz`, `tags.sam`
- **Verification**: Check `Total number of low quality reads=` and total tag counts in step 1 logs; check alignment rate in Bowtie2/BWA stderr; verify discovered SNP count in `snp_quality_stats.txt` and non-empty genotype calls in `production_snps.vcf`.

## 11. Result Interpretation and Quality-Control Guidance

A zero exit code only indicates that the JVM did not crash; it does **not** prove scientific validity (for example, an imputation run can finish cleanly with high error rates, or an `-intersect` join can drop 90% of samples due to naming mismatches). Always perform these post-run checks:

| Analysis | Suggested validation |
| :--- | :--- |
| Genotype filtering | Original vs. remaining taxa and SNP counts; missingness; MAF distributions |
| GLM / MLM | Sample counts, association statistics, unexpected inflation, model covariates |
| PCA | Variance explained, eigenvalues, missing-data treatment |
| Kinship | Matrix dimensions, taxa identifiers, missing or unexpected values |
| Imputation | Missingness before/after, genotype concordance when validation data exists |
| GBSv2 | Read counts, tag counts, alignment rates, SNP discovery counts |

## 12. Common Pitfalls and Troubleshooting

- **`OutOfMemoryError: Java heap space`**: Increase `-Xmx` at the very start of `./run_pipeline.pl` without altering input data or scientific thresholds. For GBSv2 step 1 (`GBSSeqToTagDBPlugin`), allocate roughly 2.5x `mxKmerNum * 100 bytes`.
- **Positions not sorted in VCF/HapMap**: Add `-sortPositions` immediately after the import flag, or run `-SortGenotypeFilePlugin -inputFile in.vcf -outputFile sorted.vcf -fileType VCF -endPlugin`.
- **`-input<id>` placement error**: `-input<id>` must immediately follow either `-combine<id>` or the `-endPlugin` of a plugin that consumes multiple inputs (like `-mlm` or `-KinshipPlugin`).
- **GBSv2 `GBSSeqToTagDBPlugin` processes 0 files without error**: Ensure the `-i` directory path has no trailing slash, FASTQ filenames match `flowcell_lane_fastq.txt.gz` (or `.fq.gz`/`.fastq.gz`), and `Flowcell` and `Lane` values match the key file.
- **HDF5 depth values appear capped at 127**: TASSEL encodes read depths in HDF5 as signed log-approximated bytes (`-128` to `127`). Export from TASSEL to VCF with `-exportIncludeDepth true` rather than reading raw HDF5 depth bytes externally.

## 13. Acceptance Test Checklist

Before publishing changes to this skill, test a fresh agent session against these prompts:

- [ ] Convert a VCF to HapMap and preserve the original file.
- [ ] Filter SNPs by MAF and missingness, explaining the chosen thresholds.
- [ ] Perform a GWAS using genotypes, phenotypes, and a kinship matrix.
- [ ] Build a GBSv2 workflow from FASTQ files and a reference genome.
- [ ] Diagnose a TASSEL OutOfMemoryError without changing input data.
- [ ] Handle a nonexistent plugin parameter without inventing a replacement.
