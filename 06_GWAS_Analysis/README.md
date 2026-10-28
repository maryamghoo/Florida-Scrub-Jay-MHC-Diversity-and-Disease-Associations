# GWAS phenotype preparation, association testing, annotation, and visualization

The R script identifies individuals shared between the whole-genome and MHC datasets and creates incidence and intensity phenotypes. **Its three files in `r_output/` are used directly as inputs to the Bash script.** No `bash_input/` folder or copying is required. The optional third R script creates a two-panel Manhattan plot from the completed PLINK results. For plotting, the four PLINK association files must first be copied into `r_input/` alongside `chr_map.txt`.

## Requirements

- R packages: `dplyr`, `readr` (phenotype preparation); `data.table`, `ggplot2`, `ggrepel` (Manhattan plot).
- Bash, `bcftools`, PLINK 1.9, Python 3, Perl, and `awk`/`grep`.
- ANNOVAR, the species-specific `fsj` `refGene` database, and the corresponding annotated GFF3 for functional annotation.

## Folder structure and files

```text
05_GWAS_Analysis/
├── 01_GWAS_Phenotype_Preparation.R
├── 02_GWAS_Genotype_Association_Annotation.sh
├── 03_GWAS_Manhattan_Plot.R
├── README.md
├── r_input/
│   ├── pop140.txt
│   ├── MHC_sample_alleles.csv
│   ├── chr_map.txt                 # accession -> chromosome/scaffold label
│   ├── snps_afterLD_DP3_unrelated.vcf.gz
│   ├── incidence_gwas.assoc.logistic           # copied for plotting
│   ├── incidence_gwas.assoc.logistic.adjusted  # copied for plotting
│   ├── intensity_gwas.assoc.linear             # copied for plotting
│   └── intensity_gwas.assoc.linear.adjusted    # copied for plotting
├── r_output/
│   ├── shared_sampleIDs.txt
│   ├── mhc_incidence.csv
│   └── mhc_intensity.csv
├── annotation_resources/
│   ├── table_annovar.pl
│   ├── annovar_db/                # fsj refGene database files
│   └── Merged_Annotation_FSJgenome_July2024_FINAL.braker.interpro.gff3
├── bash_output_phenotypes/
├── bash_output_vcf/
├── bash_output_genotype/
├── bash_output_incidence/
├── bash_output_intensity/
├── bash_output_annotation/
└── manhattan_output/              # PNG created by script 03
```

The `bash_output_*` folders are created automatically by the Bash script; `manhattan_output/` is created by the plotting R script. You may keep ANNOVAR and the GFF3 elsewhere instead of creating `annotation_resources/` by setting `ANNOVAR`, `ANNOVAR_DB`, and `GFF3` when running the script.

## Inputs

| File or resource | Folder | Used for |
| --- | --- | --- |
| `pop140.txt` | `r_input/` | Whole-genome sample IDs used by the R script |
| `MHC_sample_alleles.csv` | `r_input/` | MHC metadata and disease phenotypes used by the R script |
| `snps_afterLD_DP3_unrelated.vcf.gz` | `r_input/` | Starting SNP VCF for the Bash script |
| `chr_map.txt` | `r_input/` | Two whitespace-separated columns: reference accession and display name (`Chr1`, `HiC_scaffold_182`, etc.); used by the Manhattan plot |
| `incidence_gwas.assoc.logistic` | `r_input/` (copy from `bash_output_incidence/`) | Raw incidence association p-values and SNP positions for the Manhattan plot |
| `incidence_gwas.assoc.logistic.adjusted` | `r_input/` (copy from `bash_output_incidence/`) | BH-FDR adjusted incidence p-values for the Manhattan plot |
| `intensity_gwas.assoc.linear` | `r_input/` (copy from `bash_output_intensity/`) | Raw intensity association p-values and SNP positions for the Manhattan plot |
| `intensity_gwas.assoc.linear.adjusted` | `r_input/` (copy from `bash_output_intensity/`) | BH-FDR adjusted intensity p-values for the Manhattan plot |
| `shared_sampleIDs.txt` | `r_output/` | Shared IDs; used by Bash to subset the VCF |
| `mhc_incidence.csv` | `r_output/` | `SampleID`, `Incidence`; used by Bash to create PLINK phenotype |
| `mhc_intensity.csv` | `r_output/` | `SampleID`, `logIntensity`; used by Bash to create PLINK phenotype |
| `table_annovar.pl` | `annotation_resources/` (or `ANNOVAR`) | ANNOVAR software |
| `annovar_db/` | `annotation_resources/` (or `ANNOVAR_DB`) | `fsj` reference gene annotation database |
| `Merged_Annotation_FSJgenome_July2024_FINAL.braker.interpro.gff3` | `annotation_resources/` (or `GFF3`) | Gene feature and functional annotations |

## Outputs by folder

**`r_output/` — created by `01_GWAS_Phenotype_Preparation.R`**

| File | Contents |
| --- | --- |
| `shared_sampleIDs.txt` | IDs shared between whole-genome and MHC data; one per line |
| `mhc_incidence.csv` | Shared samples with nonmissing infection incidence (0/1) |
| `mhc_intensity.csv` | Shared infected samples with nonmissing log-transformed intensity |

**`bash_output_phenotypes/` — PLINK-ready phenotypes**

| File | Contents |
| --- | --- |
| `incidence_pheno.txt` | Headerless `FID IID Incidence` |
| `intensity_pheno.txt` | Headerless `FID IID logIntensity` |

**`bash_output_vcf/` — shared-sample VCF**

| File | Contents |
| --- | --- |
| `MHC_snps.vcf.gz` | SNP VCF restricted to shared samples |
| `MHC_snps.vcf.gz.csi` | VCF index (default bcftools index format) |
| `MHC_snps.stats` | Variant/sample statistics from bcftools |

**`bash_output_genotype/` — initial PLINK data**

| Files | Contents |
| --- | --- |
| `genotype_base.bed`, `genotype_base.bim`, `genotype_base.fam` | Converted binary SNP genotypes, variant information, and sample information |

**`bash_output_incidence/` — incidence QC, PCA, and GWAS**

| Files | Contents |
| --- | --- |
| `incidence_qc.bed`, `.bim`, `.fam` | SNP dataset after `--geno 0.20 --mac 5` filtering |
| `incidence_pca.eigenvec`, `.eigenval` | 10 principal components and their eigenvalues |
| `incidence_pcs.txt` | FID, IID, PC1, and PC2 covariates |
| `incidence_gwas.assoc.logistic` | SNP logistic regression estimates and p-values |
| `incidence_gwas.assoc.logistic.adjusted` | Multiple-testing-adjusted p-values, including FDR_BH |
| `incidence_FDR_significant.txt` | Incidence SNPs with FDR_BH < 0.05 |

**`bash_output_intensity/` — intensity QC, PCA, and GWAS**

| Files | Contents |
| --- | --- |
| `intensity_keep.txt` | FID and IID for individuals with intensity phenotypes |
| `intensity_qc.bed`, `.bim`, `.fam` | Infected-bird SNP dataset after `--geno 0.20 --mac 5` filtering |
| `intensity_pca.eigenvec`, `.eigenval` | 10 principal components and their eigenvalues |
| `intensity_pcs.txt` | FID, IID, PC1, and PC2 covariates |
| `intensity_gwas.assoc.linear` | SNP linear regression estimates (BETA) and p-values |
| `intensity_gwas.assoc.linear.adjusted` | Multiple-testing-adjusted p-values, including FDR_BH |
| `intensity_FDR_significant.txt` | Intensity SNPs with FDR_BH < 0.05 |

**`bash_output_annotation/` — significant intensity SNPs and gene annotation**

| File | Contents |
| --- | --- |
| `significant_intensity_snps.txt` | FDR-significant intensity SNP IDs |
| `significant_intensity_effects.txt` | Regression results (including BETA) for those SNPs |
| `significant_intensity_adjusted.txt` | Adjusted association results for those SNPs |
| `significant_intensity_snps.vcf` | Significant variants extracted from the VCF |
| `significant_intensity_snps.avinput` | ANNOVAR five-column variant input |
| `significant_intensity_snps_IDmap.txt` | SNP ID, chromosome, position, REF and ALT |
| `significant_intensity_snps_annotated.fsj_multianno.txt` | ANNOVAR gene-based annotation (when significant SNPs exist) |
| `significant_gene_ids.txt` | Gene IDs extracted from ANNOVAR |
| `significant_genes_full.gff3` | All matching GFF3 records |
| `significant_genes_only.gff3` | Only matching gene-level GFF3 records |
| `significant_genes_functional_annotation_FULL.tsv` | Detailed annotations, one row per selected GFF3 feature |
| `significant_genes_functional_annotation_SUMMARY.tsv` | Summarized annotations, one row per selected gene |

**`manhattan_output/` — intended output folder for `03_GWAS_Manhattan_Plot.R`**

| File | Contents |
| --- | --- |
| `Plasmodium_GWAS_FINAL.png` | Two-panel incidence/intensity Manhattan plot with named chromosomes, the MHC-containing scaffold, and significant SNP position labels |

The plotting script reads **all five files directly from `r_input/`**: the four copied raw/adjusted PLINK results and `chr_map.txt`. The Bash workflow still writes the original association results to `bash_output_incidence/` and `bash_output_intensity/`; copying them for plotting does not change the original analyses. It plots raw SNP p-values as `-log10(P)` and uses the existing `FDR_BH` values for significance; it **does not recalculate FDR after omitting scaffolds from the display**. It expects the MHC-containing scaffold `JBGGOT010000146.1` to appear in `chr_map.txt` as `HiC_scaffold_182` and to have retained GWAS SNPs. The R plot marks this entire **scaffold**, not verified nucleotide coordinates of the MHC locus.

PLINK also produces `.log` files and sometimes `.nosex` files in the relevant stage folders. ANNOVAR can produce other intermediate files. The annotation TSV files have headers even when no genes are identified. Functional information is included only if present in the GFF3.

## Processing and run order

1. Run the R script from the project directory:

   ```bash
   Rscript 01_GWAS_Phenotype_Preparation.R
   ```

2. Keep its three outputs in **`r_output/`**. The Bash script reads them there automatically; do not copy or move them.
3. Place the starting VCF in **`r_input/`**. Make ANNOVAR resources available in **`annotation_resources/`** or configure their paths.
4. Run the Bash script (from any working directory):

   ```bash
   bash /path/to/05_GWAS_Analysis/02_GWAS_Genotype_Association_Annotation.sh
   ```

If software or annotation data live elsewhere, configure paths without editing the Bash script:

```bash
PLINK=/path/to/plink \
BCFTOOLS=/path/to/bcftools \
ANNOVAR=/path/to/table_annovar.pl \
ANNOVAR_DB=/path/to/annovar_database \
GFF3=/path/to/genome_annotation.gff3 \
bash 02_GWAS_Genotype_Association_Annotation.sh
```

5. After the Bash script finishes, copy its four GWAS association files to **`r_input/`** (from the `05_GWAS_Analysis/` project directory):

   ```bash
   cp bash_output_incidence/incidence_gwas.assoc.logistic r_input/
   cp bash_output_incidence/incidence_gwas.assoc.logistic.adjusted r_input/
   cp bash_output_intensity/intensity_gwas.assoc.linear r_input/
   cp bash_output_intensity/intensity_gwas.assoc.linear.adjusted r_input/
   ```

   Ensure `r_input/chr_map.txt` is also present. These **five files** are the inputs for `03_GWAS_Manhattan_Plot.R`; the plotting script no longer reads the four association files directly from the `bash_output_*` directories.

6. Run the Manhattan plotting script **from the project directory** (`05_GWAS_Analysis/`):

   ```bash
   Rscript 03_GWAS_Manhattan_Plot.R
   ```

   The script uses `data_dir <- "r_input"`, which is relative to the current R working directory. Run it from the project directory (or set your R working directory there first). The PNG is written to the output directory specified by your local script; if you use `manhattan_output/`, its documented filename is `Plasmodium_GWAS_FINAL.png`.

   Install the plotting packages if needed:

   ```r
   install.packages(c("data.table", "ggplot2", "ggrepel"))
   ```

## Manhattan plot design and interpretation

- **Panels:** (A) infection incidence; (B) log-transformed infection intensity among infected birds.
- **Included sequences:** named chromosomes (`Chr1`, `Chr1A`, etc., and `ChrZ`) plus the MHC-containing `HiC_scaffold_182`, placed last for visualization. All other unplaced scaffolds are omitted **only from the figure**, not from association testing or the original multiple-testing correction.
- **Colors and labels:** alternating chromosome colors; red points mark SNPs with `FDR_BH < 0.05`, labeled by chromosome and position in Mb. Light-blue shading and a dashed blue line identify the MHC-containing scaffold.
- **X-axis scaling:** chromosome/scaffold widths are based on the maximum SNP coordinate observed in the plotted dataset, with a minimum *display* width equal to 12% of the largest observed span. SNPs retain within-sequence relative position, but inter-sequence display distances are **not** proportional to true physical genomic distances. This is a figure-design choice to keep short chromosome labels readable.
- **MHC limitation:** the blue marker identifies the whole **unplaced scaffold**, not the exact MHC interval. Positioning it after the named chromosomes does not establish its physical distance from other associated SNPs. BLAST coordinates derived from a different accession should not be transferred without verifying the assembly correspondence.
- **Genes vs. SNPs:** the plot labels significant SNP associations, not gene names. The functional annotation generated in `bash_output_annotation/` is a separate step.

Incidence is tested using additive logistic regression with PLINK's `--1` option for 0/1 phenotypes. Intensity uses additive linear regression of log-transformed parasite intensity. The analyses use PC1 and PC2, perform separate genotype filtering and PCA, and select Benjamini–Hochberg FDR-adjusted p-values below 0.05. Annotation is performed on significant intensity SNPs.
