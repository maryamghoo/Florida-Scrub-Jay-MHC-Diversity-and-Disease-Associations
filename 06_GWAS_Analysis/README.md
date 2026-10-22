# GWAS phenotype preparation, association testing, and annotation

The R script identifies individuals shared between the whole-genome and MHC datasets and creates incidence and intensity phenotypes. **Its three files in `r_output/` are used directly as inputs to the Bash script.** No `bash_input/` folder or copying is required.

## Requirements

- R packages: `dplyr`, `readr` (for the R phenotype-preparation script).
- Bash, `bcftools`, PLINK 1.9, Python 3, Perl, and `awk`/`grep`.
- ANNOVAR, the species-specific `fsj` `refGene` database, and the corresponding annotated GFF3 for functional annotation.

## Folder structure and files

```text
05_GWAS_Analysis/
├── 01_GWAS_Phenotype_Preparation.R
├── 02_GWAS_Genotype_Association_Annotation.sh
├── README.md
├── r_input/
│   ├── pop140.txt
│   ├── MHC_sample_alleles.csv
│   └── snps_afterLD_DP3_unrelated.vcf.gz
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
└── bash_output_annotation/
```

The `bash_output_*` folders are created automatically by the Bash script. You may keep ANNOVAR and the GFF3 elsewhere instead of creating `annotation_resources/` by setting `ANNOVAR`, `ANNOVAR_DB`, and `GFF3` when running the script.

## Inputs

| File or resource | Folder | Used for |
| --- | --- | --- |
| `pop140.txt` | `r_input/` | Whole-genome sample IDs used by the R script |
| `MHC_sample_alleles.csv` | `r_input/` | MHC metadata and disease phenotypes used by the R script |
| `snps_afterLD_DP3_unrelated.vcf.gz` | `r_input/` | Starting SNP VCF for the Bash script |
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

Incidence is tested using additive logistic regression with PLINK's `--1` option for 0/1 phenotypes. Intensity uses additive linear regression of log-transformed parasite intensity. The analyses use PC1 and PC2, perform separate genotype filtering and PCA, and select Benjamini–Hochberg FDR-adjusted p-values below 0.05. Annotation is performed on significant intensity SNPs.

The Bash script has been syntax-checked, but the complete pipeline has not been run against the actual genomic datasets.
