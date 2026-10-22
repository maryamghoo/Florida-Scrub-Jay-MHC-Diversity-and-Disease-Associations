# MHC genetic cluster allele analysis

This R script summarizes MHC class II allele frequencies and presence across genetic clusters and creates an allele-frequency plot. It excludes the Admixed group.

## Requirements

Install R and these packages, if needed:

```r
install.packages(c("ggplot2", "RColorBrewer"))
```

- `ggplot2`: creates the allele-frequency plot.
- `RColorBrewer`: supplies the plot colors.

## Input data

The input file, `MHC_genetic_cluster_alleles.csv`, comes from the **`01_AmpliSAS_Postprocessing`** directory. It is produced by the AmpliSAS postprocessing R script and contains one row per genetic cluster, the `SampleSize` column, and allele carrier counts for `ApcoIIB001` through `ApcoIIB038`.

Copy the file from the earlier step into the `r_input` folder for this analysis.

```text
project/
├── MHC_Genetic_Cluster_Analysis.R
├── r_input/
│   └── MHC_genetic_cluster_alleles.csv
└── r_output/
    └── genetic_cluster_analysis/     # created by the script
```

The script uses the current working directory as `project/`. Set the working directory to this folder before running it.

## Run

```bash
Rscript MHC_Genetic_Cluster_Analysis.R
```

You can also run the script in RStudio after setting the working directory to `project/`.

## What the script does

1. Reads the genetic-cluster allele counts and excludes `Admixed`.
2. Calculates each allele's relative frequency by dividing its carrier count by the sum of all allele carrier counts in that cluster. **This is not the fraction of sampled individuals carrying the allele** (`count / SampleSize`).
3. Creates an allele-frequency bar plot with the original allele colors, symbols, and genetic-cluster facet order.
4. Counts the number of non-Admixed clusters where each allele is present, and lists those clusters.

All 38 allele columns are included, including the non-functional alleles `ApcoIIB002` and `ApcoIIB022`.

## Output files

The script saves these files in `r_output/genetic_cluster_analysis/`:

| File | Contents |
| --- | --- |
| `MHC_cluster_allele_frequencies.csv` | Sample size, total allele-presence observations, and relative frequency of each allele for each cluster |
| `cluster_histogram.png` | Allele-frequency plot across genetic clusters (300 dpi) |
| `MHC_allele_cluster_presence.csv` | Each allele, number of clusters where it occurs, and cluster names |

The script does **not** use FASTA alignments or calculate nucleotide diversity.
