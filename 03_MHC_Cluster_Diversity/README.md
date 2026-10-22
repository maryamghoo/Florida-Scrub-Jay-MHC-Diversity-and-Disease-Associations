# MHC genetic-cluster nucleotide diversity

This R script calculates MHC nucleotide diversity for genetic clusters and compares it with genome-wide diversity. It does not calculate diversity for individual birds.

## Requirements

R packages: `ape`, `pegas`, `boot`, `ggplot2`, and `patchwork`.

## Input files

Put these inside the working directory:

```text
project/
├── MHC_Cluster_Nucleotide_Diversity_Commented.R
├── r_input/
│   ├── MHC_genetic_cluster_alleles.csv
│   ├── cluster_alignments/
│   │   ├── Cluster_Merritt_Island.fasta
│   │   ├── Cluster_SL.fasta
│   │   └── ...
│   └── genome_pi/
│       ├── merritt.windowed.pi
│       ├── sl.windowed.pi
│       ├── cbc.windowed.pi
│       ├── sbc.windowed.pi
│       ├── archbold.windowed.pi
│       └── ocala.windowed.pi
└── r_output/                  # Created automatically
```

The cluster allele-count CSV and cluster FASTA alignments come from the earlier `01_AmpliSAS_Postprocessing` workflow. Each cluster FASTA must contain aligned functional allele sequences, with the original gap positions preserved. The six `.windowed.pi` files provide genome-wide window-level diversity and must contain `PI` and `N_VARIANTS` columns.

## Run

From the project working directory:

```bash
Rscript MHC_Cluster_Nucleotide_Diversity_Commented.R
```

## Calculations

- **Unweighted MHC diversity:** `pegas::nuc.div()` on distinct functional alleles, each counted once.
- **Carrier-weighted MHC diversity:** pairwise DNA distances weighted by the number of birds carrying each allele. This is carrier-weighted, not allele-copy-frequency-weighted diversity.
- **Other MHC summaries:** number of functional alleles, segregating sites, proportion of segregating sites, 10,000 bootstrap replicates, and BCa intervals where available.
- **Genome comparison:** summary genome-wide π, observed heterozygosity (`Ho`), expected heterozygosity (`He`), and four Spearman correlations. `FIS` and `FROH` are not used.
- **Windowed genome diversity:** retains windows with at least 20 variants for the genome-wide boxplot.

## Outputs

All files are saved to `r_output/cluster_diversity/`:

- `MHC_cluster_diversity.csv` — MHC diversity measures for each cluster.
- `MHC_cluster_diversity.png` and `.pdf` — unweighted versus carrier-weighted diversity points.
- `pi_clusteres.jpg` — unweighted π with bootstrap confidence intervals.
- `MHC_genome_comparison.csv` and `MHC_genome_correlations.csv` — MHC and genome-wide comparisons.
- `MHC_genome_comparison.png` and `.pdf` — four-panel MHC versus genome-wide comparison.
- `MHC_genome_windowed_pi.csv` and `MHC_genome_windowed_pi_summary.csv` — filtered window-level data and summary.
- `MHC_genomewide_pi_vertical.png` and `.pdf` — **vertical two-panel boxplot**, MHC bootstrap π above genome-wide windowed π.

The last three windowed-diversity outputs are generated when all six `.windowed.pi` files are present. 
