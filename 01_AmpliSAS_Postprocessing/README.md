# MHC allele analysis

This R script filters AmpliSAS allele results, exports allele sequences, and prepares individual- and genetic-cluster-level MHC datasets. Nucleotide diversity is calculated using a sequence alignment made separately in Geneious Prime.

## Requirements

Install R and these R packages:

```r
install.packages(c("readxl", "openxlsx", "ape", "pegas"))
```

- `readxl`: reads the AmpliSAS Excel results.
- `openxlsx`: writes the filtered Excel file.
- `ape`: reads the aligned DNA sequences.
- `pegas`: calculates nucleotide diversity.

For the alignment step, use **Geneious Prime v. 2023.1.2**, with **Clustal Omega** and its **default parameters**. Alignment is done manually outside R.

## Folder setup

The script uses these folders:

```text
project/
├── MHC_Allele_Analysis.R
├── r_input/
│   ├── results.xlsx
│   ├── duplicate_samples.txt
│   ├── sampling_info.txt
│   ├── disease_info.csv
│   ├── supertypes.csv
│   └── functional_alleles_alignment.fasta  # added after Geneious alignment
└── r_output/                               # created by the script
```

**Input files**

- `results.xlsx`: AmpliSAS results. The script expects four header/metadata rows, allele data beginning on row 5, the first nine columns as allele information, and sample columns beginning at column 10.
- `duplicate_samples.txt`: sample column names to remove, one per line. Optional: if the file is absent, this step is skipped.
- `sampling_info.txt`: tab-delimited file with `SampleID`, `Region`, `Location`, `Genetic_cluster`, `Y`, `X`, `Year`, `Month`, `Sex`, `Mass`, `Age`, and `Age_rank`.
- `disease_info.csv`: file with `SampleID`, `Incidence`, and `Intensity`.
- `supertypes.csv`: mapping of functional MHC alleles to supertypes, with columns `allele` (e.g., `ApcoIIB001`) and `supertype` (`A`, `B`, or `C`). These assignments are derived from the **supertyping methods** described in the supplementary materials of Ghoojaei et al. (2026), *Immunogenetic diversity, selection, and genome-wide disease associations in the declining Florida scrub-jay (Aphelocoma coerulescens)*.
- `functional_alleles_alignment.fasta`: aligned functional allele sequences exported from Geneious Prime. This file is needed for the nucleotide diversity step, not for the initial FASTA export.

## Processing steps

### 1. Filter the AmpliSAS results

The script:

1. Keeps allele rows with `LENGTH` from **200 through 269**, inclusive, and initial `SAMPLES >= 2`.
2. Removes sample columns with `DEPTH_AMPLICON < 700`.
3. Removes any sample columns named in `duplicate_samples.txt`.
4. Recalculates `DEPTH` and `SAMPLES` from the retained sample columns.
5. Names the retained allele rows `ApcoIIB001`, `ApcoIIB002`, and so on, in their existing order.

It saves `r_output/results_filtered.xlsx`. Alleles are not filtered again based on the recalculated `SAMPLES` values.

### 2. Export the allele sequences

The script creates two FASTA files in `r_output/`:

- `MHC_all_alleles.fasta`: all retained alleles, with `ApcoIIB002` and `ApcoIIB022` labeled **with stop codon** in their FASTA headers.
- `MHC_functional_alleles.fasta`: excludes `ApcoIIB002` and `ApcoIIB022`.

These two alleles are considered non-functional based on the predefined IDs; the code does not independently test the sequences for stop codons.

### 3. Align the functional allele sequences in Geneious Prime

This is a **manual step** between FASTA export and nucleotide diversity calculation.

1. Import `r_output/MHC_functional_alleles.fasta` into **Geneious Prime v. 2023.1.2**.
2. Select all functional allele sequences and run **Clustal Omega**, using **default parameters**.
3. Export the multiple sequence alignment in **aligned FASTA** format, preserving the allele identifiers (for example, `ApcoIIB001`).
4. Save the exported file as `r_input/functional_alleles_alignment.fasta`.

Use the same alignment for all individuals. Do not remove alignment gaps or rename the allele identifiers. If using an existing previously generated alignment, place that file at the same path instead; its sequence content and gap placement will determine the nucleotide diversity values.

### 4. Create the individual dataset, calculate nucleotide diversity, and assign supertypes

The script combines the remaining sample data with `sampling_info.txt` and `disease_info.csv`, matching rows by `SampleID`. For AmpliSAS sample names, it removes the first underscore and everything after it before matching.

For each individual, the output contains:

- Sampling and disease information, including `logIntensity = log10(Intensity)` for positive intensities; zero or missing intensity produces a missing log value.
- Presence (`1`) or absence (`0`) for `ApcoIIB001` through `ApcoIIB038`.
- `Total_Alleles` and `Functional_Alleles` counts.
- `Genotype`: functional alleles present, separated by underscores, excluding `ApcoIIB002` and `ApcoIIB022`.
- `Nucleotide_Diversity`: calculated from the Geneious alignment using `pegas::nuc.div()`.
- `SupertypeIIA`, `SupertypeIIB`, and `SupertypeIIC`: presence (`1`) if at least one functional allele assigned to that supertype is present, otherwise absence (`0`), using `r_input/supertypes.csv`.
- `Supertype_count`: number of distinct supertypes present in the individual (0–3).

For nucleotide diversity, one functional allele gives `0`; no functional alleles or missing allele names in the alignment give `NA`. The function follows the default settings of `pegas::nuc.div()`.

The script assigns supertypes after calculating nucleotide diversity and before saving the individual dataset. The assignments in `supertypes.csv` come from the **supertyping methods in the supplementary materials of Ghoojaei et al. (2026)**, *Immunogenetic diversity, selection, and genome-wide disease associations in the declining Florida scrub-jay (Aphelocoma coerulescens)*. The script uses these predefined assignments; it does not perform supertyping itself. The non-functional alleles `ApcoIIB002` and `ApcoIIB022` are excluded from supertype assignment.

The script saves `r_output/MHC_sample_alleles.csv`.

### 5. Summarize MHC allele diversity by location and genetic cluster

The script also creates two summary tables from the individual allele presence/absence data:

- `r_output/MHC_location_allele_summary.csv`: for each `Location`, reports the sample size, number of distinct MHC alleles, number of distinct functional MHC alleles, average number of alleles per bird, and average number of functional alleles per bird.
- `r_output/MHC_cluster_mean_alleles.csv`: for each `Genetic_cluster`, reports the average number of alleles per bird and average number of functional alleles per bird.

A distinct allele is counted once per location even when multiple birds carry it. Functional counts exclude `ApcoIIB002` and `ApcoIIB022`. Individuals without a recorded location or genetic cluster are excluded from the respective summaries.

### 6. Summarize alleles by genetic cluster

The script saves `r_output/MHC_genetic_cluster_alleles.csv`, with one row per `Genetic_cluster`, a `SampleSize` column counting unique sample IDs, and a column for each of the 38 alleles giving the number of individuals where that allele is present. Samples without a genetic cluster are excluded.

### 7. Export aligned functional alleles for each genetic cluster

After creating the genetic-cluster summary, the additional R section uses the existing `aligned_dna` alignment and the allele presence/absence columns in `result` to generate one FASTA file per genetic cluster.

- Only alleles detected in at least one individual in that cluster are included.
- Only **functional alleles** are exported; `ApcoIIB002` and `ApcoIIB022` are excluded.
- Each distinct allele is written **once per cluster**, even if multiple individuals carry it.
- The alignment and gap positions from `r_input/functional_alleles_alignment.fasta` are preserved. The sequences are not realigned separately for each cluster.
- Clusters without a detected functional allele do not get a FASTA file.

The files are saved to `r_output/cluster_alignments/` and named `Cluster_<Genetic_cluster>.fasta` (with spaces or special characters in cluster names replaced by underscores). For example, `Cluster_1.fasta` contains the aligned functional alleles present in genetic cluster 1.

These files show **which distinct functional allele sequences occur in each cluster**. They do **not** represent allele frequencies or copy numbers; use the carrier counts in `MHC_genetic_cluster_alleles.csv` for those comparisons. The FASTA export does not itself calculate cluster-level nucleotide diversity.

## How to run

From the project folder:

```bash
Rscript MHC_Allele_Analysis.R
```

**First run:** If you have not aligned the sequences yet, the script will create the filtered Excel and both FASTA files, then stop when it tries to read the missing `r_input/functional_alleles_alignment.fasta`. This is expected with the current script. Complete Step 3 in Geneious, save the aligned FASTA in `r_input/`, and run the script again to generate the final CSV and genetic-cluster FASTA files.

**Later runs:** With the aligned FASTA already present, the full script can run in one pass. If allele sequences or allele identifiers change after filtering, repeat the alignment step so the FASTA matches the new alleles.

## Output files

All outputs are in `r_output/`:

| File | Contents |
| --- | --- |
| `results_filtered.xlsx` | Filtered AmpliSAS results, with updated allele names and read summaries |
| `MHC_all_alleles.fasta` | All allele sequences, including labeled non-functional alleles |
| `MHC_functional_alleles.fasta` | Functional allele sequences for Geneious alignment |
| `MHC_sample_alleles.csv` | Individual metadata, disease information, allele counts, genotypes, nucleotide diversity, allele presence/absence, and supertype presence/counts |
| `MHC_genetic_cluster_alleles.csv` | Sample size and allele-carrier counts by genetic cluster |
| `MHC_location_allele_summary.csv` | Distinct total and functional alleles, plus average counts per bird, by location |
| `MHC_cluster_mean_alleles.csv` | Average total and functional alleles per bird by genetic cluster |
| `cluster_alignments/Cluster_<Genetic_cluster>.fasta` | One aligned FASTA per cluster containing each distinct functional allele present, once |

The original aligned FASTA is kept in `r_input/`, because it is a manual input to both the nucleotide diversity calculation and the cluster FASTA exports.

**Note:** The cluster FASTA outputs require the additional cluster-export code to be appended after the genetic-cluster CSV export in the R script.
