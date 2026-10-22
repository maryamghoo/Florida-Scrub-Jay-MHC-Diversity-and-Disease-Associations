# FASTQ preprocessing: trimming, optional line-ending fix, and merging

These scripts prepare paired-end FASTQ reads and merge overlapping R1/R2 reads.

1. `trim_primers.R` removes bases from the start of each read using DADA2.
2. `fix_trimmed_fastq.sh` fixes Windows-style line endings, only when necessary.
3. `merge_reads.sh` merges paired-end reads using PEAR.

## Requirements

- **R** and **DADA2** for trimming. Install DADA2 in R:

  ```r
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("dada2")
  ```

- **Bash**, **SLURM** (`sbatch`), and standard Linux tools for the line-ending job.
- **PEAR** installed and available as `pear` on the cluster for the merging job. Load its module first if your cluster uses environment modules.

The R script uses `library(dada2)` to trim and filter reads. No additional R packages are required.

## Folder setup

```text
project/
├── trim_primers.R
├── fix_trimmed_fastq.sh
├── merge_reads.sh
├── README.md
├── input/
│   ├── Sample_1_R1_001.fastq
│   └── Sample_1_R2_001.fastq
└── output/
    ├── trimmed/
    │   ├── Sample_1_R1_001_trim.fastq
    │   ├── Sample_1_R2_001_trim.fastq
    │   └── trimming_summary.csv
    ├── fixed/                  # only if the fix job is run
    │   ├── Sample_1_R1_001_trim.fastq
    │   └── Sample_1_R2_001_trim.fastq
    └── merged/
        ├── Sample_1.assembled.fastq
        ├── Sample_1.unassembled.forward.fastq
        ├── Sample_1.unassembled.reverse.fastq
        └── Sample_1.discarded.fastq
```

Scripts create their own output directories. The `merged` filenames shown above are PEAR's normal outputs; some files, such as discarded reads, may be empty.

## Step 1: Trim reads

Put matching R1 and R2 FASTQ files in `input/`. Accepted input extensions are `.fastq` and `.fastq.gz`.

```bash
Rscript trim_primers.R
```

You can also run the script in RStudio after setting the working directory to `project/`.

Settings at the top of `trim_primers.R`:

```r
working_dir <- getwd()
input_dir <- file.path(working_dir, "input")
output_dir <- file.path(working_dir, "output", "trimmed")
trim_bases <- 20
```

The script removes the first **20 bases** from each R1 and R2 read by default and applies DADA2's default filters. It does not search for or confirm primer sequences. It writes uncompressed FASTQ files to `output/trimmed/` and `trimming_summary.csv` with read counts. R1 and R2 remain separate.

## Step 2: Fix FASTQ line endings (optional)

Only use this job if trimmed FASTQ files have Windows (`CRLF`) line endings. On Linux, check a file with:

```bash
file output/trimmed/Sample_1_R1_001_trim.fastq
```

If the output mentions `CRLF line terminators`, submit:

```bash
sbatch fix_trimmed_fastq.sh
```

The job reads `output/trimmed/` and saves corrected copies in `output/fixed/`, without changing the originals. It removes carriage returns at the end of lines; it does not change DNA sequences or merge reads. If the files already use Unix (`LF`) endings, **skip this job**.

**Important:** Fixed files keep exactly the same names as trimmed files, such as `Sample_1_R1_001_trim.fastq`. No manual renaming is needed before merging.

## Step 3: Merge paired-end reads

Open `merge_reads.sh` and set `source_dir` near the top:

```bash
source_dir="trimmed"  # Use "fixed" if you ran the line-ending job
```

- Use `trimmed` if you skipped Step 2.
- Use `fixed` if you completed Step 2.

Both versions have the same filename format, so PEAR can use either one directly. Then, from the project directory, submit:

```bash
sbatch merge_reads.sh
```

For each matched pair, PEAR writes files under `output/merged/`:

- `Sample_1.assembled.fastq`: reads PEAR successfully merged.
- `Sample_1.unassembled.forward.fastq`: forward reads that did not merge.
- `Sample_1.unassembled.reverse.fastq`: reverse reads that did not merge.
- `Sample_1.discarded.fastq`: reads PEAR discarded.

The `.assembled.fastq` files contain the merged sequences for downstream analysis. PEAR requires enough overlap between forward and reverse reads; not all pairs will necessarily merge. Check the job log (`merge_reads_JOBID.out` and `.err`) for the results. The original input files remain unchanged.

## Step 4: Identify MHC alleles using AmpliSAS

The merged FASTQ files in `output/merged/` (`*.assembled.fastq`) are used as inputs to the AmpliSAS web server to identify MHC alleles, as described in the Methods section of Ghoojaei et al. (2026), **“Immunogenetic diversity, selection, and genome-wide disease associations in the declining Florida scrub-jay (*Aphelocoma coerulescens*).”** This step is performed on the web server; no additional code is included here.
