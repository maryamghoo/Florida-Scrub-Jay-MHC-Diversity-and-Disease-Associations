# Library
library(dada2)  # Trims primers and filters paired-end FASTQ reads.

# Input and output folders
working_dir <- getwd()
input_dir <- file.path(working_dir, "input")
output_dir <- file.path(working_dir, "output", "trimmed")

# Number of bases to remove from the start of R1 and R2
trim_bases <- 20

if (!dir.exists(input_dir)) stop("Input folder not found: ", input_dir)

# Find forward and reverse FASTQ files
fnFs <- sort(list.files(input_dir, pattern = "_R1_001\\.fastq(\\.gz)?$", full.names = TRUE))
fnRs <- sort(list.files(input_dir, pattern = "_R2_001\\.fastq(\\.gz)?$", full.names = TRUE))
if (!length(fnFs) || !length(fnRs)) stop("No paired FASTQ files found.")

# Match R1 and R2 files by sample name
sampleFs <- sub("_R1_001\\.fastq(\\.gz)?$", "", basename(fnFs))
sampleRs <- sub("_R2_001\\.fastq(\\.gz)?$", "", basename(fnRs))
if (anyDuplicated(sampleFs) || anyDuplicated(sampleRs) || !setequal(sampleFs, sampleRs)) {
  stop("R1 and R2 files must match one-to-one.")
}
fnRs <- fnRs[match(sampleFs, sampleRs)]

# Set up output file names
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
trimFs <- file.path(output_dir, paste0(sampleFs, "_R1_001_trim.fastq"))
trimRs <- file.path(output_dir, paste0(sampleFs, "_R2_001_trim.fastq"))

# Trim the first bases of each read and apply DADA2's default filters
out <- filterAndTrim(fnFs, trimFs, fnRs, trimRs,
                     trimLeft = c(trim_bases, trim_bases), compress = FALSE)

# Save a summary of reads before and after filtering
write.csv(data.frame(sample = sampleFs, out, check.names = FALSE),
          file.path(output_dir, "trimming_summary.csv"), row.names = FALSE)
