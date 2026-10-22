library(dplyr)  # For filtering, selecting, and organizing data
library(readr)  # For reading and writing CSV files

# Set input and output directories
working_dir <- getwd()
input_dir <- file.path(working_dir, "r_input")
output_dir <- file.path(working_dir, "r_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read whole-genome IDs and MHC data
pops <- read.table(
  file.path(input_dir, "pop140.txt"),
  header = FALSE,
  stringsAsFactors = FALSE
)

mhc <- read_csv(
  file.path(input_dir, "MHC_sample_alleles.csv"),
  show_col_types = FALSE
)

# Exclude Admixed individuals
mhc <- mhc %>%
  filter(!is.na(Genetic_cluster), Genetic_cluster != "Admixed")

# Identify shared SampleIDs
shared_ids <- intersect(pops$V1, mhc$SampleID)

# Save shared SampleIDs without a header
write.table(
  shared_ids,
  file.path(output_dir, "shared_sampleIDs.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# Retain shared individuals in the same order
mhc_shared <- mhc %>%
  filter(SampleID %in% shared_ids) %>%
  slice(match(shared_ids, SampleID))

# Prepare incidence dataset
mhc_incidence <- mhc_shared %>%
  filter(!is.na(Incidence)) %>%
  select(SampleID, Incidence)

write_csv(
  mhc_incidence,
  file.path(output_dir, "mhc_incidence.csv")
)

# Prepare intensity dataset for infected individuals
mhc_intensity <- mhc_shared %>%
  filter(Incidence == 1, !is.na(logIntensity)) %>%
  select(SampleID, logIntensity)

write_csv(
  mhc_intensity,
  file.path(output_dir, "mhc_intensity.csv")
)