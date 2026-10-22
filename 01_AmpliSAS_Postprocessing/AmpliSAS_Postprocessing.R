# Libraries
library(readxl)    # Reads the AmpliSAS Excel file.
library(openxlsx)  # Saves the filtered results as an Excel file.
library(ape)    # Handles DNA sequences in DNAbin format.
library(pegas)  # Calculates nucleotide diversity.

# Input and output directories
working_dir <- getwd()
input_dir <- file.path(working_dir, "r_input")
output_dir <- file.path(working_dir, "r_output")
input_file <- file.path(input_dir, "results.xlsx")
output_file <- file.path(output_dir, "results_filtered.xlsx")

#####---------- AMPLISAS results cleaning ----------#####

# Filtering settings
min_length <- 200
max_length <- 269
min_samples <- 2
min_depth_amplicon <- 700
duplicate_file <- file.path(input_dir, "duplicate_samples.txt")

# Read the complete sheet without treating any row as a header
if (!file.exists(input_file)) stop("Input file not found: ", input_file)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
sheet_name <- excel_sheets(input_file)[1]
data <- as.data.frame(read_excel(input_file, sheet = sheet_name, col_names = FALSE,
                                 col_types = "text", .name_repair = "minimal"),
                      stringsAsFactors = FALSE)

if (nrow(data) < 5 || ncol(data) < 10) stop("Unexpected AmpliSAS worksheet layout.")
if (!identical(as.character(unlist(data[4, c(1, 3, 5)])),
               c("SEQUENCE", "LENGTH", "SAMPLES")) ||
    data[1, 9] != "DEPTH_AMPLICON") {
  stop("Expected LENGTH/SAMPLES in row 4 and DEPTH_AMPLICON in row 1.")
}

# Keep allele rows meeting both length and sample-count requirements
alleles <- data[5:nrow(data), , drop = FALSE]
lengths <- suppressWarnings(as.numeric(alleles[[3]]))
sample_counts <- suppressWarnings(as.numeric(alleles[[5]]))
keep_rows <- !is.na(lengths) & !is.na(sample_counts) &
  lengths >= min_length & lengths <= max_length &
  sample_counts >= min_samples

# Keep allele rows that pass filtering
alleles <- alleles[keep_rows, , drop = FALSE]

# Keep sample columns with sufficient amplicon depth
amplicon_depths <- suppressWarnings(as.numeric(unlist(data[1, 10:ncol(data)])))
keep_cols <- c(seq_len(9), which(!is.na(amplicon_depths) &
                                   amplicon_depths >= min_depth_amplicon) + 9)

# Retain the top four header rows and all nine allele-information columns
filtered <- rbind(data[1:4, keep_cols, drop = FALSE],
                  alleles[, keep_cols, drop = FALSE])
rownames(filtered) <- NULL

# Remove duplicate sample columns listed in a text file (if provided)
if (file.exists(duplicate_file)) {
  duplicate_samples <- trimws(readLines(duplicate_file, warn = FALSE))
  duplicate_samples <- unique(duplicate_samples[nzchar(duplicate_samples)])
  sample_names <- trimws(as.character(unlist(filtered[4, ])))
  keep_cols <- seq_len(ncol(filtered)) <= 9 | !(sample_names %in% duplicate_samples)
  filtered <- filtered[, keep_cols, drop = FALSE]
}

# Recalculate DEPTH and SAMPLES from the remaining sample columns
if (nrow(filtered) > 4) {
  sample_reads <- filtered[5:nrow(filtered), -(1:9), drop = FALSE]
  read_counts <- matrix(
    suppressWarnings(as.numeric(unlist(sample_reads))),
    nrow = nrow(sample_reads), ncol = ncol(sample_reads)
  )
  read_counts[is.na(read_counts)] <- 0
  filtered[5:nrow(filtered), 4] <- as.character(rowSums(read_counts))
  filtered[5:nrow(filtered), 5] <- as.character(rowSums(read_counts > 0))
}

# Assign allele names
filtered[5:nrow(filtered), 9] <- sprintf(
  "ApcoIIB%03d", seq_len(nrow(filtered) - 4)
)

# Write the results without adding R-generated column names
wb <- createWorkbook()
addWorksheet(wb, sheet_name)
writeData(wb, sheet_name, filtered, colNames = FALSE, rowNames = FALSE,
          keepNA = FALSE)
saveWorkbook(wb, output_file, overwrite = TRUE)

#####---------- Export allele sequences as FASTA files ----------#####

sequences <- as.character(filtered[5:nrow(filtered), 1])
allele_names <- as.character(filtered[5:nrow(filtered), 9])

# Identify non-functional alleles
non_functional <- c("ApcoIIB002", "ApcoIIB022")

# Add stop codon information to FASTA headers
fasta_names <- ifelse(
  allele_names %in% non_functional,
  paste0(allele_names, " with stop codon"),
  allele_names
)

# Function to write FASTA files
write_fasta <- function(names, sequences, file) {
  fasta_lines <- as.vector(rbind(paste0(">", names), sequences))
  writeLines(fasta_lines, file)
}

# FASTA including all alleles
write_fasta(
  fasta_names,
  sequences,
  file.path(output_dir, "MHC_all_alleles.fasta")
)

# FASTA excluding non-functional alleles
keep <- !(allele_names %in% non_functional)

write_fasta(
  fasta_names[keep],
  sequences[keep],
  file.path(output_dir, "MHC_functional_alleles.fasta")
)

#####---------- Add additional sampling and disease info to prepare dataset ----------#####

# Read sample information
sampling_info <- read.delim(
  file.path(input_dir, "sampling_info.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

sampling_info$SampleID <- trimws(sampling_info$SampleID)

# Read disease information
disease_info <- read.csv(
  file.path(input_dir, "disease_info.csv"),
  header = TRUE,
  stringsAsFactors = FALSE
)

disease_info$SampleID <- trimws(disease_info$SampleID)

# Calculate log-transformed intensity
disease_info$Intensity <- as.numeric(disease_info$Intensity)
disease_info$logIntensity <- ifelse(
  disease_info$Intensity > 0,
  log10(disease_info$Intensity),
  NA_real_
)

# Get sample names and remove everything after the first underscore
sample_names <- trimws(as.character(unlist(filtered[4, 10:ncol(filtered)])))
sample_names <- sub("_.*$", "", sample_names)
anyDuplicated(sample_names)

# Get allele names and read counts
allele_names <- as.character(filtered[5:nrow(filtered), 9])
sample_reads <- filtered[5:nrow(filtered), 10:ncol(filtered), drop = FALSE]

read_counts <- matrix(
  suppressWarnings(as.numeric(unlist(sample_reads))),
  nrow = nrow(sample_reads),
  ncol = ncol(sample_reads)
)

read_counts[is.na(read_counts)] <- 0

# Create presence/absence matrix
all_alleles <- sprintf("ApcoIIB%03d", 1:38)

pa_matrix <- matrix(
  0L,
  nrow = length(sample_names),
  ncol = length(all_alleles),
  dimnames = list(NULL, all_alleles)
)

for (i in seq_along(all_alleles)) {
  allele_index <- which(allele_names == all_alleles[i])
  
  if (length(allele_index) > 0) {
    pa_matrix[, i] <- as.integer(
      colSums(read_counts[allele_index, , drop = FALSE] > 0) > 0
    )
  }
}

# Count total and functional alleles
functional_alleles <- setdiff(all_alleles, non_functional)

allele_counts <- data.frame(
  SampleID = sample_names,
  Total_Alleles = rowSums(pa_matrix),
  Functional_Alleles = rowSums(
    pa_matrix[, functional_alleles, drop = FALSE]
  ),
  check.names = FALSE
)

# Combine counts and presence/absence
result <- cbind(
  allele_counts,
  as.data.frame(pa_matrix, check.names = FALSE)
)

# Add Region, Location, and Genetic_cluster from sampling_info
result <- merge(
  sampling_info[, c(
    "SampleID", "Region", "Location", "Genetic_cluster",
    "Y", "X", "Year", "Month", "Sex", "Mass", "Age", "Age_rank"
  )],
  result,
  by = "SampleID",
  all.y = TRUE,
  sort = FALSE
)

# Add disease information
result <- merge(
  result,
  disease_info[, c("SampleID", "Incidence", "Intensity", "logIntensity")],
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

# Sort data by SampleID
result <- result[order(result$SampleID), ]

# Create genotype using functional alleles only
result$Genotype <- apply(
  result[, functional_alleles, drop = FALSE],
  1,
  function(x) {
    present <- functional_alleles[as.numeric(x) == 1]
    
    if (length(present) == 0) return(NA_character_)
    
    paste(present, collapse = "_")
  }
)

# Read the original aligned allele sequences
aligned_dna <- ape::read.dna(
  file.path(input_dir, "functional_alleles_alignment.fasta"),
  format = "fasta"
)

# Calculate nucleotide diversity using the original alignment
calculate_nuc_div <- function(genotype) {
  
  if (is.na(genotype) || !nzchar(genotype)) return(NA_real_)
  
  alleles <- strsplit(genotype, "_", fixed = TRUE)[[1]]
  
  # No alleles
  if (length(alleles) == 0) return(NA_real_)
  
  # One allele
  if (length(alleles) == 1) return(0)
  
  # Check allele names
  if (!all(alleles %in% rownames(aligned_dna))) return(NA_real_)
  
  # Extract aligned sequences
  dna <- aligned_dna[alleles, , drop = FALSE]
  
  # Calculate nucleotide diversity
  pegas::nuc.div(dna)
}

# Calculate nucleotide diversity for each individual
result$Nucleotide_Diversity <- vapply(
  result$Genotype,
  calculate_nuc_div,
  numeric(1)
)

# Read supertype assignments
supertypes <- read.csv(
  file.path(input_dir, "supertypes.csv"),
  stringsAsFactors = FALSE
)

# Create supertype presence/absence columns
for (st in c("A", "B", "C")) {
  alleles <- supertypes$allele[
    !is.na(supertypes$supertype) & supertypes$supertype == st
  ]
  
  alleles <- intersect(alleles, functional_alleles)
  
  result[[paste0("SupertypeII", st)]] <- as.integer(
    rowSums(result[, alleles, drop = FALSE] == 1, na.rm = TRUE) > 0
  )
}

# Count distinct supertypes per individual
result$Supertype_count <- rowSums(
  result[, c("SupertypeIIA", "SupertypeIIB", "SupertypeIIC")]
)

# Save final CSV
write.csv(
  result,
  file.path(output_dir, "MHC_sample_alleles.csv"),
  row.names = FALSE,
  na = ""
)

#####---------- Summarize MHC allele diversity by location and genetic cluster ----------#####

# Summarize distinct alleles and average allele counts by location
location_data <- result[
  !is.na(result$Location) & nzchar(trimws(result$Location)),
  , drop = FALSE
]

location_summary <- do.call(rbind, lapply(split(location_data, location_data$Location), function(x) {
  data.frame(
    Location = x$Location[1],
    SampleSize = nrow(x),
    Total_Different_Alleles = sum(colSums(x[, all_alleles, drop = FALSE] == 1) > 0),
    Total_Different_Functional_Alleles = sum(colSums(x[, functional_alleles, drop = FALSE] == 1) > 0),
    Average_Alleles_Per_Bird = mean(x$Total_Alleles),
    Average_Functional_Alleles_Per_Bird = mean(x$Functional_Alleles)
  )
}))
rownames(location_summary) <- NULL

write.csv(
  location_summary,
  file.path(output_dir, "MHC_location_allele_summary.csv"),
  row.names = FALSE
)

# Summarize average allele counts by genetic cluster
cluster_mean_data <- result[
  !is.na(result$Genetic_cluster) & nzchar(trimws(result$Genetic_cluster)),
  , drop = FALSE
]

cluster_mean_summary <- aggregate(
  cbind(Total_Alleles, Functional_Alleles) ~ Genetic_cluster,
  data = cluster_mean_data,
  FUN = mean
)
names(cluster_mean_summary)[2:3] <- c(
  "Average_Alleles_Per_Bird",
  "Average_Functional_Alleles_Per_Bird"
)

write.csv(
  cluster_mean_summary,
  file.path(output_dir, "MHC_cluster_mean_alleles.csv"),
  row.names = FALSE
)

#####---------- Extract allele frequency infor for genetic clusters ----------#####

# Exclude samples without genetic cluster information
cluster_data <- result[
  !is.na(result$Genetic_cluster) &
    nzchar(trimws(result$Genetic_cluster)),
  ,
  drop = FALSE
]

# Count individuals and allele presence within each genetic cluster
cluster_summary <- aggregate(
  cluster_data[, all_alleles],
  by = list(Genetic_cluster = cluster_data$Genetic_cluster),
  FUN = sum
)

# Calculate sample size for each genetic cluster
sample_sizes <- aggregate(
  SampleID ~ Genetic_cluster,
  data = cluster_data,
  FUN = function(x) length(unique(x))
)

names(sample_sizes)[2] <- "SampleSize"

# Combine sample sizes and allele counts
cluster_summary <- merge(
  sample_sizes,
  cluster_summary,
  by = "Genetic_cluster",
  sort = FALSE
)

# Arrange genetic cluster summary columns
cluster_summary <- cluster_summary[, c(
  "Genetic_cluster",
  "SampleSize",
  all_alleles
)]

# Save CSV
write.csv(
  cluster_summary,
  file.path(output_dir, "MHC_genetic_cluster_alleles.csv"),
  row.names = FALSE
)

#####---------- Export aligned functional allele sequences for each genetic cluster ----------#####

# Create output directory
cluster_fasta_dir <- file.path(output_dir, "cluster_alignments")
dir.create(cluster_fasta_dir, recursive = TRUE, showWarnings = FALSE)

# Get genetic clusters
clusters <- unique(na.omit(result$Genetic_cluster))
clusters <- clusters[nzchar(trimws(as.character(clusters)))]

for (cluster in clusters) {
  
  # Select individuals belonging to the cluster
  cluster_data <- result[
    result$Genetic_cluster == cluster &
      !is.na(result$Genetic_cluster),
    ,
    drop = FALSE
  ]
  
  # Identify functional alleles present in the cluster
  present_alleles <- functional_alleles[
    colSums(cluster_data[, functional_alleles, drop = FALSE] == 1) > 0
  ]
  
  # Keep alleles available in the alignment
  present_alleles <- present_alleles[
    present_alleles %in% rownames(aligned_dna)
  ]
  
  if (length(present_alleles) == 0) next
  
  # Extract aligned sequences
  cluster_alignment <- aligned_dna[
    present_alleles, ,
    drop = FALSE
  ]
  
  # Save FASTA file
  output_fasta <- file.path(
    cluster_fasta_dir,
    paste0("Cluster_", gsub("[^A-Za-z0-9_-]", "_", cluster), ".fasta")
  )
  
  ape::write.dna(
    cluster_alignment,
    file = output_fasta,
    format = "fasta",
    nbcol = -1,
    colw = 1000000
  )
}