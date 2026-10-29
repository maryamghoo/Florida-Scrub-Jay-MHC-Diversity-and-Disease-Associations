# CREATE MHC BLAST INTERVALS FOR GWAS
# Output: MHC_BLAST_intervals_for_GWAS.tsv

library(data.table)

# 1. FILE PATHS

data_dir <- "r_input"

blast_files <- c(
  "chromosome 11.txt",
  "HiC_scaffold_56.txt",
  "HiC_scaffold_75.txt",
  "HiC_scaffold_97.txt",
  "HiC_scaffold_182.txt"
)

# 2. READ CHROMOSOME MAP

chr_map <- as.data.table(
  read.table(
    file.path(data_dir, "chr_map.txt"),
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE
  )
)

setnames(chr_map, c("CHR", "ChrName"))

chr_map[, CHR := trimws(as.character(CHR))]
chr_map[, ChrName := trimws(as.character(ChrName))]

chr_map <- unique(chr_map, by = "CHR")

# Verify chromosome IDs
print(
  chr_map[
    ChrName %in% c(
      "Chr11",
      "HiC_scaffold_56",
      "HiC_scaffold_75",
      "HiC_scaffold_97",
      "HiC_scaffold_182"
    )
  ]
)

# 3. EXTRACT BLAST COORDINATES FROM FASTA HEADERS

extract_intervals <- function(file) {
  
  lines <- readLines(file, warn = FALSE)
  
  # FASTA headers start with >
  headers <- grep("^>", lines, value = TRUE)
  
  # Extract accession:start-end
  pattern <- "^>([^:[:space:]]+):([0-9]+)-([0-9]+)"
  
  valid <- grepl(pattern, headers)
  headers <- headers[valid]
  
  matches <- regmatches(
    headers,
    regexec(pattern, headers)
  )
  
  data.table(
    Accession = vapply(matches, `[`, character(1), 2),
    start = as.numeric(
      vapply(matches, `[`, character(1), 3)
    ),
    end = as.numeric(
      vapply(matches, `[`, character(1), 4)
    ),
    ChrName = if (grepl("chromosome 11", basename(file),
                        fixed = TRUE)) {
      "Chr11"
    } else {
      tools::file_path_sans_ext(basename(file))
    }
  )
}

# Read all five files

regions_raw <- rbindlist(
  lapply(
    file.path(data_dir, blast_files),
    extract_intervals
  )
)

# 4. REMOVE DUPLICATE INTERVALS

regions_raw <- unique(regions_raw)

# 5. MERGE OVERLAPPING INTERVALS

setorder(regions_raw, ChrName, start, end)

merge_intervals <- function(d) {
  
  d <- copy(d)
  setorder(d, start, end)
  
  # Running maximum handles nested overlapping intervals
  d[, running_end := cummax(end)]
  
  # Start a new group only when intervals do not overlap
  d[, group := cumsum(
    start > shift(running_end, fill = -Inf)
  )]
  
  d[, .(
    start = min(start),
    end = max(end)
  ), by = group][, group := NULL][]
}

regions_merged <- regions_raw[
  ,
  merge_intervals(.SD),
  by = ChrName
]

# 6. MATCH BLAST REGIONS TO GWAS CHROMOSOME IDs

regions_final <- merge(
  regions_merged,
  chr_map,
  by = "ChrName",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(regions_final$CHR)) {
  stop("Some BLAST regions could not be matched to chr_map.txt")
}

# 7. ADD EVIDENCE LABELS

# These labels describe the original BLAST hit groups.
# They do not independently verify percent identity.

regions_final[, evidence := fifelse(
  ChrName == "Chr11",
  "partial BLAST match",
  fifelse(
    ChrName == "HiC_scaffold_182",
    "100% identity hit scaffold",
    "other BLAST match"
  )
)]

# 8. ORDER COLUMNS

setcolorder(
  regions_final,
  c("CHR", "ChrName", "start", "end", "evidence")
)

setorder(regions_final, ChrName, start)

# 9. SAVE TSV

output_file <- file.path(
  data_dir,
  "MHC_BLAST_intervals_for_GWAS.tsv"
)

fwrite(
  regions_final,
  output_file,
  sep = "\t"
)

