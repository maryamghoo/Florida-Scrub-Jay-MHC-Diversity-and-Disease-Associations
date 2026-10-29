# PLASMODIUM GWAS — MAIN MANHATTAN PLOT
# Chromosome 11 MHC-related BLAST location
# Four MHC-related unplaced scaffolds


# 1. PACKAGES

packages <- c("data.table", "ggplot2", "ggrepel")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(data.table)
library(ggplot2)
library(ggrepel)


# 2. FILE PATHS

data_dir <- "r_input"

out_dir <- "r_output"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list(
  incidence_raw = "incidence_gwas.assoc.logistic",
  incidence_adj = "incidence_gwas.assoc.logistic.adjusted",
  intensity_raw = "intensity_gwas.assoc.linear",
  intensity_adj = "intensity_gwas.assoc.linear.adjusted",
  chr_map = "chr_map.txt",
  regions = "MHC_BLAST_intervals_for_GWAS.tsv"
)

paths <- lapply(files, function(x) file.path(data_dir, x))

stopifnot(all(file.exists(unlist(paths))))


# 3. READ CHROMOSOME MAP

# chr_map.txt has an empty middle tab-separated column.

chr_map_raw <- fread(
  paths$chr_map,
  header = FALSE,
  sep = "\t",
  colClasses = "character",
  fill = TRUE
)

chr_map <- chr_map_raw[, .(
  CHR = trimws(V1),
  ChrName = trimws(V3)
)]

chr_map <- chr_map[
  !is.na(CHR) & !is.na(ChrName) &
    CHR != "" & ChrName != ""
]

chr_map <- unique(chr_map, by = "CHR")


# 4. READ BLAST INTERVALS

regions <- fread(paths$regions)

stopifnot(
  all(c("CHR", "ChrName", "start", "end") %in% names(regions))
)

regions[, CHR := as.character(CHR)]
regions[, start := as.numeric(start)]
regions[, end := as.numeric(end)]

regions <- regions[
  is.finite(start) & is.finite(end) &
    start > 0 & end >= start
]


# 5. SELECT CHROMOSOMES AND FOUR SCAFFOLDS

main_chr <- chr_map[
  grepl("^Chr([0-9]+[A-Z]?|Z)$", ChrName)
]

target_scaffolds <- c(
  "HiC_scaffold_56",
  "HiC_scaffold_75",
  "HiC_scaffold_97",
  "HiC_scaffold_182"
)

mhc_scaffolds <- chr_map[
  ChrName %in% target_scaffolds
]

selected_chr <- unique(
  rbindlist(list(main_chr, mhc_scaffolds)),
  by = "CHR"
)

stopifnot(
  all(c("Chr11", target_scaffolds) %in% selected_chr$ChrName)
)


# 6. READ GWAS RESULTS

read_gwas <- function(raw_file, adj_file, phenotype) {
  
  raw <- fread(raw_file)
  adj <- fread(adj_file)
  
  stopifnot(
    all(c("CHR", "SNP", "BP", "TEST", "P") %in% names(raw)),
    all(c("CHR", "SNP", "FDR_BH") %in% names(adj))
  )
  
  raw <- raw[
    TEST == "ADD",
    .(
      CHR = as.character(CHR),
      SNP = as.character(SNP),
      BP = as.numeric(BP),
      P = as.numeric(P)
    )
  ]
  
  adj <- adj[, .(
    CHR = as.character(CHR),
    SNP = as.character(SNP),
    FDR_BH = as.numeric(FDR_BH)
  )]
  
  if (
    anyDuplicated(raw[, .(CHR, SNP)]) ||
    anyDuplicated(adj[, .(CHR, SNP)])
  ) {
    stop("Duplicate chromosome/SNP combinations found.")
  }
  
  d <- merge(
    raw,
    adj,
    by = c("CHR", "SNP"),
    all.x = TRUE,
    sort = FALSE
  )
  
  d <- d[
    is.finite(P) &
      P >= 0 & P <= 1 &
      is.finite(BP) & BP > 0 &
      is.finite(FDR_BH)
  ]
  
  d[, logP := -log10(pmax(P, .Machine$double.xmin))]
  
  d[, Significant := FDR_BH < 0.05]
  
  d[, Phenotype := phenotype]
  
  d
}

incidence <- read_gwas(
  paths$incidence_raw,
  paths$incidence_adj,
  "A. Infection incidence"
)

intensity <- read_gwas(
  paths$intensity_raw,
  paths$intensity_adj,
  "B. Infection intensity"
)

gwas <- rbindlist(list(incidence, intensity))


# 7. SELECT GWAS DATA FOR PLOTTING

plot_gwas <- merge(
  gwas[CHR %in% selected_chr$CHR],
  selected_chr,
  by = "CHR",
  sort = FALSE
)

if (nrow(plot_gwas) == 0) {
  stop("No GWAS data found for selected chromosomes.")
}


# 8. CALCULATE CHROMOSOME DISPLAY WIDTHS

plot_chr <- plot_gwas[, .(
  chr_length = max(BP)
), by = .(CHR, ChrName)]

plot_chr[, chr_index := match(CHR, selected_chr$CHR)]

setorder(plot_chr, chr_index)

min_width <- max(plot_chr$chr_length) * 0.12

plot_chr[, display_width := pmax(chr_length, min_width)]

plot_chr[, display_offset :=
           shift(cumsum(display_width), fill = 0)
]

plot_chr[, display_center :=
           display_offset + display_width / 2
]


# 9. CALCULATE SNP DISPLAY POSITIONS

plot_gwas <- merge(
  plot_gwas,
  plot_chr[, .(
    CHR,
    chr_index,
    chr_length,
    display_width,
    display_offset
  )],
  by = "CHR",
  sort = FALSE
)

plot_gwas[, plot_x :=
            display_offset +
            (BP / chr_length) * display_width
]


# 10. CHROMOSOME COLORS

chromosome_colors <- c(
  "#355C7D",
  "#19A7A0",
  "#7B61A8",
  "#E6A04B"
)

plot_gwas[, PlotColor :=
            chromosome_colors[
              ((chr_index - 1) %% length(chromosome_colors)) + 1
            ]
]


# 11. FDR-SIGNIFICANT SNPS AND LABELS

plot_sig <- copy(plot_gwas[Significant == TRUE])

plot_sig[, SNP_label := paste0(
  ChrName,
  ": ",
  sprintf("%.2f", BP / 1e6),
  " Mb"
)]


# 12. MHC SCAFFOLD SHADING

mhc_display <- plot_chr[
  ChrName %in% target_scaffolds
]

mhc_display[, `:=`(
  xmin = display_offset,
  xmax = display_offset + display_width
)]


# 13. CHROMOSOME 11 BLAST LOCATION

# Two BLAST regions on chromosome 11.
# Use the midpoint of their combined genomic span
# as one approximate marker in the main figure.

chr11_intervals <- regions[ChrName == "Chr11"]

if (nrow(chr11_intervals) == 0) {
  stop("No Chr11 BLAST intervals found.")
}

chr11_start <- min(chr11_intervals$start)
chr11_end <- max(chr11_intervals$end)

chr11_mhc_midpoint <- (chr11_start + chr11_end) / 2

chr11_display <- copy(
  plot_chr[ChrName == "Chr11"]
)

stopifnot(nrow(chr11_display) == 1)

chr11_display[, mhc_plot_x :=
                display_offset +
                (chr11_mhc_midpoint / chr_length) * display_width
]

cat(
  "Chromosome 11 BLAST span:",
  chr11_start, "-", chr11_end, "\n"
)


# 14. PHENOTYPE ORDER

phenotype_order <- c(
  "A. Infection incidence",
  "B. Infection intensity"
)

plot_gwas[, Phenotype := factor(
  Phenotype,
  levels = phenotype_order
)]

plot_sig[, Phenotype := factor(
  Phenotype,
  levels = phenotype_order
)]


# 15. CREATE MAIN MANHATTAN PLOT

p_final <- ggplot() +
  
  # Subtle blue shading on four MHC-related scaffolds.
  # This highlights scaffold identity, not gene boundaries.
  
  geom_rect(
    data = mhc_display,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "#8AC9EC",
    alpha = 0.20,
    inherit.aes = FALSE
  ) +
  
  # Original GWAS SNPs
  
  geom_point(
    data = plot_gwas,
    aes(
      x = plot_x,
      y = logP,
      color = PlotColor
    ),
    size = 0.3,
    alpha = 0.7
  ) +
  
  scale_color_identity() +
  
  # One dashed line on Chr11 marking approximate
  # midpoint of the two BLAST hit regions.
  
  geom_vline(
    data = chr11_display,
    aes(xintercept = mhc_plot_x),
    color = "#0077B6",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # FDR-significant SNPs
  
  geom_point(
    data = plot_sig,
    aes(x = plot_x, y = logP),
    color = "#D62828",
    shape = 1,
    size = 2.3,
    stroke = 0.8
  ) +
  
  # Significant SNP labels
  
  geom_text_repel(
    data = plot_sig,
    aes(
      x = plot_x,
      y = logP,
      label = SNP_label
    ),
    color = "#B91C1C",
    size = 2.8,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    seed = 123
  ) +
  
  facet_grid(
    Phenotype ~ .,
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    breaks = plot_chr$display_center,
    labels = plot_chr$ChrName,
    expand = expansion(mult = c(0.01, 0.01)),
    guide = guide_axis(
      angle = 90,
      check.overlap = FALSE
    )
  ) +
  
  labs(
    title = "Genome-wide associations with Plasmodium infection",
    
    subtitle = paste(
      "Blue shading = MHC-related BLAST scaffolds |",
      "Dashed line = approximate Chr11 BLAST location |",
      "Red circles = FDR-significant SNPs (BH < 0.05)"
    ),
    
    x = "Chromosome / MHC-related scaffolds",
    
    y = expression(-log[10](italic(P)))
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 9),
    
    axis.text.x = element_text(
      size = 9,
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      color = "black"
    ),
    
    strip.background = element_rect(
      fill = "#E8EEF5",
      color = NA
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    panel.grid.major.y = element_line(
      color = "gray92",
      linewidth = 0.2
    ),
    
    legend.position = "none"
  )


# 16. DISPLAY AND SAVE

print(p_final)

ggsave(
  filename = file.path(
    out_dir,
    "Plasmodium_GWAS_MAIN_MHC.png"
  ),
  plot = p_final,
  width = 20,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)


# 17. PRINT SUMMARY

cat("\nMain Manhattan plot saved successfully.\n")
cat("Output folder:", out_dir, "\n")
cat("Chr11 BLAST midpoint:", chr11_mhc_midpoint, "bp\n")
cat("MHC-related scaffolds:\n")
print(mhc_display[, .(ChrName, CHR)])
