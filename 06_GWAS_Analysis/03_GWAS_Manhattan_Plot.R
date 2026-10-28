# PLASMODIUM GWAS MANHATTAN PLOT
# Named chromosomes + MHC-containing scaffold
# Unequal-width bins with minimum display width

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
# All input files are located in the r_input folder
data_dir <- "r_input"

files <- list(
  incidence_raw = "incidence_gwas.assoc.logistic",
  incidence_adj = "incidence_gwas.assoc.logistic.adjusted",
  intensity_raw = "intensity_gwas.assoc.linear",
  intensity_adj = "intensity_gwas.assoc.linear.adjusted",
  chr_map = "chr_map.txt"
)

paths <- lapply(files, function(f) file.path(data_dir, f))
out_dir <- file.path(data_dir, "GWAS_Manhattan_Figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# 3. READ CHROMOSOME MAP
chr_map <- as.data.table(
  read.table(
    paths$chr_map,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
)

setnames(chr_map, c("CHR", "ChrName"))

chr_map[, CHR := trimws(as.character(CHR))]
chr_map[, ChrName := trimws(as.character(ChrName))]

chr_map <- chr_map[
  !is.na(CHR) & !is.na(ChrName) &
    CHR != "" & ChrName != ""
]

chr_map <- unique(chr_map, by = "CHR")


# 4. DEFINE CHROMOSOMES TO DISPLAY
# MHC-containing scaffold from chromosome map
mhc_chr <- "JBGGOT010000146.1"
mhc_name <- "HiC_scaffold_182"

stopifnot(
  nrow(chr_map[
    CHR == mhc_chr & ChrName == mhc_name
  ]) == 1
)

# Keep named chromosomes only
main_chr <- chr_map[
  grepl("^Chr([0-9]+[A-Z]?|Z)$", ChrName)
]

# Add MHC scaffold at the end for visualization
selected_chr <- rbindlist(list(
  main_chr[, .(CHR, ChrName)],
  data.table(CHR = mhc_chr, ChrName = mhc_name)
))

selected_chr <- unique(selected_chr, by = "CHR")


# 5. READ GWAS RESULTS
read_gwas <- function(raw_file, adj_file, phenotype) {
  
  raw <- fread(raw_file)
  adj <- fread(adj_file)
  
  required_raw <- c("CHR", "SNP", "BP", "TEST", "P")
  required_adj <- c("CHR", "SNP", "FDR_BH")
  
  stopifnot(
    all(required_raw %in% names(raw)),
    all(required_adj %in% names(adj))
  )
  
  # Additive SNP association results only
  raw <- raw[TEST == "ADD", .(CHR, SNP, BP, P)]
  adj <- adj[, .(CHR, SNP, FDR_BH)]
  
  raw[, CHR := as.character(CHR)]
  adj[, CHR := as.character(CHR)]
  
  if (
    anyDuplicated(raw[, .(CHR, SNP)]) ||
    anyDuplicated(adj[, .(CHR, SNP)])
  ) {
    stop("Duplicate chromosome/SNP combinations found.")
  }
  
  # Match raw and adjusted p-values
  d <- merge(
    raw, adj,
    by = c("CHR", "SNP"),
    all.x = TRUE,
    sort = FALSE
  )
  
  # Keep valid association results
  d <- d[
    is.finite(P) &
      P >= 0 & P <= 1 &
      is.finite(BP) &
      BP > 0 &
      is.finite(FDR_BH)
  ]
  
  # Manhattan Y-axis uses unadjusted p-values
  d[, logP := -log10(
    pmax(P, .Machine$double.xmin)
  )]
  
  # Significance uses original full-GWAS FDR
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

# Combine both analyses before visualization filtering
gwas <- rbindlist(list(incidence, intensity))

rm(incidence, intensity)
gc()


# 6. FILTER FOR VISUALIZATION
# Exclude other unplaced scaffolds from the figure only.
# Original FDR-adjusted p-values remain unchanged.

plot_gwas <- gwas[CHR %in% selected_chr$CHR]

plot_gwas <- merge(
  plot_gwas,
  selected_chr,
  by = "CHR",
  all.x = TRUE,
  sort = FALSE
)


# 7. CALCULATE CHROMOSOME DISPLAY WIDTHS
# Maximum observed SNP position per chromosome/scaffold
plot_chr <- plot_gwas[, .(
  chr_length = max(BP)
), by = .(CHR, ChrName)]

# Preserve selected chromosome order
plot_chr[, chr_index := match(CHR, selected_chr$CHR)]
setorder(plot_chr, chr_index)

# Minimum width = 12% of largest observed chromosome span
min_width <- max(plot_chr$chr_length) * 0.12

# Larger chromosomes retain proportional widths.
# Smaller chromosomes receive the minimum width.
plot_chr[, display_width := pmax(chr_length, min_width)]

# Calculate cumulative display positions
plot_chr[, display_offset :=
           shift(cumsum(display_width), fill = 0)
]

plot_chr[, display_center :=
           display_offset + display_width / 2
]


# 8. CALCULATE SNP POSITIONS
plot_gwas <- merge(
  plot_gwas,
  plot_chr[, .(
    CHR, chr_index, chr_length,
    display_width, display_offset
  )],
  by = "CHR",
  all.x = TRUE,
  sort = FALSE
)

# SNP positions retain their relative location within
# each chromosome/scaffold.
plot_gwas[, plot_x :=
            display_offset + (BP / chr_length) * display_width
]


# 9. CHROMOSOME COLORS
chromosome_colors <- c(
  "#355C7D",  # Navy
  "#19A7A0",  # Teal
  "#7B61A8",  # Purple
  "#E6A04B"   # Orange
)

plot_gwas[, PlotColor := chromosome_colors[
  ((chr_index - 1) %% length(chromosome_colors)) + 1
]]


# 10. SIGNIFICANT SNPS
plot_sig <- plot_gwas[Significant == TRUE]

# Label each significant SNP with chromosome and Mb position
plot_sig[, SNP_label := paste0(
  ChrName, ": ",
  sprintf("%.2f", BP / 1e6),
  " Mb"
)]

# Display significant SNPs
print(
  plot_sig[, .(
    Phenotype, ChrName, SNP, BP, P, FDR_BH
  )]
)


# 11. MHC-CONTAINING SCAFFOLD
plot_mhc <- plot_chr[CHR == mhc_chr]

if (nrow(plot_mhc) != 1) {
  stop("MHC scaffold has no SNPs in the plotted data.")
}


# 12. FORMAT PHENOTYPE PANELS
phenotype_order <- c(
  "A. Infection incidence",
  "B. Infection intensity"
)

plot_gwas[, Phenotype := factor(
  Phenotype, levels = phenotype_order
)]

plot_sig[, Phenotype := factor(
  Phenotype, levels = phenotype_order
)]


# 13. CREATE MANHATTAN PLOT
p_final <- ggplot() +
  
  # Blue shading: MHC-containing scaffold
  geom_rect(
    data = plot_mhc,
    aes(
      xmin = display_offset,
      xmax = display_offset + display_width,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "#8AC9EC",
    alpha = 0.4,
    inherit.aes = FALSE
  ) +
  
  # Background SNPs
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
  
  # Significant SNPs in red
  geom_point(
    data = plot_sig,
    aes(x = plot_x, y = logP),
    color = "#D62828",
    size = 2
  ) +
  
  # Significant SNP labels: chromosome and Mb
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
  
  # Dashed line marking MHC-containing scaffold
  geom_vline(
    data = plot_mhc,
    aes(xintercept = display_center),
    color = "#0077B6",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  
  # Incidence and intensity panels
  facet_grid(
    Phenotype ~ .,
    scales = "free_y"
  ) +
  
  # Chromosome/scaffold names from chr_map.txt
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
      "Unequal-width chromosome bins (minimum display width) |",
      "Red = FDR-significant SNPs (BH < 0.05) |",
      "Blue = MHC-containing scaffold"
    ),
    x = "Chromosome / MHC-containing scaffold",
    y = expression(-log[10](italic(P)))
  ) +
  
  theme_classic(base_size = 11) +
  
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 9
    ),
    axis.text.x = element_text(
      size = 9,
      color = "black",
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 15)
    ),
    axis.title.y = element_text(size = 11),
    strip.background = element_rect(
      fill = "#E8EEF5",
      color = NA
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    panel.spacing = grid::unit(1, "lines"),
    panel.grid.major.y = element_line(
      color = "gray92",
      linewidth = 0.2
    ),
    legend.position = "none"
  )


# 14. DISPLAY AND SAVE PNG
ggsave(
  filename = file.path(
    out_dir, "Plasmodium_GWAS_FINAL.png"
  ),
  plot = p_final,
  width = 20,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
