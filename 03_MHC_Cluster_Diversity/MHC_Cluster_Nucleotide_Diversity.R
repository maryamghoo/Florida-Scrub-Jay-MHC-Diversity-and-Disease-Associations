# Libraries
library(ape)
library(pegas)
library(boot)
library(ggplot2)

# Input and output folders
working_dir <- getwd()
input_dir <- file.path(working_dir, "r_input")
output_dir <- file.path(working_dir, "r_output")
counts_file <- file.path(input_dir, "MHC_genetic_cluster_alleles.csv")
alignment_dir <- file.path(input_dir, "cluster_alignments")
genome_pi_dir <- file.path(input_dir, "genome_pi")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Analysis settings
non_functional <- c("ApcoIIB002", "ApcoIIB022")
bootstrap_replicates <- 10000
set.seed(123)

# Read cluster allele counts and exclude Admixed
counts <- read.csv(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
counts <- counts[!is.na(counts$Genetic_cluster) &
                   nzchar(trimws(counts$Genetic_cluster)) &
                   trimws(counts$Genetic_cluster) != "Admixed", , drop = FALSE]

cluster_order <- c("Merritt Island", "SL", "CBC", "SBC", "Archbold", "Ocala")
cluster_order <- c(intersect(cluster_order, counts$Genetic_cluster),
                   setdiff(counts$Genetic_cluster, cluster_order))
counts <- counts[match(cluster_order, counts$Genetic_cluster), , drop = FALSE]

# Calculate unweighted and carrier-weighted diversity for each cluster
cluster_results <- lapply(seq_len(nrow(counts)), function(i) {
  cluster <- counts$Genetic_cluster[i]
  filename <- paste0("Cluster_", gsub("[^A-Za-z0-9_-]", "_", cluster), ".fasta")
  dna <- ape::read.dna(file.path(alignment_dir, filename), format = "fasta")
  allele_ids <- intersect(rownames(dna), names(counts))
  allele_ids <- setdiff(allele_ids, non_functional)
  carrier_counts <- as.numeric(counts[i, allele_ids, drop = TRUE])
  carrier_counts[is.na(carrier_counts)] <- 0
  names(carrier_counts) <- allele_ids
  allele_ids <- allele_ids[carrier_counts > 0]
  carrier_counts <- carrier_counts[allele_ids]
  dna <- dna[allele_ids, , drop = FALSE]
  n_alleles <- length(allele_ids)
  total_carriers <- sum(carrier_counts)

  # Nucleotide diversity and segregating sites
  original <- if (n_alleles >= 2)
    pegas::nuc.div(dna, variance = TRUE, pairwise.deletion = FALSE)
  else c(NA_real_, NA_real_)
  segregating_sites <- if (n_alleles >= 2) length(ape::seg.sites(dna)) else 0L

  # Carrier-weighted pairwise diversity
  weighted_pi <- NA_real_
  if (n_alleles >= 2 && total_carriers >= 2) {
    distances <- as.matrix(ape::dist.dna(dna, model = "raw", pairwise.deletion = FALSE))
    weights <- outer(carrier_counts, carrier_counts)
    weighted_pi <- sum(distances[upper.tri(distances)] *
                       weights[upper.tri(weights)]) / choose(total_carriers, 2)
  }

  # Bootstrap confidence intervals for unweighted diversity
  ci_lower <- ci_upper <- NA_real_
  bootstrap_values <- numeric(0)
  if (n_alleles >= 3) {
    b <- boot::boot(
      data = seq_len(n_alleles),
      statistic = function(data, indices) {
        pegas::nuc.div(dna[data[indices], , drop = FALSE], pairwise.deletion = FALSE)
      },
      R = bootstrap_replicates
    )
    bootstrap_values <- as.numeric(b$t)
    ci <- tryCatch(boot::boot.ci(b, type = "bca"), error = function(e) NULL)
    if (!is.null(ci) && !is.null(ci$bca)) {
      ci_lower <- ci$bca[4]
      ci_upper <- ci$bca[5]
    }
  }

  list(
    stats = data.frame(
      Genetic_cluster = cluster,
      SampleSize = counts$SampleSize[i],
      Distinct_Functional_Alleles = n_alleles,
      Total_Allele_Carrier_Observations = total_carriers,
      Alignment_Length = ncol(dna),
      Segregating_Sites = segregating_sites,
      Proportion_Segregating = segregating_sites / ncol(dna),
      Pi_Unweighted = unname(original[1]),
      Unweighted_Variance = unname(original[2]),
      Unweighted_CI_Lower = ci_lower,
      Unweighted_CI_Upper = ci_upper,
      Pi_Carrier_Weighted = weighted_pi
    ),
    bootstrap = data.frame(Genetic_cluster = cluster, Pi = bootstrap_values)
  )
})

# Save cluster diversity summary
summary <- do.call(rbind, lapply(cluster_results, `[[`, "stats"))
write.csv(summary, file.path(output_dir, "MHC_cluster_diversity.csv"), row.names = FALSE)

# Plot unweighted and carrier-weighted diversity
plot_data <- rbind(
  data.frame(Genetic_cluster = summary$Genetic_cluster,
             Method = "Unweighted", Pi = summary$Pi_Unweighted),
  data.frame(Genetic_cluster = summary$Genetic_cluster,
             Method = "Carrier-weighted", Pi = summary$Pi_Carrier_Weighted)
)
plot_data$Genetic_cluster <- factor(plot_data$Genetic_cluster,
                                    levels = summary$Genetic_cluster)
plot_diversity <- ggplot(plot_data, aes(Genetic_cluster, Pi, shape = Method)) +
  geom_point(size = 3, position = position_dodge(width = 0.35)) +
  theme_minimal(base_size = 14) +
  labs(x = "Genetic Cluster", y = expression(italic(pi)),
       title = "MHC nucleotide diversity by genetic cluster") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "MHC_cluster_diversity.png"), plot_diversity,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(output_dir, "MHC_cluster_diversity.pdf"), plot_diversity,
       width = 10, height = 5)

# Plot unweighted diversity with confidence intervals
original_df <- data.frame(
  Cluster = factor(summary$Genetic_cluster, levels = summary$Genetic_cluster),
  Pi = summary$Pi_Unweighted,
  LowerCI = summary$Unweighted_CI_Lower,
  UpperCI = summary$Unweighted_CI_Upper
)

plot_original <- ggplot(original_df, aes(x = Cluster, y = Pi)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.1, linewidth = 0.8) +
  theme_minimal(base_size = 16) +
  labs(
    title = expression("Nucleotide Diversity (" * italic(pi) * ") and Confidence Intervals by Genetic Clusters"),
    x = "Genetic Cluster",
    y = expression(italic(pi))
  ) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

ggsave(file.path(output_dir, "pi_clusteres.jpg"), plot_original,
       width = 10, height = 5, units = "in", dpi = 300)

# Prepare MHC bootstrap data and cluster colors
boot_df <- do.call(rbind, lapply(cluster_results, `[[`, "bootstrap"))
boot_df <- boot_df[is.finite(boot_df$Pi), , drop = FALSE]
boot_df$Genetic_cluster <- factor(boot_df$Genetic_cluster,
                                 levels = summary$Genetic_cluster)
unit_colors <- c(
  "Merritt Island" = "#d62728", "SL" = "royalblue1", "CBC" = "navy",
  "SBC" = "chartreuse3", "Archbold" = "mediumorchid4", "Ocala" = "darkorange"
)
library(patchwork)

# Genome-wide summary values
genome_info <- data.frame(
  Genetic_cluster = c("Archbold", "CBC", "Merritt Island", "Ocala", "SBC", "SL"),
  Genome_pi = c(0.919e-3, 0.818e-3, 0.903e-3, 0.915e-3, 0.810e-3, 0.867e-3),
  Ho = c(0.329, 0.287, 0.272, 0.318, 0.289, 0.435),
  He = c(0.338, 0.298, 0.284, 0.321, 0.308, 0.477)
)

# Combine MHC and genome-wide summaries
comparison <- merge(summary, genome_info, by = "Genetic_cluster", sort = FALSE)
comparison <- comparison[match(summary$Genetic_cluster, comparison$Genetic_cluster), ]
comparison$MHC_avg_alleles <- comparison$Total_Allele_Carrier_Observations /
  comparison$SampleSize
comparison$Genetic_cluster <- factor(comparison$Genetic_cluster,
                                     levels = summary$Genetic_cluster)
write.csv(comparison, file.path(output_dir, "MHC_genome_comparison.csv"), row.names = FALSE)

# Calculate Spearman correlations
correlation_pairs <- list(
  c("Pi_Unweighted", "Genome_pi"),
  c("Pi_Unweighted", "Ho"),
  c("Distinct_Functional_Alleles", "Genome_pi"),
  c("MHC_avg_alleles", "Ho")
)
correlation_results <- do.call(rbind, lapply(correlation_pairs, function(pair) {
  x <- comparison[[pair[1]]]
  y <- comparison[[pair[2]]]
  valid <- is.finite(x) & is.finite(y)
  test <- suppressWarnings(cor.test(x[valid], y[valid], method = "spearman", exact = FALSE))
  data.frame(MHC_Measure = pair[1], Genome_Measure = pair[2],
             N_Clusters = sum(valid), Spearman_rho = unname(test$estimate),
             P_value = test$p.value)
}))
write.csv(correlation_results, file.path(output_dir, "MHC_genome_correlations.csv"),
          row.names = FALSE)

theme_pub <- theme_classic(base_size = 12) +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank(), legend.text = element_text(size = 14))

# Create MHC versus genome-wide scatterplots
comparison_scatter <- function(x, y, x_label, y_label) {
  ggplot(comparison, aes(x = .data[[x]], y = .data[[y]], fill = Genetic_cluster)) +
    geom_point(shape = 21, size = 3.5, colour = "black", stroke = 0.5) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE,
                color = "black", linewidth = 0.8) +
    scale_fill_manual(values = unit_colors) +
    labs(x = x_label, y = y_label) + theme_pub
}

p1 <- comparison_scatter("Genome_pi", "Pi_Unweighted",
                         expression("Genome-wide nucleotide diversity (" * pi * ")"),
                         expression("MHC nucleotide diversity (" * pi * ")"))
p2 <- comparison_scatter("Ho", "Pi_Unweighted",
                         expression("Observed heterozygosity (" * H[O] * ")"),
                         expression("MHC nucleotide diversity (" * pi * ")"))
p3 <- comparison_scatter("Genome_pi", "Distinct_Functional_Alleles",
                         expression("Genome-wide nucleotide diversity (" * pi * ")"),
                         "Number of functional MHC alleles")
p4 <- comparison_scatter("Ho", "MHC_avg_alleles",
                         expression("Observed heterozygosity (" * H[O] * ")"),
                         "Mean functional MHC alleles per individual")

supplemental_plot <- (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", legend.direction = "horizontal",
        plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(output_dir, "MHC_genome_comparison.png"),
       supplemental_plot, width = 10, height = 8, dpi = 600, bg = "white")
ggsave(file.path(output_dir, "MHC_genome_comparison.pdf"),
       supplemental_plot, width = 10, height = 8, bg = "white")

# Read genome-wide windowed nucleotide diversity
window_files <- c(
  "Merritt Island" = "merritt.windowed.pi", "SL" = "sl.windowed.pi",
  "CBC" = "cbc.windowed.pi", "SBC" = "sbc.windowed.pi",
  "Archbold" = "archbold.windowed.pi", "Ocala" = "ocala.windowed.pi"
)

missing_window_files <- window_files[!file.exists(file.path(genome_pi_dir, window_files))]
if (length(missing_window_files) == 0) {
  genome_pi_df <- do.call(rbind, lapply(names(window_files), function(cluster) {
    windows <- read.table(file.path(genome_pi_dir, window_files[[cluster]]), header = TRUE)
    windows <- windows[!is.na(windows$N_VARIANTS) & windows$N_VARIANTS >= 20 &
                         is.finite(windows$PI), , drop = FALSE]
    data.frame(Genetic_cluster = cluster, Pi = windows$PI)
  }))
  genome_pi_df$Genetic_cluster <- factor(genome_pi_df$Genetic_cluster,
                                        levels = summary$Genetic_cluster)
  write.csv(genome_pi_df, file.path(output_dir, "MHC_genome_windowed_pi.csv"),
            row.names = FALSE)

  genome_window_summary <- aggregate(Pi ~ Genetic_cluster, data = genome_pi_df,
                                     FUN = function(x) c(Mean = mean(x), SD = sd(x),
                                                        N_Windows = length(x)))
  genome_window_summary <- data.frame(
    Genetic_cluster = genome_window_summary$Genetic_cluster,
    genome_window_summary$Pi, row.names = NULL)
  write.csv(genome_window_summary,
            file.path(output_dir, "MHC_genome_windowed_pi_summary.csv"), row.names = FALSE)

  # Plot MHC bootstrap and genome-wide windowed diversity
  shared_theme <- theme_minimal(base_size = 16) +
    theme(legend.position = "none", panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 16, face = "bold"),
          axis.text = element_text(size = 14),
          axis.text.x = element_text(angle = 45, hjust = 1))

  mhc_panel <- ggplot(boot_df, aes(Genetic_cluster, Pi, fill = Genetic_cluster)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, color = "black", linewidth = 0.6) +
    geom_point(data = original_df, aes(x = Cluster, y = Pi), inherit.aes = FALSE,
               size = 2.8, shape = 21, fill = "white", color = "black", stroke = 0.8) +
    scale_fill_manual(values = unit_colors) +
    labs(title = expression("MHC Nucleotide Diversity (" * italic(pi) * ")"),
         x = "Genetic Cluster", y = expression("MHC " * italic(pi))) + shared_theme

  genome_panel <- ggplot(genome_pi_df, aes(Genetic_cluster, Pi, fill = Genetic_cluster)) +
    geom_boxplot(width = 0.65, outlier.size = 1.2,
                 color = "black", linewidth = 0.6) +
    scale_fill_manual(values = unit_colors) +
    labs(title = expression("Genome-wide Nucleotide Diversity (" * italic(pi) * ")"),
         x = "Genetic Cluster", y = expression("Genome-wide " * italic(pi))) + shared_theme

  # Save the vertical two-panel figure
  plot_vertical <- mhc_panel / genome_panel + patchwork::plot_layout(ncol = 1) +
    plot_annotation(tag_levels = "A")

  ggsave(file.path(output_dir, "MHC_genomewide_pi_vertical.png"),
         plot_vertical, width = 9, height = 11, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, "MHC_genomewide_pi_vertical.pdf"),
         plot_vertical, width = 9, height = 11, bg = "white")
}
