# Libraries
library(ggplot2)  # Plots allele frequencies.
library(RColorBrewer) # Supplies allele plot colors.

# Input and output directories
working_dir <- getwd()
input_dir <- file.path(working_dir, "r_input")
output_dir <- file.path(working_dir, "r_output")
input_file <- file.path(input_dir, "MHC_genetic_cluster_alleles.csv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read genetic-cluster allele counts
cluster_data <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
# Exclude admixed group
cluster_data <- cluster_data[!is.na(cluster_data$Genetic_cluster) &
                               cluster_data$Genetic_cluster != "Admixed", , drop = FALSE]

all_alleles <- sprintf("ApcoIIB%03d", 1:38)

# Convert allele counts to numbers
cluster_data[all_alleles] <- lapply(cluster_data[all_alleles], as.numeric)
# Calculate allele frequencies relative to total allele-presence counts per cluster
# This is NOT the proportion of individuals carrying an allele (count / SampleSize).
allele_totals <- rowSums(cluster_data[all_alleles])
frequency <- cluster_data[c("Genetic_cluster", "SampleSize")]
frequency$Total_Allele_Observations <- allele_totals
for (a in all_alleles) {
  frequency[[a]] <- ifelse(allele_totals > 0, cluster_data[[a]] / allele_totals, NA_real_)
}
write.csv(frequency, file.path(output_dir, "MHC_cluster_allele_frequencies.csv"), row.names = FALSE)

# Long format for plotting
long_freq <- data.frame(
  Location = rep(frequency$Genetic_cluster, each = length(all_alleles)),
  variable = rep(all_alleles, times = nrow(frequency)),
  value = as.vector(t(as.matrix(frequency[all_alleles])))
)
long_freq$Location <- factor(long_freq$Location,
                             levels = c("Merritt Island", "SL", "CBC", "SBC", "Archbold", "Ocala"))
long_freq$variable <- factor(long_freq$variable, levels = all_alleles)

# Use the original custom allele labels
custom_labels <- function(x) {
  x <- as.character(x)
  x[x == "ApcoIIB016"] <- "bolditalic('★*016')"
  x[x == "ApcoIIB006"] <- "bolditalic('*006')"
  x[x == "ApcoIIB004"] <- "bolditalic('▲*004')"
  x[x == "ApcoIIB008"] <- "bolditalic('■*008')"
  x[x == "ApcoIIB002"] <- "bolditalic('⬤*002')"
  x[x == "ApcoIIB022"] <- "bolditalic('⬤*022')"
  x[x == "ApcoIIB001"] <- "bolditalic('*001')"
  x[x == "ApcoIIB003"] <- "bolditalic('★*003')"
  x[x == "ApcoIIB005"] <- "bolditalic('*005')"
  x[x == "ApcoIIB007"] <- "bolditalic('*007')"
  x[x == "ApcoIIB009"] <- "bolditalic('*009')"
  x[x == "ApcoIIB010"] <- "bolditalic('*010')"
  x[x == "ApcoIIB011"] <- "bolditalic('*011')"
  x[x == "ApcoIIB012"] <- "bolditalic('*012')"
  x[x == "ApcoIIB013"] <- "bolditalic('*013')"
  x[x == "ApcoIIB014"] <- "bolditalic('*014')"
  x[x == "ApcoIIB015"] <- "bolditalic('*015')"
  x[x == "ApcoIIB017"] <- "bolditalic('*017')"
  x[x == "ApcoIIB018"] <- "bolditalic('*018')"
  x[x == "ApcoIIB019"] <- "bolditalic('*019')"
  x[x == "ApcoIIB020"] <- "bolditalic('*020')"
  x[x == "ApcoIIB021"] <- "bolditalic('*021')"
  x[x == "ApcoIIB023"] <- "bolditalic('*023')"
  x[x == "ApcoIIB024"] <- "bolditalic('*024')"
  x[x == "ApcoIIB025"] <- "bolditalic('*025')"
  x[x == "ApcoIIB026"] <- "bolditalic('*026')"
  x[x == "ApcoIIB027"] <- "bolditalic('*027')"
  x[x == "ApcoIIB028"] <- "bolditalic('*028')"
  x[x == "ApcoIIB029"] <- "bolditalic('*029')"
  x[x == "ApcoIIB030"] <- "bolditalic('*030')"
  x[x == "ApcoIIB031"] <- "bolditalic('*031')"
  x[x == "ApcoIIB032"] <- "bolditalic('*032')"
  x[x == "ApcoIIB033"] <- "bolditalic('*033')"
  x[x == "ApcoIIB034"] <- "bolditalic('*034')"
  x[x == "ApcoIIB035"] <- "bolditalic('*035')"
  x[x == "ApcoIIB036"] <- "bolditalic('*036')"
  x[x == "ApcoIIB037"] <- "bolditalic('*037')"
  x[x == "ApcoIIB038"] <- "bolditalic('*038')"
  return(x)
}

my_colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(38)

# Original genetic cluster plot (without Admixed)
plot_cluster <- ggplot(long_freq, aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = 'stack', width = 0.8) +
  facet_wrap(~Location, ncol = 1, strip.position = "right") +
  scale_y_continuous(limits = c(0, 0.27), breaks = seq(0, 0.30, by = 0.1), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.015, 0.015), labels = function(x) parse(text = custom_labels(x))) +
  scale_fill_manual(values = my_colors) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 12),
    axis.text.y = element_text(size = 12),
    axis.ticks.y = element_line(size = 0.36),
    strip.text.y = element_text(size = 14, color = "black", angle = 0, hjust = 0.5, vjust = 0.5),
    panel.spacing.y = grid::unit(0, "lines"),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 14, face = "bold", margin = margin(r = 10)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Alleles", y = "Frequency") +
  guides(fill = "none")

print(plot_cluster)

ggsave(file.path(output_dir, "cluster_histogram.png"), plot = plot_cluster, dpi = 300)

# Count clusters containing each allele
presence <- data.frame(
  Allele = all_alleles,
  ClustersPresent = vapply(all_alleles, function(a) sum(cluster_data[[a]] > 0), integer(1)),
  Genetic_clusters = vapply(all_alleles, function(a) {
    paste(as.character(cluster_data$Genetic_cluster[cluster_data[[a]] > 0]), collapse = ", ")
  }, character(1))
)
write.csv(presence, file.path(output_dir, "MHC_allele_cluster_presence.csv"), row.names = FALSE)
