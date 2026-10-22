# Libraries
library(binom)   # Wilson confidence intervals for prevalence
library(dplyr)   # Data filtering and summaries
library(ggplot2) # Plots
library(ggpubr)  # Combines plots
library(see)     # Half-violin plot

# Input and output directories
working_dir <- getwd()
input_dir <- file.path(working_dir, "r_input")
output_dir <- file.path(working_dir, "r_output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Cluster order and colors
custom_order <- c("Merritt Island", "SL", "CBC", "SBC", "Archbold", "Ocala")
unit_colors <- c(
  "Merritt Island" = "#d62728",
  "SL" = "royalblue1",
  "CBC" = "navy",
  "SBC" = "chartreuse3",
  "Archbold" = "mediumorchid4",
  "Ocala" = "darkorange"
)

# Read data and keep birds with disease incidence information
data <- read.csv(file.path(input_dir, "MHC_sample_alleles.csv"), stringsAsFactors = FALSE)
data <- data %>%
  filter(!is.na(Incidence), trimws(as.character(Incidence)) != "",
         Genetic_cluster %in% custom_order) %>%
  mutate(
    Incidence = as.numeric(Incidence),
    Intensity = as.numeric(Intensity),
    logIntensity = as.numeric(logIntensity)
  ) %>%
  filter(Incidence %in% c(0, 1))

# Calculate prevalence and Wilson 95% confidence intervals
prevalence_df <- data %>%
  group_by(Genetic_cluster) %>%
  summarise(
    SampleSize = n(),
    Infected = sum(Incidence == 1),
    Prevalence = Infected / SampleSize,
    .groups = "drop"
  )
ci <- binom.confint(prevalence_df$Infected, prevalence_df$SampleSize,
                    methods = "wilson")
prevalence_df$LowerCI <- ci$lower
prevalence_df$UpperCI <- ci$upper
prevalence_df$Cluster <- factor(prevalence_df$Genetic_cluster, levels = custom_order)

# Summarize intensity among infected birds
withmalaria_data <- data %>% filter(Incidence == 1)
summary_data <- withmalaria_data %>%
  group_by(Genetic_cluster) %>%
  summarise(
    MeanIntensity = mean(Intensity, na.rm = TRUE),
    Mean = mean(logIntensity, na.rm = TRUE),
    N_log = sum(!is.na(logIntensity)),
    SD_log = sd(logIntensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    MarginError = qt(0.975, df = N_log - 1) * SD_log / sqrt(N_log),
    LowerCI = Mean - MarginError,
    UpperCI = Mean + MarginError
  )

# Prepare data for the half-violin plot
plot_data <- withmalaria_data %>%
  left_join(summary_data, by = "Genetic_cluster") %>%
  mutate(Genetic_cluster = factor(Genetic_cluster, levels = custom_order))

# Plot infection prevalence
prevalenceplot <- ggplot(
  prevalence_df,
  aes(x = Cluster, y = Prevalence, fill = Cluster)
) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = expression("Prevalence of " * italic("Plasmodium") * " Infection by Genetic Cluster"),
    y = "Prevalence"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

# Plot log-intensity among infected birds
violinplot <- ggplot(
  plot_data,
  aes(Genetic_cluster, logIntensity, fill = Genetic_cluster)
) +
  geom_violinhalf(side = "left", trim = FALSE, color = "black") +
  geom_point(aes(y = Mean), color = "white", size = 3, shape = 18) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, color = "black") +
  scale_fill_manual(values = unit_colors) +
  labs(
    title = "Log-transformed Infection Intensity by Genetic Cluster",
    y = "Log-transformed Intensity",
    x = "Genetic Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

# Combine and save the two-panel figure
combined_plot <- ggarrange(prevalenceplot, violinplot,
                           labels = c("a", "b"), ncol = 1, nrow = 2)
ggsave(file.path(output_dir, "prevalence_logIntensity.jpg"),
       plot = combined_plot, device = "jpeg", width = 10, height = 12, dpi = 300)

# Save cluster-level disease summary
summary_table <- prevalence_df %>%
  select(Genetic_cluster, SampleSize, Infected, Prevalence,
         Prevalence_LowerCI = LowerCI, Prevalence_UpperCI = UpperCI) %>%
  left_join(
    summary_data %>%
      select(Genetic_cluster, MeanIntensity, MeanlogIntensity = Mean,
             logIntensity_LowerCI = LowerCI, logIntensity_UpperCI = UpperCI),
    by = "Genetic_cluster"
  ) %>%
  mutate(Genetic_cluster = factor(Genetic_cluster, levels = custom_order)) %>%
  arrange(Genetic_cluster)
write.csv(summary_table,
          file.path(output_dir, "MHC_disease_summary_by_cluster.csv"),
          row.names = FALSE)
