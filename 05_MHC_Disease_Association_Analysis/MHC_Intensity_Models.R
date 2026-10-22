# Libraries
library(dplyr)  # Data preparation.
library(lme4)   # Mixed-effects models.
library(car)    # Variance inflation factors.
library(MuMIn)  # AIC model comparison.
library(lmtest) # Breusch-Pagan diagnostics.

# Input and output directories
working_dir <- getwd()
input_dir <- file.path(working_dir, "r_input")
output_dir <- file.path(working_dir, "r_output", "intensity_models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read and filter disease data
data <- read.csv(file.path(input_dir, "MHC_sample_alleles.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
data <- data %>%
  filter(!is.na(Incidence), Incidence %in% c(0, 1),
         !is.na(Genetic_cluster), Genetic_cluster != "", Genetic_cluster != "Admixed")
data <- data %>% filter(Incidence == 1, is.finite(as.numeric(logIntensity)))

# Prepare predictors
data$Genotype <- as.character(data$Genotype)
genotype_counts <- table(data$Genotype)
rare_genotypes <- names(genotype_counts)[genotype_counts < 5]
data$Genotype[data$Genotype %in% rare_genotypes] <- "Rare"
for (v in intersect(c("Genetic_cluster", "Month", "Year", "Genotype"), names(data))) {
  data[[v]] <- factor(data[[v]])
}
allele_cols <- grep("^ApcoIIB[0-9]{3}$", names(data), value = TRUE)
supertype_cols <- grep("^SupertypeII", names(data), value = TRUE)
for (v in c(allele_cols, supertype_cols)) data[[v]] <- as.numeric(as.character(data[[v]]))
for (v in intersect(c("Functional_Alleles", "Nucleotide_Diversity", "Supertype_count",
                      "logIntensity"), names(data))) {
  data[[v]] <- as.numeric(data[[v]])
}

# Allele and supertype counts
allele_counts <- data.frame(Allele = allele_cols,
                            Carriers = colSums(data[allele_cols] == 1, na.rm = TRUE))
write.csv(allele_counts, file.path(output_dir, "allele_carrier_counts.csv"), row.names = FALSE)
if (length(supertype_cols) > 0) {
  supertype_counts <- data.frame(Supertype = supertype_cols,
                                 Carriers = colSums(data[supertype_cols] == 1, na.rm = TRUE))
  write.csv(supertype_counts, file.path(output_dir, "supertype_carrier_counts.csv"), row.names = FALSE)
}

# Predictor correlations
numeric_cols <- names(data)[vapply(data, is.numeric, logical(1))]
numeric_cols <- numeric_cols[vapply(data[numeric_cols], function(x)
  sum(is.finite(x)) >= 3 && sd(x, na.rm = TRUE) > 0, logical(1))]
cor_matrix <- cor(data[numeric_cols], use = "pairwise.complete.obs")
write.csv(cor_matrix, file.path(output_dir, "correlation_matrix.csv"))

# VIF screening: original full set, reduced set, and count/genotype set
screen_sets <- list(Full = c("Month", "Genetic_cluster", "Nucleotide_Diversity", "Genotype", "Functional_Alleles", "Supertype_count", "SupertypeIIC", "ApcoIIB001", "ApcoIIB003", "ApcoIIB004", "ApcoIIB005", "ApcoIIB006", "ApcoIIB007", "ApcoIIB008", "ApcoIIB009", "ApcoIIB010", "ApcoIIB011", "ApcoIIB012", "ApcoIIB016"),
                    Reduced = c("Month", "Genetic_cluster", "Nucleotide_Diversity", "SupertypeIIC", "ApcoIIB001", "ApcoIIB003", "ApcoIIB004", "ApcoIIB005", "ApcoIIB006", "ApcoIIB007", "ApcoIIB008", "ApcoIIB009", "ApcoIIB010", "ApcoIIB011", "ApcoIIB012", "ApcoIIB016"),
                    Count_Genotype = c("Month", "Genetic_cluster", "Genotype", "Functional_Alleles", "Supertype_count"))
vif_rows <- list()
for (group in names(screen_sets)) {
  terms <- intersect(screen_sets[[group]], names(data))
  terms <- terms[vapply(data[terms], function(x) length(unique(na.omit(x))) > 1, logical(1))]
  if (!length(terms)) next
  formula <- reformulate(terms, response = "logIntensity")
  screened <- tryCatch(glm(formula, data = data, family = gaussian()), error = function(e) NULL)
  if (is.null(screened)) next
  scores <- tryCatch(car::vif(screened), error = function(e) NULL)
  if (is.null(scores)) next
  if (is.matrix(scores)) {
    vals <- scores[, if ("GVIF^(1/(2*Df))" %in% colnames(scores)) "GVIF^(1/(2*Df))" else 1]
  } else vals <- scores
  vif_rows[[group]] <- data.frame(Group = group, Predictor = names(vals), VIF = as.numeric(vals))
}
if (length(vif_rows)) write.csv(do.call(rbind, vif_rows),
  file.path(output_dir, "vif_results.csv"), row.names = FALSE)

# Original candidate models (Year random intercept)
model_terms <- list(
  mod1 = c("Genetic_cluster"),
  mod2 = c("Month"),
  mod3 = c("Nucleotide_Diversity"),
  mod4 = c("SupertypeIIC"),
  mod5 = c("Genetic_cluster", "Month"),
  mod6 = c("Genetic_cluster", "Nucleotide_Diversity"),
  mod7 = c("Genetic_cluster", "SupertypeIIC"),
  mod8 = c("Month", "Nucleotide_Diversity"),
  mod9 = c("Month", "SupertypeIIC"),
  mod10 = c("Nucleotide_Diversity", "SupertypeIIC"),
  mod11 = c("Genetic_cluster", "Month", "Nucleotide_Diversity"),
  mod12 = c("Genetic_cluster", "Month", "SupertypeIIC"),
  mod13 = c("Genetic_cluster", "Nucleotide_Diversity", "SupertypeIIC"),
  mod14 = c("Month", "Nucleotide_Diversity", "SupertypeIIC"),
  mod15 = c("Genetic_cluster", "Month", "Nucleotide_Diversity", "SupertypeIIC"),
  mod16 = c("Genotype"),
  mod17 = c("Functional_Alleles"),
  mod18 = c("Supertype_count"),
  mod19 = c("Month", "Genotype"),
  mod20 = c("Month", "Functional_Alleles"),
  mod21 = c("Month", "Supertype_count"),
  mod22 = c("Genetic_cluster", "Genotype"),
  mod23 = c("Genetic_cluster", "Functional_Alleles"),
  mod24 = c("Genetic_cluster", "Supertype_count"),
  mod25 = c("Month", "Genetic_cluster", "Genotype"),
  mod26 = c("Month", "Genetic_cluster", "Functional_Alleles"),
  mod27 = c("Month", "Genetic_cluster", "Supertype_count"),
  mod28 = c("ApcoIIB003"),
  mod29 = c("ApcoIIB004"),
  mod30 = c("ApcoIIB016"),
  mod31 = c("Month"),
  mod32 = c("Genetic_cluster"),
  mod33 = c("Nucleotide_Diversity"),
  mod34 = c("ApcoIIB003", "ApcoIIB004"),
  mod35 = c("ApcoIIB003", "ApcoIIB016"),
  mod36 = c("ApcoIIB004", "ApcoIIB016"),
  mod37 = c("ApcoIIB003", "Month"),
  mod38 = c("ApcoIIB003", "Genetic_cluster"),
  mod39 = c("ApcoIIB003", "Nucleotide_Diversity"),
  mod40 = c("ApcoIIB004", "Month"),
  mod41 = c("ApcoIIB004", "Genetic_cluster"),
  mod42 = c("ApcoIIB004", "Nucleotide_Diversity"),
  mod43 = c("ApcoIIB016", "Month"),
  mod44 = c("ApcoIIB016", "Genetic_cluster"),
  mod45 = c("ApcoIIB016", "Nucleotide_Diversity"),
  mod46 = c("Month", "Genetic_cluster"),
  mod47 = c("Month", "Nucleotide_Diversity"),
  mod48 = c("Genetic_cluster", "Nucleotide_Diversity"),
  mod49 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016"),
  mod50 = c("ApcoIIB003", "ApcoIIB004", "Month"),
  mod51 = c("ApcoIIB003", "ApcoIIB004", "Genetic_cluster"),
  mod52 = c("ApcoIIB003", "ApcoIIB004", "Nucleotide_Diversity"),
  mod53 = c("ApcoIIB003", "ApcoIIB016", "Month"),
  mod54 = c("ApcoIIB003", "ApcoIIB016", "Genetic_cluster"),
  mod55 = c("ApcoIIB003", "ApcoIIB016", "Nucleotide_Diversity"),
  mod56 = c("ApcoIIB004", "ApcoIIB016", "Month"),
  mod57 = c("ApcoIIB004", "ApcoIIB016", "Genetic_cluster"),
  mod58 = c("ApcoIIB004", "ApcoIIB016", "Nucleotide_Diversity"),
  mod59 = c("ApcoIIB003", "Month", "Genetic_cluster"),
  mod60 = c("ApcoIIB003", "Month", "Nucleotide_Diversity"),
  mod61 = c("ApcoIIB003", "Genetic_cluster", "Nucleotide_Diversity"),
  mod62 = c("ApcoIIB004", "Month", "Genetic_cluster"),
  mod63 = c("ApcoIIB004", "Month", "Nucleotide_Diversity"),
  mod64 = c("ApcoIIB004", "Genetic_cluster", "Nucleotide_Diversity"),
  mod65 = c("ApcoIIB016", "Month", "Genetic_cluster"),
  mod66 = c("ApcoIIB016", "Month", "Nucleotide_Diversity"),
  mod67 = c("ApcoIIB016", "Genetic_cluster", "Nucleotide_Diversity"),
  mod68 = c("Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod69 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016", "Month"),
  mod70 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016", "Genetic_cluster"),
  mod71 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016", "Nucleotide_Diversity"),
  mod72 = c("ApcoIIB003", "ApcoIIB004", "Month", "Genetic_cluster"),
  mod73 = c("ApcoIIB003", "ApcoIIB004", "Month", "Nucleotide_Diversity"),
  mod74 = c("ApcoIIB003", "ApcoIIB004", "Genetic_cluster", "Nucleotide_Diversity"),
  mod75 = c("ApcoIIB003", "ApcoIIB016", "Month", "Genetic_cluster"),
  mod76 = c("ApcoIIB003", "ApcoIIB016", "Month", "Nucleotide_Diversity"),
  mod77 = c("ApcoIIB003", "ApcoIIB016", "Genetic_cluster", "Nucleotide_Diversity"),
  mod78 = c("ApcoIIB004", "ApcoIIB016", "Month", "Genetic_cluster"),
  mod79 = c("ApcoIIB004", "ApcoIIB016", "Month", "Nucleotide_Diversity"),
  mod80 = c("ApcoIIB004", "ApcoIIB016", "Genetic_cluster", "Nucleotide_Diversity"),
  mod81 = c("ApcoIIB003", "Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod82 = c("ApcoIIB004", "Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod83 = c("ApcoIIB016", "Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod84 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016", "Month", "Genetic_cluster"),
  mod85 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016", "Month", "Nucleotide_Diversity"),
  mod86 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016", "Genetic_cluster", "Nucleotide_Diversity"),
  mod87 = c("ApcoIIB003", "ApcoIIB004", "Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod88 = c("ApcoIIB003", "ApcoIIB016", "Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod89 = c("ApcoIIB004", "ApcoIIB016", "Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod90 = c("ApcoIIB003", "ApcoIIB004", "ApcoIIB016", "Month", "Genetic_cluster", "Nucleotide_Diversity")
)

# Fit candidates with available predictors and common complete cases for AIC
available_models <- model_terms[vapply(model_terms, function(x)
  all(x %in% names(data)), logical(1))]
all_terms <- unique(unlist(available_models, use.names = FALSE))
analysis_cols <- unique(c("logIntensity", "Year", all_terms))
model_data <- data[complete.cases(data[analysis_cols]), , drop = FALSE]
model_data <- droplevels(model_data)
models <- list()
for (model_name in names(available_models)) {
  terms <- available_models[[model_name]]
  if (any(vapply(model_data[terms], function(x) length(unique(x)) < 2, logical(1)))) next
  formula <- as.formula(paste("logIntensity", "~", paste(c(terms, "(1|Year)"), collapse = " + ")))
  fitted <- tryCatch(lme4::lmer(formula, data = model_data, REML = FALSE), error = function(e) NULL)
  if (!is.null(fitted)) models[[model_name]] <- fitted
}

# Compare fitted models
if (length(models)) {
  selection <- MuMIn::model.sel(models)
  write.csv(as.data.frame(selection),
            file.path(output_dir, "model_selection_AIC.csv"), row.names = TRUE)
}

# Save available model coefficients
if (length(models)) {
  coefficients <- do.call(rbind, lapply(names(models), function(id) {
    estimate <- lme4::fixef(models[[id]])
    data.frame(Model = id, Term = names(estimate), Estimate = as.numeric(estimate))
  }))
  coefficients$FoldChange <- exp(coefficients$Estimate)
  write.csv(coefficients, file.path(output_dir, "model_coefficients.csv"), row.names = FALSE)
}

# Automatically select best-supported models (delta AICc < 2)
best_model_names <- intersect(
  rownames(selection)[!is.na(selection$delta) & selection$delta < 2],
  names(models)
)

# Diagnostics for best-supported models
for (id in best_model_names) {
  model <- models[[id]]
  
  png(file.path(output_dir, paste0(id, "_QQ.png")))
  qqnorm(residuals(model))
  qqline(residuals(model))
  dev.off()
  
  png(file.path(output_dir, paste0(id, "_residuals.png")))
  plot(fitted(model), residuals(model))
  abline(h = 0, lty = 2)
  dev.off()
  
  test <- lmtest::bptest(lm(residuals(model) ~ fitted(model)))
  
  write.csv(
    data.frame(
      Model = id,
      BPStatistic = unname(test$statistic),
      P = test$p.value,
      Singular = lme4::isSingular(model)
    ),
    file.path(output_dir, paste0(id, "_diagnostics.csv")),
    row.names = FALSE
  )
}

# Export summaries of best-supported models
capture.output(
  for (id in best_model_names) {
    cat("\n", id, "\n\n")
    print(summary(models[[id]]))
  },
  file = file.path(output_dir, "best_models_summary.txt")
)
