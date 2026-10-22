# Libraries
library(dplyr)  # Data preparation.
library(lme4)   # Mixed-effects models.
library(car)    # Variance inflation factors.
library(MuMIn)  # AIC model comparison.
library(pROC)   # ROC curves and AUC.

# Input and output directories
working_dir <- getwd()
input_dir <- file.path(working_dir, "r_input")
output_dir <- file.path(working_dir, "r_output", "incidence_models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read and filter disease data
data <- read.csv(file.path(input_dir, "MHC_sample_alleles.csv"),
                 stringsAsFactors = FALSE, check.names = FALSE)
data <- data %>%
  filter(!is.na(Incidence), Incidence %in% c(0, 1),
         !is.na(Genetic_cluster), Genetic_cluster != "", Genetic_cluster != "Admixed")

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
screen_sets <- list(Full = c("Month", "Genetic_cluster", "Nucleotide_Diversity", "Genotype", "Functional_Alleles", "Supertype_count", "SupertypeIIC", "ApcoIIB001", "ApcoIIB003", "ApcoIIB004", "ApcoIIB005", "ApcoIIB006", "ApcoIIB007", "ApcoIIB008", "ApcoIIB009", "ApcoIIB010", "ApcoIIB011", "ApcoIIB012", "ApcoIIB013", "ApcoIIB014", "ApcoIIB015", "ApcoIIB016", "ApcoIIB018", "ApcoIIB020", "ApcoIIB021", "ApcoIIB024", "ApcoIIB027", "ApcoIIB028"),
                    Reduced = c("Month", "Genetic_cluster", "Nucleotide_Diversity", "SupertypeIIC", "ApcoIIB001", "ApcoIIB003", "ApcoIIB004", "ApcoIIB005", "ApcoIIB006", "ApcoIIB007", "ApcoIIB008", "ApcoIIB009", "ApcoIIB010", "ApcoIIB011", "ApcoIIB012", "ApcoIIB013", "ApcoIIB014", "ApcoIIB015", "ApcoIIB016", "ApcoIIB018", "ApcoIIB020", "ApcoIIB021", "ApcoIIB024", "ApcoIIB027", "ApcoIIB028"),
                    Count_Genotype = c("Month", "Genetic_cluster", "Genotype", "Functional_Alleles", "Supertype_count"))
vif_rows <- list()
for (group in names(screen_sets)) {
  terms <- intersect(screen_sets[[group]], names(data))
  terms <- terms[vapply(data[terms], function(x) length(unique(na.omit(x))) > 1, logical(1))]
  if (!length(terms)) next
  formula <- reformulate(terms, response = "Incidence")
  screened <- tryCatch(glm(formula, data = data, family = binomial(link = "logit")), error = function(e) NULL)
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
  mod1 = c("Month"),
  mod2 = c("Genetic_cluster"),
  mod3 = c("Genotype"),
  mod4 = c("Functional_Alleles"),
  mod5 = c("Supertype_count"),
  mod6 = c("Month", "Genetic_cluster"),
  mod7 = c("Month", "Genotype"),
  mod8 = c("Month", "Functional_Alleles"),
  mod9 = c("Month", "Supertype_count"),
  mod10 = c("Genetic_cluster", "Genotype"),
  mod11 = c("Genetic_cluster", "Functional_Alleles"),
  mod12 = c("Genetic_cluster", "Supertype_count"),
  mod13 = c("Month", "Genetic_cluster", "Genotype"),
  mod14 = c("Month", "Genetic_cluster", "Functional_Alleles"),
  mod15 = c("Month", "Genetic_cluster", "Supertype_count"),
  mod16 = c("Nucleotide_Diversity"),
  mod17 = c("SupertypeIIC"),
  mod18 = c("Month", "Nucleotide_Diversity"),
  mod19 = c("Month", "SupertypeIIC"),
  mod20 = c("Genetic_cluster", "Nucleotide_Diversity"),
  mod21 = c("Genetic_cluster", "SupertypeIIC"),
  mod22 = c("Nucleotide_Diversity", "SupertypeIIC"),
  mod23 = c("Month", "Genetic_cluster", "Nucleotide_Diversity"),
  mod24 = c("Month", "Genetic_cluster", "SupertypeIIC"),
  mod25 = c("Month", "Nucleotide_Diversity", "SupertypeIIC"),
  mod26 = c("Genetic_cluster", "Nucleotide_Diversity", "SupertypeIIC"),
  mod27 = c("Month", "Genetic_cluster", "Nucleotide_Diversity", "SupertypeIIC"),
  mod28 = c("ApcoIIB008"),
  mod29 = c("ApcoIIB008", "Month"),
  mod30 = c("ApcoIIB008", "Genetic_cluster"),
  mod31 = c("ApcoIIB008", "Nucleotide_Diversity"),
  mod32 = c("ApcoIIB008", "Month", "Genetic_cluster"),
  mod33 = c("ApcoIIB008", "Month", "Nucleotide_Diversity"),
  mod34 = c("ApcoIIB008", "Genetic_cluster", "Nucleotide_Diversity"),
  mod35 = c("ApcoIIB008", "Month", "Genetic_cluster", "Nucleotide_Diversity")
)

# Fit candidates with available predictors and common complete cases for AIC
available_models <- model_terms[vapply(model_terms, function(x)
  all(x %in% names(data)), logical(1))]
all_terms <- unique(unlist(available_models, use.names = FALSE))
analysis_cols <- unique(c("Incidence", "Year", all_terms))
model_data <- data[complete.cases(data[analysis_cols]), , drop = FALSE]
model_data <- droplevels(model_data)
models <- list()
for (model_name in names(available_models)) {
  terms <- available_models[[model_name]]
  if (any(vapply(model_data[terms], function(x) length(unique(x)) < 2, logical(1)))) next
  formula <- as.formula(paste("Incidence", "~", paste(c(terms, "(1|Year)"), collapse = " + ")))
  fitted <- tryCatch(lme4::glmer(formula, data = model_data, family = binomial(link = "logit")), error = function(e) NULL)
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
  coefficients$OddsRatio <- exp(coefficients$Estimate)
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
  
  pearson <- residuals(model, type = "pearson")
  rdf <- df.residual(model)
  chisq <- sum(pearson^2)
  
  diag <- data.frame(
    Model = id,
    PearsonChiSq = chisq,
    Ratio = chisq / rdf,
    DF = rdf,
    P = pchisq(chisq, df = rdf, lower.tail = FALSE),
    Singular = lme4::isSingular(model)
  )
  
  write.csv(
    diag,
    file.path(output_dir, paste0(id, "_diagnostics.csv")),
    row.names = FALSE
  )
  
  prediction <- predict(model, type = "response")
  actual <- model.response(model.frame(model))
  
  roc_obj <- tryCatch(
    pROC::roc(actual, prediction, quiet = TRUE),
    error = function(e) NULL
  )
  
  if (!is.null(roc_obj)) {
    png(file.path(output_dir, paste0(id, "_ROC.png")))
    plot(roc_obj)
    dev.off()
    
    write.csv(
      data.frame(Model = id, AUC = as.numeric(pROC::auc(roc_obj))),
      file.path(output_dir, paste0(id, "_AUC.csv")),
      row.names = FALSE
    )
  }
}

# Export summaries of best-supported models
capture.output(
  for (id in best_model_names) {
    cat("\n", id, "\n\n")
    print(summary(models[[id]]))
  },
  file = file.path(output_dir, "best_models_summary.txt")
)
