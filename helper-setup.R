# Shared set-up for the tests. testthat runs files starting with "helper" before the tests.
# Run all tests from the app folder with:  testthat::test_dir("tests/testthat")

suppressPackageStartupMessages({
  library(testthat)
  library(shiny)
  library(tidyverse)
  library(limma)
  library(DescTools)
  library(Biobase)
  library(data.table)
})

app_dir <-"[directory pathway where GEO files have been saved eg. C:/Users/natkady]"
source(file.path(app_dir, "diff_cell.R"))

# A small simulated GEO series with known answers: 200 genes x 12 samples on the log2 scale.
# Genes 1-10 are 3 log2 units higher in the "treated" group. The treated samples come first,
# so the order the groups appear in differs from alphabetical order (control, treated).
make_fake_series <- function(n_genes = 200, groups = rep(c("treated", "control"), each = 6),
                             de_group = "treated", seed = 1) {
  set.seed(seed)
  n <- length(groups)
  ex <- matrix(rnorm(n_genes * n, mean = 8, sd = 0.5), nrow = n_genes,
               dimnames = list(paste0("gene", seq_len(n_genes)), sprintf("GSM%03d", seq_len(n))))
  ex[1:10, groups == de_group] <- ex[1:10, groups == de_group] + 3
  pheno <- data.frame(title = paste("sample", seq_len(n)),
                      geo_accession = colnames(ex),
                      status = "Public on Jan 01 2024",
                      source_name_ch1 = "tissue",
                      organism_ch1 = "Homo sapiens",
                      characteristics_ch1 = paste("group:", groups),
                      last_update_date = "Jan 01 2024",
                      `group:ch1` = groups,
                      `age:ch1` = as.character(seq(30, by = 2, length.out = n)),
                      check.names = FALSE, row.names = colnames(ex))
  list(ex = ex, pheno = pheno, groups = groups)
}

# Wrap the simulated series the way GEOquery::getGEO() returns it
fake_getGEO <- function(GEO, ...) {
  d <- make_fake_series()
  eset <- Biobase::ExpressionSet(assayData = d$ex, phenoData = Biobase::AnnotatedDataFrame(d$pheno))
  setNames(list(eset), paste0(GEO, "_series_matrix.txt.gz"))
}

# Reference analysis written directly from the GEO2R R script (limma), independent of the app code
geo2r_reference <- function(ex, groups) {
  gs <- factor(make.names(groups))
  design <- model.matrix(~ gs + 0)
  colnames(design) <- levels(gs)
  fit <- lmFit(ex, design)
  cts <- apply(combn(levels(gs), 2), 2, paste, collapse = "-")
  fit2 <- eBayes(contrasts.fit(fit, makeContrasts(contrasts = cts, levels = design)), 0.01)
  list(cts = cts, fit2 = fit2, tT = topTable(fit2, adjust.method = "fdr", sort.by = "B", number = Inf))
}
