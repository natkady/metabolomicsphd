# Unit tests for the analysis functions in diff_cell.R, using simulated data with known answers.

test_that("diff_cell removes GEO administrative columns and sets column names and types", {
  d <- make_fake_series()
  clean <- diff_cell(d$pheno)
  expect_setequal(colnames(clean), c("group", "age"))
  expect_s3_class(clean$group, "factor")
  expect_type(clean$age, "double")
  expect_equal(rownames(clean), rownames(d$pheno))  # sample IDs kept, needed for PLS-DA
})

test_that("diff_cell works when listed administrative columns are absent from a dataset", {
  d <- make_fake_series()
  pheno <- d$pheno[, setdiff(colnames(d$pheno), c("status", "organism_ch1"))]
  expect_no_error(clean <- diff_cell(pheno))
  expect_setequal(colnames(clean), c("group", "age"))
})

test_that("two groups: contrast name and fold-change direction are correct when groups are not in alphabetical order", {
  d <- make_fake_series()
  clean <- diff_cell(d$pheno)
  res <- top_gene_exp(clean, as.data.frame(t(d$ex)), as.data.frame(d$ex), "group")
  expect_equal(res$cts, "control-treated")
  expect_null(res$tT)
  # logFC for "control-treated" must equal mean(control) - mean(treated), computed directly from the data
  direct <- rowMeans(d$ex[, d$groups == "control"]) - rowMeans(d$ex[, d$groups == "treated"])
  expect_equal(res$tT.logfc[names(direct), "logFC"], unname(direct), tolerance = 1e-8)
  expect_true(all(res$tT.logfc[paste0("gene", 1:10), "logFC"] < -2))
})

test_that("two groups: the true differentially expressed genes are found and few false positives occur", {
  d <- make_fake_series()
  res <- top_gene_exp(diff_cell(d$pheno), as.data.frame(t(d$ex)), as.data.frame(d$ex), "group")
  expect_setequal(rownames(res$tT.logfc)[1:10], paste0("gene", 1:10))
  expect_true(all(res$tT.logfc[paste0("gene", 1:10), "adj.P.Val"] < 0.05))
  expect_lt(mean(res$tT.logfc[paste0("gene", 11:200), "adj.P.Val"] < 0.05), 0.05)
})

test_that("two groups: results are identical to an analysis written from the GEO2R script", {
  d <- make_fake_series()
  res <- top_gene_exp(diff_cell(d$pheno), as.data.frame(t(d$ex)), as.data.frame(d$ex), "group")
  ref <- geo2r_reference(d$ex, d$groups)
  expect_equal(res$cts, ref$cts)
  cols <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  expect_equal(res$tT.logfc[rownames(ref$tT), cols], ref$tT[, cols], tolerance = 1e-8)
})

test_that("three groups: all pairwise contrasts are computed, each matching the GEO2R-style analysis", {
  groups <- rep(c("B", "C", "A"), each = 4)
  d <- make_fake_series(groups = groups, de_group = "C")
  res <- top_gene_exp(diff_cell(d$pheno), as.data.frame(t(d$ex)), as.data.frame(d$ex), "group")
  ref <- geo2r_reference(d$ex, groups)
  expect_equal(res$cts, c("A-B", "A-C", "B-C"))
  expect_true("F" %in% colnames(res$tT))
  expect_equal(res$tT[rownames(ref$tT), c("F", "adj.P.Val")], ref$tT[, c("F", "adj.P.Val")], tolerance = 1e-8)
  for (i in seq_along(res$cts)) {
    g <- strsplit(res$cts[i], "-")[[1]]
    direct <- rowMeans(d$ex[, groups == g[1]]) - rowMeans(d$ex[, groups == g[2]])
    expect_equal(res$tT.logfc[[i]][names(direct), "logFC"], unname(direct), tolerance = 1e-8)
  }
})

test_that("unlogged data are log2-transformed and already-logged data are left unchanged", {
  d <- make_fake_series()
  clean <- diff_cell(d$pheno)
  logged <- top_gene_exp(clean, as.data.frame(t(d$ex)), as.data.frame(d$ex), "group")
  raw <- 2^d$ex
  unlogged <- top_gene_exp(clean, as.data.frame(t(raw)), as.data.frame(raw), "group")
  expect_equal(unlogged$tT.logfc[rownames(logged$tT.logfc), "logFC"], logged$tT.logfc$logFC, tolerance = 1e-8)
})

test_that("volcano and mean difference tables contain only genes below the FDR threshold", {
  d <- make_fake_series()
  res <- top_gene_exp(diff_cell(d$pheno), as.data.frame(t(d$ex)), as.data.frame(d$ex), "group")
  for (tab in list(gene_vol_table(res$cts, res$tT.logfc, 0.05), gene_md_table(res$cts, res$tT.logfc, 0.05))) {
    expect_true(all(tab$adj.P.Val < 0.05))
    expect_true(all(tab$diffexpressed == ifelse(tab$logFC > 0, "UP", "DOWN")))
    expect_true(all(paste0("gene", 1:10) %in% rownames(tab)))
  }
})

test_that("plot functions run without error", {
  d <- make_fake_series()
  res <- top_gene_exp(diff_cell(d$pheno), as.data.frame(t(d$ex)), as.data.frame(d$ex), "group")
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(gene_hist(NULL, NULL, res$tT.logfc, res$tT))
  expect_no_error(gene_qqplot(res$fit2))
  expect_no_error(gene_densities(as.data.frame(d$ex), res$gs))
  expect_no_error(gene_meanvar(res$fit2, "GSE000001"))
  expect_s3_class(gene_volcano(res$cts, res$tT.logfc, 0.05), "ggplot")
  expect_s3_class(gene_md(res$cts, res$tT.logfc, 0.05), "ggplot")
})
