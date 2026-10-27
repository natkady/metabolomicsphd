# Validation against GEO2R itself, on the dataset used in the thesis (GSE148426).
#
# This test is skipped until a GEO2R results file is provided. To create it:
#   1. Open https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE148426
#   2. Define one group per level of 'estrogen receptor status' (Negative, Positive, Unknown,
#      Weak Positive / Borderline) and assign every sample to its group.
#   3. Keep the default options and click 'Analyze'.
#   4. Click 'Download full table' and save the file in this folder as GEO2R_GSE148426.tsv
#
# The test downloads the series from GEO (internet needed, roughly 3.5 GB of RAM), runs the app's
# functions with the same grouping, and checks that the F statistics and adjusted p-values agree.

test_that("results agree with GEO2R for GSE148426 (estrogen receptor status)", {
  ref_file <- test_path("GEO2R_GSE148426.tsv")
  skip_if_not(file.exists(ref_file), "GEO2R results file not provided (see instructions at top of file)")
  skip_if_offline()

  Sys.setenv(VROOM_CONNECTION_SIZE = 10000001)
  gse <- GEOquery::getGEO("GSE148426", GSEMatrix = TRUE, AnnotGPL = FALSE)[["GSE148426_series_matrix.txt.gz"]]
  pheno <- diff_cell(gse@phenoData@data)
  assay <- as.data.frame(gse@assayData[["exprs"]])
  res <- top_gene_exp(pheno, as.data.frame(t(assay)), assay, "estrogen receptor status")

  ref <- read.delim(ref_file, check.names = FALSE)
  ids <- intersect(as.character(ref$ID), rownames(res$tT))
  expect_gt(length(ids), 0.9 * nrow(ref))
  ref <- ref[match(ids, ref$ID), ]
  expect_equal(res$tT[ids, "F"], ref$F, tolerance = 1e-3)
  expect_equal(res$tT[ids, "adj.P.Val"], ref$adj.P.Val, tolerance = 1e-3)
})
