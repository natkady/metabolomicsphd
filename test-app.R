# Tests of the Shiny app's server logic, using simulated data in place of a GEO download.

load_app <- function() {
  app_env <- new.env(parent = globalenv())
  withr::with_dir(app_dir, sys.source("geo_analysis.R", envir = app_env))
  app_env$getGEO <- fake_getGEO  # no internet needed: the app reads the simulated series
  app_env
}

test_that("loading a dataset shows its phenotype data and records it for the About tab", {
  app_env <- load_app()
  testServer(app_env$server, {
    session$setInputs(geoacc = "GSE000001", button = 1)
    expect_equal(nrow(phenotab()), 12)
    expect_equal(geo_info()$accession, "GSE000001")
    about <- about_text()
    expect_true(any(grepl("GSE000001", about)))
    expect_true(any(grepl(app_env$app_version, about, fixed = TRUE)))
    expect_true(any(grepl("Most recent sample update date in GEO: Jan 01 2024", about, fixed = TRUE)))
  })
})

test_that("running the analysis produces results, and loading new data clears them", {
  app_env <- load_app()
  testServer(app_env$server, {
    session$setInputs(geoacc = "GSE000001", button = 1)
    session$setInputs(meta = "group", meta_button = 1)
    session$setInputs(ctschoice = "control-treated")
    hist_before <- output$genehist$src
    expect_false(is.null(jsonlite::fromJSON(output$genetables)$x))
    expect_false(is.null(output$genects))

    session$setInputs(geoacc = "GSE000002", button = 2)
    expect_lt(nchar(output$genehist$src), nchar(hist_before) / 2)  # now an empty plot area
    expect_null(jsonlite::fromJSON(output$genetables)$x)          # table emptied
    expect_null(output$genects)                                   # contrast choice removed
    expect_equal(basename(output$histsave), "No_Results.txt")     # old downloads no longer available
    expect_equal(geo_info()$accession, "GSE000002")
  })
})

test_that("the pathway tab explains the problem instead of failing when Java is missing", {
  app_env <- load_app()
  app_env$java_available <- FALSE
  testServer(app_env$server, {
    session$setInputs(geoacc = "GSE000001", button = 1)
    session$setInputs(meta = "group", meta_button = 1)
    expect_match(as.character(output$geneset$html), "Java was not found")
  })
})
