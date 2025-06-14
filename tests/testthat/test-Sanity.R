library(testthat)
suppressPackageStartupMessages(library(SingleCellExperiment))

set.seed(2025)
sce <- simulate_independent_cells(N_cell = 100L, N_gene = 10L)
sf <- colSums(counts(sce))
sizeFactors(sce) <- sf/mean(sf)


# Test Errors & Warnings --------------------------------------------------

test_that("Sanity() enforces positive vmin", {
    expect_error(Sanity(counts(sce), vmin = 0),
                 regexp = "Minimum variance must be positive")
})

test_that("Sanity() enforces vmax > vmin", {
    expect_error(Sanity(counts(sce), vmin = 0.5, vmax = 0.125),
                 regexp = "Maximum variance must be greater than the minimum")
})

test_that("Sanity() enforces non-negative Gamma parameters", {
    expect_error(Sanity(counts(sce), a = -1, b = 1),
                 regexp = "Both Gamma parameters must be positive")
    expect_error(Sanity(counts(sce), a = 1, b = -1),
                 regexp = "Both Gamma parameters must be positive")
})


# Test Form of Results ----------------------------------------------------

test_that("Sanity() returns a list with expected components for a matrix input", {
    res <- Sanity(counts(sce))
    expected_names <- c("mu", "var_mu", "var", "delta", "var_delta", "likelihood")
    expect_true(is.list(res))
    expect_true(all(expected_names %in% names(res)))
    expect_equal(length(res[["mu"]]), nrow(counts(sce)))
    expect_equal(length(res[["var_mu"]]), nrow(counts(sce)))
    expect_equal(dim(res[["delta"]]), dim(counts(sce)))
    expect_equal(dim(res[["var_delta"]]), dim(counts(sce)))
    expect_equal(dim(res[["likelihood"]]), c(nrow(counts(sce)), 160L))
})

test_that("Sanity() output components are numeric and finite", {
  res <- Sanity(counts(sce))
  for (val in seq_along(res)) {
      expect_true(all(is.numeric(val)))
      expect_true(all(is.finite(val)))
  }
})

test_that("Sanity() method for SummarizedExperiment modifies rowData and adds assays", {
    res <- Sanity(as(sce, "SummarizedExperiment"))
    new_cols <- c("sanity_log_activity_mean", "sanity_log_activity_mean_sd", "sanity_activity_sd")
    expect_true(all(new_cols %in% colnames(rowData(res))))
    expect_true("logcounts" %in% assayNames(res))
    expect_true("logcounts_sd" %in% assayNames(res))
})

test_that("Sanity() names the assays correctly", {
    res <- Sanity(as(sce, "SummarizedExperiment"), name = "sanity")
    expect_equal(assayNames(res), c(assayNames(sce), "sanity", "sanity_sd"))
})


# Test Logic of Results ---------------------------------------------------

test_that("Sanity() variances within the given limits", {
    res <- Sanity(counts(sce), vmin = 0.01, vmax = 1.5, nbin = 20)
    expect_true(all(res[["var"]] >= 0.01))
    expect_true(all(res[["var"]] <= 1.5))
    expect_true(all(res[["var_mu"]] >= 0))
    expect_true(all(res[["var_delta"]] >= 0))
})

test_that("Sanity() returns a valid posterior for gene variance", {
    res <- Sanity(counts(sce), vmin = 0.01, vmax = 1.5, nbin = 20)
    expect_true(all(res[["likelihood"]] >= 0))
    expect_equal(rowSums(res[["likelihood"]]), rep(1, nrow(sce), ignore_attr = TRUE))
})

test_that("Sanity() produced LTQ that some to 1", {
    # requires size factors to work so use the SCE object
    res <- sce |> Sanity() |> logcounts() |> exp() |> colSums()
    expect_equal(res, rep(1, length(res)), tolerance = .01, ignore_attr = TRUE)
})

test_that("Sanity() returns nearly zero cell-specific differences for uniform input", {
    mu <- exp(rowData(sce)$ltq_mean)
    sf <- colData(sce)[["cell_size"]]
    rate <- outer(mu, sf, FUN = "*")
    mat <- matrix(rpois(n = prod(dim(rate)), lambda = rate), nrow = nrow(sce))
    res <- Sanity(mat, size.factors = sf / mean(sf), vmin = 1e-5, vmax = 10)
    p <- with(res, abs(delta) / sqrt(var_delta)) |>
        pnorm(lower.tail = FALSE) |>
        p.adjust()
    expect_true(min(p) > 0.05)
})
