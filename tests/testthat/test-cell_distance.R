library(testthat)
library(SummarizedExperiment)
library(BiocParallel)

set.seed(2025)
sce <- simulate_branched_random_walk(N_gene=12, N_path=3, length_path=5)
sce <- suppressWarnings(Sanity(sce))

# calculateSNR tests --------------------------------------------------------

test_that("calculateSNR returns a named numeric vector of correct length", {
    snr <- calculateSNR(sce)
    expect_true(is.numeric(snr))
    expect_equal(length(snr), nrow(sce))
    expect_equal(names(snr), rownames(sce))
})

test_that("calculateSNR values are non-negative", {
    snr <- calculateSNR(sce)
    expect_true(all(snr >= 0))
})

test_that("calculateSNR respects subset.row", {
    snr_all <- calculateSNR(sce)
    snr_sub <- calculateSNR(sce, subset.row=1:6)
    expect_equal(length(snr_sub), 6L)
    expect_equal(snr_sub, snr_all[1:6])
})

test_that("calculateSNR warns on negative variances", {
    obj <- sce
    assay(obj, "logcounts_sd")[,] <- 0
    snr <- expect_warning(calculateSNR(obj), "negative activity variances")
    expect_true(all(is.finite(snr)))
    expect_true(all(snr >= 0))
})

test_that("calculateSNR SNR matches the filter used in calculateSanityDistance", {
    snr <- calculateSNR(sce)
    snr_cutoff <- 1
    keep <- which(snr >= snr_cutoff)

    # Distance using subset.row based on SNR should match using snr_cutoff
    d_snr <- calculateSanityDistance(sce, snr_cutoff=snr_cutoff, nbin=10L,
                                     BPPARAM=SerialParam())
    d_sub <- calculateSanityDistance(sce, snr_cutoff=0, subset.row=keep,
                                     nbin=10L, BPPARAM=SerialParam())
    expect_equal(as.numeric(d_snr), as.numeric(d_sub))
})

# calculateSanityDistance tests ---------------------------------------------

test_that("calculateSanityDistance returns a valid dist object", {
    dist_obj <- calculateSanityDistance(
        sce,
        assay="logcounts",
        assay.sd="logcounts_sd",
        gene_sd="sanity_activity_sd",
        gene_mu="sanity_log_activity_mean",
        mu_sd="sanity_log_activity_mean_sd",
        snr_cutoff=0,
        nbin=10L,
        BPPARAM=SerialParam()
    )

    expect_true(inherits(dist_obj, "dist"))
    expect_equal(attr(dist_obj, "Size"), ncol(sce))
    expect_equal(length(dist_obj), choose(ncol(sce), 2))
})

test_that("calculateSanityDistance works with subset.row", {
    # Select only the first 6 genes.
    subset_rows <- 1:6
    dist_obj_subset <- calculateSanityDistance(
        sce,
        assay="logcounts",
        assay.sd="logcounts_sd",
        gene_sd="sanity_activity_sd",
        gene_mu="sanity_log_activity_mean",
        mu_sd="sanity_log_activity_mean_sd",
        snr_cutoff=0.0,
        nbin=10L,
        subset.row=subset_rows,
        BPPARAM=SerialParam()
    )

    expect_true(inherits(dist_obj_subset, "dist"))
    expect_equal(attr(dist_obj_subset, "Size"), ncol(sce))
    expect_equal(length(dist_obj_subset), choose(ncol(sce), 2))

    # Compare with full data calculation.
    dist_obj_all <- calculateSanityDistance(
        sce,
        assay="logcounts",
        assay.sd="logcounts_sd",
        gene_sd="sanity_activity_sd",
        gene_mu="sanity_log_activity_mean",
        mu_sd="sanity_log_activity_mean_sd",
        snr_cutoff=0.0,
        nbin=10L,
        BPPARAM=bpparam()
    )
    expect_false(identical(as.numeric(dist_obj_subset), as.numeric(dist_obj_all)))
})

test_that("calculateSanityDistance warns on negative variances", {
    obj <- sce
    assay(obj, "logcounts_sd")[,] <- 0
    expect_warning(
        calculateSanityDistance(
            obj,
            assay="logcounts",
            assay.sd="logcounts_sd",
            gene_sd="sanity_activity_sd",
            gene_mu="sanity_log_activity_mean",
            mu_sd="sanity_log_activity_mean_sd",
            snr_cutoff=0,
            nbin=5L,
            BPPARAM=SerialParam()
        ),
        "negative activity variances"
    )
})

test_that("calculateSanityDistance errors with negative snr", {
    expect_error(calculateSanityDistance(sce, snr_cutoff=-1),
                 "SNR cutoff cannot be negative")
})

test_that("calculateSanityDistance errors with missing gene variance", {
    bad <- sce
    rowData(bad)$sanity_activity_sd[1] <- NA_real_
    expect_error(
        calculateSanityDistance(bad, snr_cutoff=0, nbin=5L),
        "Gene variance cannot be negative"
    )
})

test_that("calculateSanityDistance returns zero for identical cells", {
    tmp <- sce[,c(1L, 1L)]
    mu_sd <- "sanity_log_activity_mean_sd"
    rowData(tmp)[[mu_sd]] <- assay(sce, "logcounts_sd")[, 1L]
    d <- calculateSanityDistance(tmp, snr_cutoff=0, nbin=5L)
    expect_equal(as.numeric(d), rep(0, length(d)))
})
