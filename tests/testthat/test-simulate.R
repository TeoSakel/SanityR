library(testthat)
library(SingleCellExperiment)

set.seed(42)

test_that("simulate_independent_cells returns a valid SCE", {
    sce <- simulate_independent_cells(N_cell = 20, N_gene = 15)

    # Check that the returned object is a SingleCellExperiment
    expect_s4_class(sce, "SingleCellExperiment")

    # Check that required assays exist and have correct dimensions
    expect_true(all(c("counts", "logFC") %in% assayNames(sce)))
    expect_true(all(rowSums(counts(sce)) > 0))
    expect_lte(nrow(sce), 15)
    expect_equal(ncol(sce), 20)

    # Check that rowData and colData contain the expected columns
    expect_true(all(c("ltq_mean", "ltq_var") %in% colnames(rowData(sce))))
    expect_true(all("cell_size" %in% colnames(colData(sce))))
})

test_that("simulate_branched_random_walk returns a valid SCE with lineage graph", {
    sce <- simulate_branched_random_walk(N_gene = 10, N_path = 3, length_path = 5)

    # Check that the returned object is a SingleCellExperiment
    expect_s4_class(sce, "SingleCellExperiment")

    # Check that required assays exist
    expect_true(all(c("counts", "logFC") %in% assayNames(sce)))
    expect_true(all(rowSums(counts(sce)) > 0))
    expect_lte(nrow(sce), 10)
    expect_equal(ncol(sce), 3 * 5)

    # Check that colData contains neighbor info
    expect_true(all(c("cell_size", "predecessor") %in% colnames(colData(sce))))
})

test_that("simulate_independent_cells uses provided sizes", {
    cell_size <- c(100, 200, 300)
    N_cell <- 100L
    gene_size <- c(6, 60)
    N_gene <- 30L  # keep it large to ensure both are selected by chance
    sce <- simulate_independent_cells(
        cell_size=cell_size,
        gene_size=gene_size,
        N_cell=N_cell,
        N_gene=N_gene
    )
    expect_true(all(colData(sce)$cell_size %in% cell_size))
    expect_equal(ncol(sce), N_cell)
    ltq_mean <- rowData(sce)$ltq_mean
    expect_true(length(unique(ltq_mean)) == length(gene_size))
    expect_equal(nrow(sce), N_gene)
})

test_that("simulate_independent_cells works with single cell and gene", {
    sce <- simulate_independent_cells(N_cell=1, N_gene=1)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(dim(counts(sce)), c(1L, 1L))
})

test_that("simulate_branched_random_walk respects path settings", {
    sce <- simulate_branched_random_walk(N_path=2, length_path=3, N_gene=4)
    expect_equal(ncol(sce), 6)
    expect_lte(nrow(sce), 4)
    expect_true("predecessor" %in% colnames(colData(sce)))
})

test_that("simulate_branched_random_walk handles single step", {
    sce <- simulate_branched_random_walk(N_path=1, length_path=1, N_gene=2)
    expect_equal(ncol(sce), 1)
    expect_true(all(colData(sce)$predecessor %in% 0L))
})
