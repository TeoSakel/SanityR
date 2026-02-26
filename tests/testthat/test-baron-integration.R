## Baron et al. pancreas dataset integration test
##
## This test runs SanityR on the Baron et al. (2016) pancreas UMI count matrix
## (521 subsampled genes x 1,937 cells) and compares the results to the
## reference output from the Breda et al. (2021) paper using z-score tests.

# Load reference data --------------------------------------------------------

counts <- readRDS(system.file("extdata", "baron_UMI_counts.rds", package = "SanityR"))
ref_ltq <- readRDS(system.file("extdata", "baron_Sanity_normalized_counts.rds",
                               package = "SanityR"))
gene_data <- read.table(
    system.file("extdata", "baron_gene_data.txt.gz", package = "SanityR"),
    header = TRUE
)
cell_data <- read.table(
    system.file("extdata", "baron_cell_data.txt.gz", package = "SanityR"),
    header = TRUE
)

# Run Sanity with paper parameters -------------------------------------------
# Parameters match the original Sanity publication defaults:
#   vmin = 0.001, vmax = 50, nbin = 160, a = 1, b = 0
#   size.factors = colSums(x)  (cell library sizes)

vmin  <- 0.001
vmax  <- 50
nbin  <- 160L

res <- Sanity(counts, cell_size = cell_data$libsize,
              vmin = vmin, vmax = vmax, nbin = nbin, a = 1, b = 0)

# Reconstruct the log-uniform variance grid (mirrors C++ loop in Sanity.cpp) -
dv       <- log(vmax / vmin) / (nbin - 1L)
var_grid <- vmin * exp(dv * seq(0L, nbin - 1L))  # length nbin

# Align gene_data to counts row order:
gene_data <- gene_data[match(rownames(counts), gene_data$symbol), ]

# Test 1: Gene-level variance agreement via z-score --------------------------
# Under the Sanity model, P(v_g | n_g) is the normalized likelihood stored in
# res$likelihood (rows = genes, columns = variance bins).
# The posterior mean is E[v_g] = res$var[g].
# The posterior SD SE_var_g = sqrt(E[v_g^2] - E[v_g]^2) provides the natural
# uncertainty scale for comparing our estimates to the paper's sanity_var.

test_that("Sanity() variance estimates agree with paper on Baron dataset", {
    lik    <- res[["likelihood"]]             # G x nbin matrix
    E_v    <- drop(lik %*% var_grid)          # = res$var (posterior mean)
    E_v2   <- drop(lik %*% var_grid^2)        # second moment
    SE_var <- sqrt(pmax(E_v2 - E_v^2, 0))     # posterior SD

    z_var  <- (res[["var"]] - gene_data[["sanity_var"]]) / SE_var
    p_var  <- 2 * pnorm(abs(z_var), lower.tail = FALSE)
    p_adj  <- p.adjust(p_var, method = "BH")

    expect_gt(
        min(p_adj),
        0.001,
        label = "min BH-adjusted p-value for variance z-test"
    )
})

# Test 2: Per-gene LTQ agreement via z-score ----------------------------------
# We compute per-cell z-scores z_gc = d_gc / se_gc (d_gc = our - ref),

test_that("Sanity() LTQ estimates agree with paper on Baron dataset", {
    ltq_ours   <- res[["mu"]] + res[["delta"]]   # G x C
    se_gc      <- sqrt(res[["var_delta"]] + res[["var_mu"]])  # G x C
    d_gc       <- ltq_ours - ref_ltq                          # G x C
    z_gc       <- as.vector(d_gc / se_gc)                     # flatten to vector

    # Each z_gc ~ N(0,1) under H0. Apply BH over all ~1M pairs:
    p_adj <- p.adjust(2 * pnorm(abs(z_gc), lower.tail = FALSE), method = "BH")

    expect_gt(
        min(p_adj),
        0.001,
        label = "min BH-adjusted p-value for per-cell LTQ z-test"
    )
})
