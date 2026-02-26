# Helper Functions
write_gz_table <- function(x, path) {
    gz <- gzfile(path, open = "w")
    tryCatch(
        write.table(x, gz, col.names = TRUE, row.names = FALSE, quote = FALSE),
        finally = close(gz)
    )
}

# Download UMI counts
# md5:ff15753e2af3cf5271e902f373bc6f0a
# licence: CC-4.0
baron_url <- url("https://zenodo.org/records/4009187/files/Baron_UMI_counts.txt.gz")
baron_counts <- read.table(gzcon(baron_url, text = TRUE), header = TRUE, row.names = 1L)
baron_counts <- as.matrix(baron_counts)

# Cell data
# md5:789dea1233294bf81a0b542564e01af2
# licence: GPL-3 (https://github.com/jmbreda/Sanity)
baron_ct_url <- url("https://raw.githubusercontent.com/jmbreda/Sanity/refs/heads/master/reproducibility/data/Baron_Celltype.txt")
cell_data <- data.frame(
    cell_id  = colnames(baron_counts),
    celltype = scan(baron_ct_url, character(), quiet = TRUE),
    libsize  = colSums(baron_counts)
)
write_gz_table(cell_data, "inst/extdata/baron_cell_data.txt.gz")

# Gene data
# Sanity variance estimates md5: b8774cfac8866a0a04fe9eceb0ff4c87
# licence: CC-4.0
baron_var_url <- url("https://zenodo.org/records/4009187/files/Baron_Sanity_variance.txt.gz")
baron_var <- read.table(gzcon(baron_var_url, text = TRUE), header = TRUE)

gene_data <- data.frame(
    symbol = rownames(baron_counts),
    sum    = rowSums(baron_counts),
    sanity_var = baron_var[[1]]
)
write_gz_table(gene_data, "inst/extdata/baron_gene_data.txt.gz")

# Create subset of baron_counts to test agreement with paper results.
# Subsample genes using a 32x32 log-log grid on sum and sanity_var.
set.seed(42L)
n_bins   <- 32L
bins <- with(gene_data, {
    sum_brks <- seq(log(min(sum)),        log(max(sum)),        length.out = n_bins + 1L)
    var_brks <- seq(log(min(sanity_var)), log(max(sanity_var)), length.out = n_bins + 1L)
    sum_bin  <- cut(log(sum),        breaks = sum_brks, include.lowest = TRUE, labels = FALSE)
    var_bin  <- cut(log(sanity_var), breaks = var_brks, include.lowest = TRUE, labels = FALSE)
    split(seq_along(symbol), paste(sum_bin, var_bin, sep = "_"))
})

# One gene per bin (coverage pass)
selected <- vapply(bins, function(i) sample(i, 1L), 0L)
selected <- sort(selected)

# Data to normalize
baron_counts <- baron_counts[selected, ]
saveRDS(baron_counts, "inst/extdata/baron_UMI_counts.rds", compress = "xz")

# Download results from original paper
# Normalized Counts md5: bc2a148c8ea9b4896efce3a968330eaf
# licence: CC-4.0
baron_norm_url <- url("https://zenodo.org/records/4009187/files/Baron_Sanity_normalization.txt.gz")
baron_norm <- read.table(gzcon(baron_norm_url, text = TRUE), header = TRUE, row.names = 1L)
baron_norm <- as.matrix(baron_norm[selected, ])
saveRDS(baron_norm, "inst/extdata/baron_Sanity_normalized_counts.rds", compress = "xz")

