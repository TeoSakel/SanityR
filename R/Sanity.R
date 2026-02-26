#' Estimate gene-level expression using the Sanity model
#'
#' This function provides a user-friendly interface to the Sanity model for gene
#' expression analysis.
#'
#' @param x A numeric matrix of counts (rows = features, columns = cells).
#'
#'   Alternatively, a \linkS4class{SummarizedExperiment} or a
#'   \linkS4class{SingleCellExperiment} containing such counts.
#' @param vmin The minimum value for the gene-level variance (must be > 0).
#' @param vmax The maximum value for the gene-level variance.
#' @param nbin Number of variance bins to use.
#' @param a,b Gamma prior parameters (see Details).
#' @param assay.type A string specifying the assay of `x` containing the counts.
#' @param size.factors A numeric vector of cell-specific size factors.
#'   Values are internally scaled to have mean equal to 1, so only relative differences matter.
#'   Defaults to `sizeFactors(x)` if `x` is a `SingleCellExperiment`, or `colSums(x)` if `x` is a matrix.
#' @param name Name of assay to store the normalized values.
#' @param subset.row A vector specifying the subset of rows of `x` to process.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object specifying whether
#'   the calculations should be parallelized.
#' @param ... For the generic, further arguments to pass to each method.
#'
#'   For the `SummarizedExperiment` method, further arguments to pass to the
#'   `ANY` method.
#'
#'   For the `SingleCellExperiment` method, further arguments to pass to the
#'   `SummarizedExperiment` method.
#'
#' @return For `matrix`-like object it returns a named list with the following
#' elements (symbols as defined in the Supplementary Text of the publication):
#' \describe{
#'   \item{mu}{Posterior mean of log expression across cells \eqn{\mu_g}.}
#'   \item{var_mu}{Posterior variance of the mean expression
#'                 \eqn{\left(\delta \mu_g\right)^2}.}
#'   \item{var}{Posterior variance of expression across cells
#'              \eqn{\langle v_g \rangle}.}
#'   \item{delta}{Vector of log fold-changes for each cell relative to
#'                \eqn{\delta_{gc}}.}
#'   \item{var_delta}{Posterior variance of the cell-level fold-changes
#'                    \eqn{\epsilon_{gc}^2}.}
#'   \item{cell_size}{Vector of cell-specific size factors used in the model.}
#'   \item{likelihood}{Normalized likelihood across the evaluated variance grid
#'              \eqn{P\left(v_g \mid n_g \right)} for diagnostics.}
#' }
#'
#' If called on a \linkS4class{SingleCellExperiment} or
#' \linkS4class{SummarizedExperiment} it appends the following columns to the
#' `rowData` slot:
#' \describe{
#'   \item{sanity_log_activity_mean}{`mu`}
#'   \item{sanity_log_activity_mean_sd}{`sqrt(var_mu)`}
#'   \item{sanity_activity_sd}{`sqrt(var)`}
#' }
#' and appends the following assays (assuming `name="logcounts"`):
#' \describe{
#'   \item{assay(x, "logcounts")}{`mu + delta`}
#'   \item{assay(x, "logcounts_sd")}{`sqrt(var_mu + var_delta)`}
#' }
#'
#' @details
#'
#' The method models gene activity using a Bayesian framework, assuming a Gamma
#' prior on expression and integrating over cell-level variability. It returns
#' posterior estimates for mean expression (`mu`), cell-specific deviations
#' (`delta`), and their variances, as well as expression variance (`var`).
#' *Expected* log-normalized counts are computed by combining mean expression
#' and cell-specific log-fold changes. The *standard deviation* of log-counts is
#' computed by summing the variances of the components.
#'
#' ## `size.factors` behavior
#'
#' The `size.factors` provided are used to compute the `cell_size` by first
#' normalizing them to have a mean of 1 and then scaling them by the average
#' library size (mean of column sums of `x`). The `cell_size` parameter quantifies
#' the expected library size (total number of UMIs) in each cell. The default
#' value for matrix method is consistent with the original publication and
#' treats any variation in library size as real so counts are partially
#' explained by them. However, in the case of `x` being a `SingleCellExperiment`,
#' we defer to `sizeFactors(x)`, which allows users to apply their own
#' normalization strategy consistenly across their analysis. In this case if
#' `sizeFactors(x)` is missing (`NULL`), we will assume all cells have the same
#' library size treating any variation as noise which will increase the
#' uncertainty in the estimates of gene expression (`delta`).
#'
#' ## Gamma Prior:
#'
#' The model adopts a Bayesian framework by placing a Gamma prior `Gamma(a, b)`
#' over the gene activity, where `a` is the shape and `b` the rate parameter,
#' respectively. This allows for flexible regularization and uncertainty
#' modeling. The posterior likelihood is estimated by integrating over possible
#' values of the variance in expression.
#'
#' Intuitively:
#' \itemize{
#'   \item `a` acts as a pseudo-count added to the total count of the gene.
#'   \item `b` acts as a pseudo-count penalizing deviations from the average.
#'   expression — i.e., it regularizes the total number of UMIs that differ from
#'   the expected value.
#' }
#'
#' Setting `a=1` and `b=0` corresponds to an uninformative (uniform) prior,
#' which was used in the original Sanity model publication.
#'
#' @references Breda, J., Zavolan, M., & van Nimwegen, E. (2021). Bayesian
#' inference of gene expression states from single-cell RNA-seq data.
#' *Nature Biotechnology*, 39, 1008–1016.
#' \url{https://doi.org/10.1038/s41587-021-00875-x}
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' sce <- simulate_independent_cells(N_cell=500, N_gene=100)
#'
#' # Standard Sanity normalization
#' sce_norm <- Sanity(sce)
#' logcounts(sce_norm)[1:5,1:5]
#'
#' # Using size factors
#' sf <- colSums(counts(sce))
#' sizeFactors(sce) <- sf / mean(sf)
#' sce_norm2 <- Sanity(sce)
#' logcounts(sce_norm2)[1:5,1:5]
#'
#' @name Sanity
NULL

#' @export
#' @rdname Sanity
setGeneric("Sanity", function(x, ...) standardGeneric("Sanity"))

#' @export
#' @rdname Sanity
#' @importFrom BiocParallel bpparam bplapply
setMethod("Sanity", "ANY", function(x, size.factors=colSums(x),
                                    vmin=0.001, vmax=50, nbin=160L,
                                    a=1, b=0, BPPARAM=bpparam()) {
    stopifnot("Minimum variance must be positive"=vmin > 0)
    stopifnot("Maximum variance must be greater than the minimum"=vmax > vmin)
    stopifnot("Both Gamma parameters must be positive"=a >= 0 & b >= 0)

    # Extract matrix stats
    C <- ncol(x)
    G <- nrow(x)

    # compute `cell_size` from `size.factors` and the average library size
    if (is.null(size.factors)) {
        # If no size factors are provided, assume all cells have the same size
        # This branch is probably called when x is `SingleCellExperiment` and `sizeFactors(x)` is NULL`
        size.factors <- rep(1, C)
    } else {
        stopifnot("Length of size factors must match number of cells"=length(size.factors) == C)
        stopifnot("All size.factors must be positive" = all(size.factors > 0))
        size.factors <- size.factors / mean(size.factors)
    }
    cell_size <- mean(colSums(x)) * size.factors

    # Parallel processing over genes
    results <- bplapply(
        seq_len(G),
        function(g) get_gene_expression_level(
            counts=x[g, ],
            cell_size=cell_size,
            vmin=vmin,
            vmax=vmax,
            numbin=nbin,
            a=a, b=b
        ),
        BPPARAM=BPPARAM
    )

    list(
        # Gene-metrics: Vector(G)
        mu=vapply(results, "[[", 0, "mu"),
        var_mu=vapply(results, "[[", 0, "var_mu"),
        var=vapply(results, "[[", 0, "var"),
        # Cell-specific: Matrix(G, C)
        delta=t(vapply(results, "[[", rep(0, C), "delta")),
        var_delta=t(vapply(results, "[[", rep(0, C), "var_delta")),
        cell_size=cell_size,  # Vector(C)
        # Variance Likelihood: Matrix(G, nbin)
        likelihood=t(vapply(results, "[[", rep(0, nbin), "lik"))
    )
})

#' @export
#' @rdname Sanity
#' @importFrom methods callNextMethod
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment rowData rowData<- assay assay<-
#' @importFrom scuttle .subset2index
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
setMethod(
    "Sanity",
    "SummarizedExperiment",
    function(x, ..., assay.type="counts", name="logcounts", subset.row=NULL) {
    if (!is.null(subset.row)) {
        subset.row <- .subset2index(subset.row, x, byrow=TRUE)
        x <- x[subset.row, ]
    }

    res <- callNextMethod(x=assay(x, assay.type), ...)

    # Check for genes with bad fit
    # check also: https://github.com/jmbreda/Sanity/issues/19#issuecomment-1217563095
    lik <- res[["likelihood"]]
    res[["entropy"]] <- -rowSums(lik * log2(lik), na.rm=TRUE) / log2(ncol(lik))
    problematic <- sum(res[["entropy"]] > .9)
    if (problematic > 0) {
        msg_fmt <- paste0(
            "There are %d genes whose posterior distribution of gene variance ",
            "has high entropy (> 0.9).\nConsider using different range or ",
            "dropping them from downstream analysis."
        )
        msg <- sprintf(msg_fmt, problematic)
        warning(msg)
    }

    # Append gene-level metrics to rowData
    rd <- DataFrame(
        sanity_log_activity_mean=res[["mu"]],
        sanity_log_activity_mean_sd=sqrt(res[["var_mu"]]),
        sanity_activity_sd=sqrt(res[["var"]]),
        sanity_entropy=res[["entropy"]]
    )
    rowData(x) <- cbind(rowData(x), rd)

    # Add cell-level metrics as new assays
    assay(x, name, withDimnames=FALSE) <- with(res, mu + delta)
    mu_sd <- with(res, sqrt(var_mu + var_delta))
    assay(x, paste0(name, "_sd"), withDimnames=FALSE) <- mu_sd
    return(x)
    }
)

#' @export
#' @rdname Sanity
#' @importFrom BiocGenerics sizeFactors
#' @importFrom methods callNextMethod
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setMethod(
    "Sanity",
    "SingleCellExperiment",
    function(x, size.factors=sizeFactors(x), ...) {
        callNextMethod(x=x, size.factors=size.factors, ...)
    }
)
