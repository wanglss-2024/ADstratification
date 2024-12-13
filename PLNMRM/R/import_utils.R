## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##  INTERNAL FUNCTIONS ---------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Internal function to find the most comprehensive set of common samples between a count table and a covariates data.frame
common_samples <- function(counts, covariates) {
  ## Have default samples names been created when matching samples?
  default_names <- FALSE
  ## Sanity checks:
  name_warning <- paste("There are no matching names in the count matrix and the covariates data.frame.",
                        "Function will proceed assuming:",
                        "- samples are in the same order;",
                        "- samples are rows of the abundance matrix.", sep = "\n")
  ## no sample names in covariates: create sample names
  if (is.null(rownames(covariates))) {
    warning(name_warning)
    if (nrow(counts) != nrow(covariates)) {
      stop("Incompatible dimensions")
    }
    if (is.null(rownames(counts))) rownames(counts) <- paste0("Sample_", 1:nrow(counts))
    rownames(covariates) <- rownames(counts)
    default_names <- TRUE
  }
  ## Attempt name matching between covariates and count data
  count_names <- unlist(dimnames(counts))
  if (is.null(count_names) || !any(count_names %in% rownames(covariates))) {
    # If name matching is impossible, abort if
    ## - dimension are incompatible or
    ## - row (samples) names are conflicting
    warning(name_warning)
    if (nrow(counts) != nrow(covariates)) {
      stop("Incompatible dimensions")
    }
    if (!is.null(rownames(counts))) {
      stop("Conflicting samples names in count matrix and covariates data frames")
    }
    rownames(counts) <- rownames(covariates)
    default_names <- TRUE
  }
  ## Check whether samples are stored as columns in the abundance matrix
  ## and transpose if that's the case
  ## Based on a heuristic of matching names
  sample_are_cols <- any(colnames(counts) %in% rownames(covariates))
  if (sample_are_cols) counts <- t(counts)
  ## Ensure consistency by using only common samples
  common_samples <- intersect(rownames(counts), rownames(covariates))
  if (length(common_samples) < nrow(counts)) {
    message(paste0(nrow(counts) - length(common_samples), " samples were dropped from the abundance matrix for lack of associated covariates."))
  }
  return(list(transpose_counts = sample_are_cols,
              common_samples   = common_samples,
              default_names    = default_names))
}

## scaling functions --------

## Sanitize offset to ensure consistency with count matrix
sanitize_offset <- function(counts, offset, ...) {
  p <- ncol(counts) ## number of features
  ## Sanity check: transform vector offset and column-matrices to full matrices
  if (is.vector(offset) || (is.matrix(offset) && ncol(offset) == 1)) {
    offset_samples <- if (is.vector(offset)) {
      names(offset)
    } else {
      rownames(offset)
    }
    offset <- matrix(rep(offset, p),
                     ncol = p,
                     dimnames = list(offset_samples, colnames(counts)))
  }
  ## Sanity check: rownames
  if (is.null(rownames(offset))) {
    stop("Rownames are used for sample matching.\nPlease specify them in the offset vector/matrix.")
  }
  ## Sanity checks: offsets are available for all samples
  if (anyNA(ii <- match(rownames(counts), rownames(offset)))) {
    stop(paste("Sample(s) "),
         paste(rownames(counts)[is.na(ii)], collapse = " and "),
         " from the count table lack an offset.\nConsider checking your offset (orientation, rownames).")
  }
  offset <- offset[rownames(counts), , drop = FALSE]
  ## Sanity checks: offset are available for all species
  if (ncol(offset) != p) {
    stop(paste("There should be one offset per feature in the count table.\nYou have",
               p,
               "features but",
               ncol(offset),
               "offsets."))
  }
  offset
}



## Numeric offset
offset_numeric <- function(counts, offset, ...) {
  sanitize_offset(counts, offset, ...)
}

## No offset
offset_none <- function(counts) {
  return(NULL)
}

## Total Sum Scaling offset
offset_tss <- function(counts) {
  rowSums(counts)
}

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## EXPORTED FUNCTIONS ---------------------
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @title Prepare data for use in PLN models
#' @name prepare_data
#'
#' @description Prepare data in proper format for use in PLN model and its variants. The function (i) merges a count table and
#' a covariate data frame in the most comprehensive way and (ii) computes offsets from the count table using one of several normalization schemes (TSS, CSS, RLE, GMPR, Wrench, etc). The function fails with informative messages when the heuristics used for sample matching fail.
#'
#' @param counts Required. An abundance count table, preferably with dimensions names and species as columns.
#' @param covariates Required. A covariates data frame, preferably with row names.
#' @param offset Optional. Normalization scheme used to compute scaling factors used as offset during PLN inference. Available schemes are "TSS" (Total Sum Scaling, default), "CSS" (Cumulative Sum Scaling, used in metagenomeSeq), "RLE" (Relative Log Expression, used in DESeq2), "GMPR" (Geometric Mean of Pairwise Ratio, introduced in Chen et al., 2018), Wrench (introduced in Kumar et al., 2018) or "none". Alternatively the user can supply its own vector or matrix of offsets (see note for specification of the user-supplied offsets).
#' @param ... Additional parameters passed on to [compute_offset()]
#'
#' @return A data.frame suited for use in [PLN()] and its variants with two specials components: an abundance count matrix (in component "Abundance") and an offset vector/matrix (in component "Offset", only if offset is not set to "none")
#' @note User supplied offsets should be either vectors/column-matrices or have the same number of column as the original count matrix and either (i) dimension names or (ii) the same dimensions as the count matrix. Samples are trimmed in exactly the same way to remove empty samples.

prepare_data <- function(counts, covariates, offset = "TSS", ...) {
  ## Convert counts and covariates to expected format
  counts     <- data.matrix(counts, rownames.force = TRUE)
  covariates <- as.data.frame(covariates)
  ## sanitize abundance matrix and covariates data.frame
  common <- common_samples(counts, covariates)
  samples <- common$common_samples
  ## sanitize offset
  if (is.numeric(offset) && is.vector(offset)) {
    offset <- matrix(offset, ncol = 1, dimnames = list(names(offset), NULL))
  }
  if (common$transpose_counts) counts <- t(counts)
  if (common$default_names) {
    rownames(counts) <- rownames(covariates) <- samples
    if (is.numeric(offset)) rownames(offset) <- samples
  }
  counts <- counts[samples, , drop = FALSE]
  ## Replace NA with 0s
  if (any(is.na(counts))) {
    counts[is.na(counts)] <- 0
    warning("NA values in count table replaced with 0.")
  }
  ## filter out empty samples
  empty_samples <- which(rowSums(counts) == 0)
  if (length(empty_samples)) {
    warning(paste0("Sample(s) ",
                   paste(samples[empty_samples], collapse = " and "),
                   " dropped for lack of positive counts."))
    samples <- samples[-empty_samples]
    counts <- counts[samples, ,drop = FALSE]
  }
  covariates <- covariates[samples, , drop = FALSE]
  if (is.null(names(covariates))) names(covariates) <- paste0("Variable", seq_along(covariates))
  ## compute offset
  offset     <- compute_offset(counts, offset, ...)
  ## prepare data for PLN
  result <- data.frame(Abundance = NA, ## placeholder for Abundance, to avoid using I() and inheriting "AsIs" class
                       covariates,
                       Offset    = NA ## placeholder for Offset, to avoid using I() and inheriting "AsIs" class
                       )
  result$Abundance <- counts
  result$Offset <- offset
  result
}


#' @title Compute offsets from a count data using one of several normalization schemes
#' @name compute_offset
#'
#' @description Computes offsets from the count table using one of several normalization schemes (TSS, CSS, RLE, GMPR, etc) described in the literature.
#'
#' @inheritParams prepare_data
#' @param ... Additional parameters passed on to specific methods (for now CSS and RLE)
#' @inherit prepare_data references
#' @return If `offset = "none"`, `NULL` else a vector of length `nrow(counts)` with one offset per sample.
#'
#' @importFrom stats mad median quantile
#' @export
#'
#' ## User supplied offsets
#' my_offset <- setNames(rep(1, nrow(counts)), rownames(counts))
#' compute_offset(counts, offset = my_offset)
compute_offset <- function(counts, offset = c("TSS", "none"), ...) {
  ## special behavior for data.frame
  if (inherits(offset, "data.frame")) {
    stop(
  "You supplied a data.frame to compute_offset(). Did you mean to supply a numeric matrix?
  Try converting your data.frame to a matrix with as.matrix()."
  )
  }
  ## special behavior for numeric offset
  if (is.numeric(offset)) {
    return(offset_numeric(counts, offset, ...))
  }
  ## Choose offset function
  offset <- match.arg(offset)
  offset_function <- switch(offset,
                            "TSS"    = offset_tss,
                            "none"   = offset_none
  )
  ## Ensure that counts is a matrix
  counts <- counts %>% data.matrix()
  ## Compute offset (with optional parameters)
  offset_function(counts, ...)
}
