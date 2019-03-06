#===============================================================================
# callimbalance.R
#===============================================================================

# Imports ====================================================================

#' @import allelicimbalance
#' @import betabernsum
#' @import chenimbalance
#' @import data.table
#' @import npbin




# Functions ====================================================================

#' @title Estimate Null Parameters
#'
#' @description estimate the shape parameters of a beta-binomial null dist
#'
#' @param data_frame a data frame containing allele counts
#' @param minimum_coverage filter variants to those with at least this coverage
#' @param n_breaks n_breaks for npbin
#' @param spline_order spline order for npbin
#' @param cores number of cores to use
#' @return list of null parameters
#' @export
estimate_null_parameters <- function(
  data_frame,
  minimum_coverage = 10,
  n_breaks = 11,
  spline_order = 4,
  cores = 1
) {
  dt.ct <- convert_to_data_table(
    npbin_preprocess_counts(data_frame),
    minimum_coverage = minimum_coverage
  )

  n <- nrow(dt.ct)
  breaks <- seq(0, 1, length.out = n_breaks)
  pi_init <- initialize_weights(dt.ct, n_breaks, spline_order)

  overall_model_estimate <- emBinBspl(
    dt.ct[, xm],
    dt.ct[, m],
    breaks = breaks,
    k = spline_order,
    pi.init = pi_init,
    ncores = cores,
    err.max = 1e-3,
    iter.max = 200
  )
  null_model_estimate <- estNull(
    dt.ct[, xm],
    dt.ct[, m],
    overall_model_estimate,
    init = NULL,
    iter.max = 200,
    ncores = cores,
    ub = rep(log(1e4), 2),
    err.max = 1e-4
  )
  null_model_estimate[["coef.null"]]
}

#' @title Compute Imbalance Statistics
#'
#' @description compute per-variant imbalance statistics, including p-values
#'   and effect sizes
#'
#' @param counts_frame data frame containing allele counts
#' @param null_shape1 first shape parameter of the null distribution
#'   (from npbin)
#' @param null_shape2 second shape parameter of the null distribution
#'   (from npbin)
#' @param lsse_shape1 first shape parameter of the LSSE estimate
#' @param lsse_shape2 second shape parameter of the LSSE estimate
#' @param minimum_coverage minimum coverage level for filtering variants
#' @return data frame, including allele counts and imbalance statistics
#' @export
compute_imbalance_statistics <- function(
  counts_frame,
  null_shape1,
  null_shape2,
  lsse_shape1,
  lsse_shape2,
  minimum_coverage = 10
) {
  counts_frame[["binom_pval"]] <- apply(
    counts_frame[c("coverage", "ref_count")],
    1,
    function(row) {
      if (is.na(row[[1]])) return(NA)
      if (row[[1]] < minimum_coverage) return(NA)
      binom.test(
        row[[2]],
        row[[1]],
        p = null_shape1 / (null_shape1 + null_shape2),
        alternative = "two.sided"
      )[["p.value"]]
    }
  )
  counts_frame[["binom_fdr"]] <- p.adjust(
    counts_frame[["binom_pval"]],
    method = "BH"
  )
  counts_frame[["beta_binom_pval"]] <- apply(
    counts_frame[c("coverage", "ref_count")],
    1,
    function(row) {
      if (is.na(row[[1]])) return(NA)
      if (row[[1]] < minimum_coverage) return(NA)
      beta_binom_pval(row[[2]], row[[1]], null_shape1, null_shape2)
    }
  )
  counts_frame[["beta_binom_fdr"]] <- p.adjust(
    counts_frame[["beta_binom_pval"]],
    method = "BH"
  )
  counts_frame[["log_posterior_allelic_fold_change"]] <- apply(
    counts_frame[c("coverage", "ref_count")],
    1,
    function(row) {
      if (is.na(row[[1]])) return(NA)
      if (row[[1]] < minimum_coverage) return(NA)
      log_posterior_allelic_fold_change(
        row[[2]],
        row[[1]] - row[[2]],
        lsse_shape1,
        lsse_shape2
      )[["lpafc"]]
    }
  )
  counts_frame[["lpafc_lower"]] <- apply(
    counts_frame[c("coverage", "ref_count")],
    1,
    function(row) {
      if (is.na(row[[1]])) return(NA)
      if (row[[1]] < minimum_coverage) return(NA)
      log_posterior_allelic_fold_change(
        row[[2]],
        row[[1]] - row[[2]],
        lsse_shape1,
        lsse_shape2
      )[["lower"]]
    }
  )
  counts_frame[["lpafc_upper"]] <- apply(
    counts_frame[c("coverage", "ref_count")],
    1,
    function(row) {
      if (is.na(row[[1]])) return(NA)
      if (row[[1]] < minimum_coverage) return(NA)
      log_posterior_allelic_fold_change(
        row[[2]],
        row[[1]] - row[[2]],
        lsse_shape1,
        lsse_shape2
      )[["upper"]]
    }
  )
  counts_frame
}


#' @title coverage frame
#'
#' @description data frame of coverage levels
#'
#' @param intermediate_frames list of intermediate data frames
#' @return data frame, the coverage levels
#' @export
coverage_frame <- function(intermediate_frames) {
  as.data.frame(
    lapply(
      intermediate_frames,
      function(data_frame) data_frame[["coverage"]]
    )
  )
}

#' @title ref frame
#'
#' @description data frame of reference allele counts
#'
#' @param intermediate_frames list of intermediate data frames
#' @return data frame, the ref counts
#' @export
ref_frame <- function(intermediate_frames) {
  as.data.frame(
    lapply(
      intermediate_frames,
      function(data_frame) data_frame[["ref_count"]]
    )
  )
}

#' @title compute total coverage
#'
#' @description compute total coverage level across samples
#'
#' @param intermediate_frames list of intermediate data frames
#' @return integer, the total coverage level
#' @export
compute_total_coverage <- function(intermediate_frames) {
  coverage <- coverage_frame(intermediate_frames)
  rowSums(replace(coverage, is.na(coverage), 0))
}

#' @title compute total ref count
#'
#' @description compute total reference allele count across samples
#'
#' @param intermediate_frames list of intermediate data frames
#' @return integer, the total reference allel count
#' @export
compute_total_ref_count <- function(intermediate_frames) {
  ref_count <- ref_frame(intermediate_frames)
  rowSums(replace(ref_count, is.na(ref_count), 0))
}

#' @title meta analyze imbalance statistics
#'
#' @description meta analyze imbalance statistics
#'
#' @param row a row of imbalance data
#' @return meta-analyzed effect size
#' @export
meta_analyze <- function(row) {
  if (is.na(row[[3]])) return(row[[4]])
  if (is.na(row[[4]])) return(row[[3]])
  (row[[1]] * row[[3]] + row[[2]] * row[[4]]) / (row[[1]] + row[[2]])
}
