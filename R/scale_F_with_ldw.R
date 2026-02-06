#' LD-scaled F statistic with permutation-based quantile transformation
#'
#' @export
scale_F_with_ldw <- function(F_vals, ld_w, n_rep = 10, n_inds, full = TRUE) {

  ## Coerce F_vals to matrix (n_snp x n_method)
  F_mat <- as.matrix(F_vals)
  n <- nrow(F_mat)

  if (length(ld_w) != n) stop("`ld_w` must have length equal to nrow(F_vals).")
  if (!is.finite(n_rep) || n_rep <= 0) stop("`n_rep` must be a positive integer.")
  n_rep <- as.integer(n_rep)

  if (!is.finite(n_inds) || n_inds < 3) stop("`n_inds` must be >= 3.")
  df2 <- as.integer(n_inds - 2)

  ld_w <- as.numeric(ld_w)
  ld_w[is.na(ld_w)] <- 0

  ## Helper: take n equally spaced order statistics from a sorted vector
  take_sorted_quantiles <- function(sorted_x, n_out) {
    m <- length(sorted_x)
    if (m == 0L) return(rep(NA_real_, n_out))
    ## positions 1..m mapped onto 1..n_out
    pos <- round(seq(1, m, length.out = n_out))
    sorted_x[pos]
  }

  ## Process each method/column independently
  res <- lapply(seq_len(ncol(F_mat)), function(j) {

    Fval <- F_mat[, j]
    F_mean <- mean(Fval, na.rm = TRUE)

    ## Observed LD-scaled statistic, rescaled to original mean
    F_prime_obs <- Fval * ld_w
    denom <- mean(F_prime_obs, na.rm = TRUE)

    ## If denom is 0 (e.g., ld_w all zeros), we cannot rescale; fall back gracefully.
    if (!is.finite(denom) || denom == 0) {
      F_prime_obs <- rep(F_mean, length(Fval))
      denom <- F_mean
    } else {
      F_prime_obs <- F_prime_obs / denom * F_mean
    }

    ## Rank bookkeeping
    ord <- order(F_prime_obs)
    F_prime_obs_sorted <- F_prime_obs[ord]

    ## Build permutation null by circularly shifting ld_w
    ## (preserves LD autocorrelation along SNP order)
    null_all <- vector("numeric", n * n_rep)

    for (k in seq_len(n_rep)) {
      start <- sample.int(n, 1)
      shifted <- c(ld_w[start:n], ld_w[1:(start - 1)])

      null_k <- Fval * shifted
      denom_k <- mean(null_k, na.rm = TRUE)

      if (!is.finite(denom_k) || denom_k == 0) {
        null_k <- rep(F_mean, n)
      } else {
        null_k <- null_k / denom_k * F_mean
      }

      null_all[((k - 1) * n + 1):(k * n)] <- null_k
    }

    null_all <- null_all[is.finite(null_all)]
    null_perm_sorted <- sort(null_all)

    ## Empirical null quantiles sampled at n ranks
    null_perm <- take_sorted_quantiles(null_perm_sorted, n_out = n)

    ## Theoretical F null quantiles (same as your code)
    null_true <- stats::qf(rev(stats::ppoints(n)), df1 = 1, df2 = df2, lower.tail = FALSE)

    ## Quantile transformation + monotonicity correction
    q_F_prime_cor <- F_prime_obs_sorted - (null_perm - null_true)
    q_F_prime_cor_adj <- cummax(q_F_prime_cor)

    ## Back to original SNP order
    F_prime <- numeric(n)
    F_prime[ord] <- q_F_prime_cor_adj

    p_prime <- stats::pf(F_prime, df1 = 1, df2 = df2, lower.tail = FALSE)
    q_prime <- stats::p.adjust(p_prime, "fdr")

    if (!full) return(q_prime)

    qq_data <- data.table::data.table(
      null_true         = null_true,
      F_prime_obs       = F_prime_obs_sorted,
      F_obs_sorted      = sort(Fval),
      q_F_prime_cor     = q_F_prime_cor,
      q_F_prime_cor_adj = q_F_prime_cor_adj,
      null_perm         = null_perm
    )

    qq_data <- data.table::melt(
      qq_data,
      measure.vars = c("F_prime_obs", "F_obs_sorted", "q_F_prime_cor",
                       "null_perm", "q_F_prime_cor_adj")
    )

    list(
      qq_data  = qq_data,
      F_prime  = F_prime,
      p_prime  = p_prime,
      q_prime  = q_prime
    )
  })

  ## Return shape: if full=FALSE -> matrix of q-values; else list per method
  if (!full) {
    qmat <- do.call(cbind, res)
    colnames(qmat) <- colnames(F_mat)
    return(qmat)
  } else {
    names(res) <- colnames(F_mat)
    if (length(res) == 1L) return(res[[1]])
    return(res)
  }
}