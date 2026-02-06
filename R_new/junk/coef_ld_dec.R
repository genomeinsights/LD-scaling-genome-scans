#' Fit LD-decay coefficients within a chromosome window
#'
#' Fits a simple LD-decay model to binned LD–distance summaries using non-linear
#' least squares:
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#' where \eqn{b} is a fixed background LD level, \eqn{c} is the short-distance LD
#' intercept, and \eqn{a} controls the rate of decay with distance \eqn{d}.
#'
#' The function bins SNP-pair distances into windows of width \code{dist_unit} and
#' summarizes LD within each bin by the \code{q}-quantile of \code{r2}. Bins are
#' weighted by the number of SNP pairs contributing to the summary.
#'
#' @param el_ld A \code{data.table} (or object coercible to \code{data.table})
#'   with columns \code{pos1}, \code{pos2} (base-pair positions), and \code{r2}
#'   (pairwise LD; typically \eqn{r^2}) for SNP pairs.
#' @param q Numeric in \eqn{(0,1)}. Quantile of the \eqn{r^2} distribution within
#'   each distance bin used to summarize LD. Defaults to \code{0.95}.
#' @param dist_unit Positive numeric. Bin width for distance (in base pairs) when
#'   aggregating LD. Defaults to \code{5000}.
#' @param b Numeric in \eqn{[0,1]}. Fixed background LD level used in the decay model.
#'   Defaults to \code{0.05}.
#'
#' @return A named numeric vector with elements \code{c} and \code{a}. If the fit fails
#'   or insufficient data are available, returns \code{c(c = NA_real_, a = NA_real_)}.
#'
#' @examples
#' \dontrun{
#' coefs <- coef_ld_dec(el_ld, q = 0.95, dist_unit = 5000, b = 0.05)
#' }
#'
#' @export
coef_ld_dec <- function(el_ld, q = 0.95, dist_unit = 5000, b = 0.05) {
  
  # Defensive checks (lightweight, package-friendly)
  if (!is.finite(q) || q <= 0 || q >= 1) stop("`q` must be in (0, 1).")
  if (!is.finite(dist_unit) || dist_unit <= 0) stop("`dist_unit` must be > 0.")
  if (!is.finite(b) || b < 0 || b > 1) stop("`b` must be in [0, 1].")
  
  # Bin summaries: q-quantile of r2 per bin + bin weights (#pairs) + mean distance
  bin_dt <- el_ld[
    ,
    .(
      r2_q      = stats::quantile(r2, probs = q, na.rm = TRUE, names = FALSE),
      w         = .N,
      dist_mean = mean(dist, na.rm = TRUE)
    ),
    by = dist_bin
  ][order(dist_bin)]
  
  # Need enough bins to fit a 2-parameter curve robustly
  if (nrow(bin_dt) < 5L) return(c(c = NA_real_, a = NA_real_))
  
  # Starting values:
  # c ~ max observed bin quantile, a ~ inverse of typical distance scale
  d0 <- stats::median(bin_dt$dist_mean[bin_dt$dist_mean > 0], na.rm = TRUE)
  if (!is.finite(d0)) return(c(c = NA_real_, a = NA_real_))
  
  starts <- list(
    c = max(bin_dt$r2_q, na.rm = TRUE),
    a = 1 / (d0 + 1e-6)
  )
  
  # Fit. Use 'port' for bound constraints.
  fit <- tryCatch(
    stats::nls(
      r2_q ~ b + (c - b) / (1 + a * dist_mean),
      data = bin_dt,
      start = starts,
      weights = w,
      algorithm = "port",
      lower = c(c = 0, a = 0),
      upper = c(c = 1, a = Inf)
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(c(c = NA_real_, a = NA_real_))
  
  co <- stats::coef(fit)
  
  # Ensure stable, named return (even if nls returns in different order)
  out <- c(c = unname(co["c"]), a = unname(co["a"]))
  out
}
