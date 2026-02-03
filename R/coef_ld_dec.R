
#' Fit LD-decay coefficients within a chromosome window
#'
#' Fits a simple LD-decay model of the form
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#' to binned LD–distance data, using non-linear least squares.
#'
#' @param el_ld A \code{data.table} with columns \code{pos1}, \code{pos2}, and \code{r2}
#'   (pairwise LD) for SNP pairs.
#' @param q Numeric, quantile of the \eqn{r^2} distribution within distance bins used
#'   when summarizing LD. Defaults to \code{0.95}.
#' @param dist_unit Numeric, bin width for distance (in base pairs) when aggregating LD.
#' @param b Numeric, fixed background LD level used in the decay model.
#'
#' @return A named numeric vector with elements \code{c} and \code{a}, the fitted
#'   LD-decay parameters. If the fit fails, returns \code{NA} values.
#'
#' @examples
#' \dontrun{
#' coefs <- coef_ld_dec(el_ld, q = 0.95, dist_unit = 5000, b = 0.05)
#' }
#'
#' @export

coef_ld_dec <- function(el_ld, q = 0.95, dist_unit = 5000, b = 0.05) {
  
  el_ld[, dist := abs(pos2 - pos1)]
  el_ld[, dist_bin := floor(dist/dist_unit)]
  bin_dt <- el_ld[, .(r2_q = quantile(r2, na.rm = TRUE,prob=q),
                      w = .N,
                      dist_mean = mean(dist)),
                  by = dist_bin][order(dist_bin)]
  
  if (nrow(bin_dt) < 5) return(c(c = NA, a = NA))
  
  starts<- c(c=bin_dt[,max(r2_q, na.rm = TRUE)],a=bin_dt[,1 / (median(dist_mean[dist_mean > 0], na.rm = TRUE) + 1e-6)])
  
  fit <- tryCatch({
    fit <- nls(r2_q ~ b + (c-b) / (1 + a * dist_mean),
               data = bin_dt,
               start = starts,
               w = w,
               algorithm = "port", ## this is necessary to get convergence
               lower = c(0,  0),
               upper = c(1, Inf))
    
  }, error = function(e) NULL)
  coef(fit)
}
"grey"