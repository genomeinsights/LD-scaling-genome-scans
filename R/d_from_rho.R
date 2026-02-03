
#' Convert LD-decay parameter and quantile to a distance threshold
#'
#' For the LD-decay model
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d},}
#' this function returns the physical distance \eqn{d_\rho} at which the LD
#' curve has decayed to a given quantile position \eqn{\rho}.
#'
#' @param a Numeric scalar or vector of LD-decay rates.
#' @param rho Numeric scalar or vector of quantile positions, with values in \code{(0,1)}.
#'
#' @return A numeric vector of the same length as \code{a} and \code{rho} (after recycling),
#'   giving the distance threshold \eqn{d_\rho}.
#'
#' @examples
#' d_075 <- d_from_rho(a = 1e-6, rho = 0.75)
#'
#' @export

d_from_rho <- function(a, rho) {
  # d_ρ = (1/a)*(1/(1-ρ) - 1)
  (1 / a) * (1 / (1 - rho) - 1)
}