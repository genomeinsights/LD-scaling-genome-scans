
#' LD-scaled F statistic with permutation-based quantile transformation
#'
#' This function takes one or more columns of F-statistics (e.g. from LFMM or EMMAX),
#' scales them by a local LD measure (\code{ld_w}), and then applies a permutation-based
#' quantile transformation to obtain an LD-adjusted statistic \eqn{F'}.
#'
#' The input \code{F_vals} is coerced to a matrix; each column is processed independently.
#' Internally, \code{ld_w} is circularly permuted to generate a null distribution that
#' preserves the autocorrelation structure of LD along the genome.
#'
#' @param F_vals Numeric vector or matrix of F-statistics, with length \eqn{n_\mathrm{SNP}}
#'   (if a vector) or dimension \eqn{n_\mathrm{SNP} \times n_\mathrm{method}} (if a matrix).
#'   Each column is assumed to correspond to a different association method or model.
#' @param ld_w Numeric vector of length \eqn{n_\mathrm{SNP}}, giving the local LD measure
#'   (e.g. median \eqn{r^2} within a window) for each SNP.
#' @param n_rep Integer, number of LD permutations used to approximate the null distribution
#'   of the LD-scaled statistic. Defaults to \code{10}.
#' @param n_inds Integer, number of individuals in the association test. Used to set the
#'   denominator degrees of freedom \code{df2 = n_inds - 2} when converting \eqn{F'} to
#'   p-values under an F-distribution.
#' @param full Logical. If \code{TRUE} (default), returns both the transformed statistics
#'   and diagnostic Q–Q data for each method. If \code{FALSE}, returns only FDR-adjusted
#'   q-values for each method.
#'
#' @details
#' For each method/column in \code{F_vals}, the function:
#' \enumerate{
#'   \item Multiplies \code{F_vals} by \code{ld_w} and rescales the product to have the
#'         same mean as the original F-statistics.
#'   \item Generates \code{n_rep} circular permutations of \code{ld_w} and recomputes the
#'         LD-weighted statistic at each permutation, building an empirical null.
#'   \item Approximates the null quantiles by interpolation and aligns them with the
#'         corresponding theoretical F-distribution quantiles.
#'   \item Performs a quantile transformation by subtracting the deviation between the
#'         empirical and theoretical null distributions from the observed LD-scaled
#'         statistics.
#'   \item Applies a monotonicity correction via the cumulative maximum of the transformed
#'         quantiles to ensure that the final statistic is non-decreasing in rank.
#' }
#'
#' @return
#' If \code{full = TRUE}, a list of length equal to the number of columns in
#' \code{F_vals}, each element being a list with components:
#' \itemize{
#'   \item \code{F_prime}: Numeric vector of LD-scaled, quantile-adjusted F-statistics
#'         (\eqn{F'}).
#'   \item \code{p_prime}: Raw p-values computed from \eqn{F'} under \code{F(1, df2)}.
#'   \item \code{q_prime}: FDR-adjusted q-values (Benjamini–Hochberg).
#'   \item \code{qq_data}: A \code{data.table} containing empirical and theoretical
#'         quantiles for diagnostic Q–Q plots.
#' }
#' If \code{full = FALSE}, a matrix of q-values with one column per method.
#'
#' @examples
#' \dontrun{
#' # Example with a single method:
#' F_obs  <- rchisq(1000, df = 1)  # pseudo F-like statistic
#' ld_w   <- runif(1000, 0, 0.5)
#' res    <- scale_F_with_ldw(F_vals = F_obs, ld_w = ld_w, n_rep = 50, n_inds = 120)
#'
#' # Example with two methods:
#' F_mat  <- cbind(method1 = F_obs, method2 = F_obs * 1.2)
#' res2   <- scale_F_with_ldw(F_vals = F_mat, ld_w = ld_w, n_rep = 50, n_inds = 120)
#' }
#'
#' @export

scale_F_with_ldw <- function(F_vals, ld_w, n_rep = 10, n_inds,full=TRUE) {
  
  stopifnot(nrow(F_vals) == length(ld_w))
  ld_w[is.na(ld_w)] <- 0
  df2 = n_inds - 2
  
  out <- apply(as.matrix(F_vals),2,function(Fval){
    F_mean <- mean(Fval, na.rm = TRUE)
    F_prime_obs <- Fval * ld_w
    F_prime_obs <- F_prime_obs / mean(F_prime_obs, na.rm = TRUE) * F_mean
    n <- length(F_prime_obs)
    
    # Keep original order
    original_order <- order(F_prime_obs)
    F_prime_obs_sorted <- sort(F_prime_obs)
    
    # Step 3: permutation-based null distribution
    
    exp_perm <- unlist(replicate(n_rep, {
      # circular shift
      start <- sample(seq_len(n), 1)
      shifted <- c(ld_w[start:n], ld_w[1:(start-1)])
      
      null <- Fval * shifted[seq_along(Fval)]
      null / mean(null,na.rm = TRUE) * F_mean
      
    }, simplify = FALSE))
    
    
    null_perm <- approx(x=1:length(exp_perm), y=sort(exp_perm),n = n)$y
    
    # Step 5: theoretical F null quantiles
    null_true <- qf(rev(ppoints(n)), df1 = 1, df2 = df2, lower.tail = FALSE)
    
    # Step 6: quantile transformation
    q_F_prime_cor <- F_prime_obs_sorted - (null_perm - null_true)
    
    q_F_prime_cor_adj <- cummax(q_F_prime_cor)
    
    F_prime = q_F_prime_cor_adj[order(original_order)]
    p_prime = pf(F_prime,1,df2 = df2,lower.tail = FALSE)
    q_prime = p.adjust(p_prime,"fdr")
    
    # Return results
    if(full){
      
      qq_data <- data.table(
        null_true         = null_true,
        F_prime_obs       = F_prime_obs_sorted,
        F_obs_sorted      = sort(Fval),
        q_F_prime_cor     = q_F_prime_cor,
        q_F_prime_cor_adj = q_F_prime_cor_adj,
        null_perm         = null_perm
      )
      
      qq_data <- melt(
        qq_data,
        measure.vars = c("F_prime_obs", "F_obs_sorted", "q_F_prime_cor",
                         "null_perm", "q_F_prime_cor_adj")
      )
      
      return(list(
        qq_data = qq_data,
        F_prime=F_prime,
        p_prime=p_prime,
        q_prime=q_prime
      ))  
      
    }else{
      return(q_prime)  
    }
  })
  
  if(length(out)==1){
    out <- unlist(out,recursive = FALSE)  
  }
  return(out)
}