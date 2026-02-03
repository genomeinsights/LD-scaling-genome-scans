
#' Estimate background LD between chromosomes
#'
#' Estimates the upper quantile of linkage disequilibrium (LD), measured as
#' \eqn{r^2}, for SNP pairs located on different chromosomes. This serves as a
#' background LD level (\eqn{b}) for downstream thresholding.
#'
#' @param gds An open GDS object.
#' @param n_sub Integer, number of SNPs to sample for estimating background LD.
#'   Defaults to \code{5000}.
#' @param q Numeric, quantile of the \eqn{r^2} distribution to use, typically
#'   \code{0.95}.
#'
#' @return A single numeric value giving the \code{q}-th quantile of \eqn{r^2}
#'   among inter-chromosomal SNP pairs.
#'
#' @examples
#' \dontrun{
#' b <- get_bg_ld(gds, n_sub = 5000, q = 0.95)
#' }
#'
#' @export

get_bg_ld <- function(gds, n_sub = 5000, q = 0.95) {
  ids <- read_gds_ids(gds)
  
  snp_pool <- sample(seq_along(ids$snp_id), min(n_sub, length(ids$snp_id)))
  #print(snp_pool)
  ld <- snpgdsLDMat(gds, snp.id = ids$snp_id[snp_pool], method = "r", slide = -1, verbose = FALSE)
  
  el <- as.data.table(reshape2::melt(ld$LD^2, value.name = "r2"))
  el[, Chr1 := ids$snp_chr[snp_pool][Var1]]
  el[, Chr2 := ids$snp_chr[snp_pool][Var2]]
  el[Chr1 != Chr2, quantile(r2, q, na.rm = TRUE)]
}