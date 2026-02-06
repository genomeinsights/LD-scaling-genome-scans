#' Compute SNP-level consistency scores (C) from OR draws
#'
#' Given a set of outlier region (OR) draws produced e.g. by \code{\link{get_ORs_for_draws}} across many parameter combinations and
#' methods, this function computes, for each SNP and method, a consistency score
#' \eqn{C} defined as the proportion of draws in which that SNP appears in at least
#' one outlier region.
#'
#' @param draws A \code{data.table} of OR draws, typically as returned by
#'   \code{get_draws()} or a similar function. It must contain:
#'   \itemize{
#'     \item \code{method}: factor or character vector indicating the outlier method.
#'     \item \code{ORs}: a list-column, where each element is a list or vector of SNP IDs
#'           belonging to the ORs discovered under that parameter draw and method.
#'   }
#' @param markers Character vector of SNP IDs giving the reference ordering for SNPs.
#'   C-scores will be returned in this order. SNPs that never appear in any OR for
#'   a given method receive a score of \code{0}.
#'
#' @details
#' For each method, the function:
#' \enumerate{
#'   \item Flattens all ORs across parameter draws into a single vector of SNP IDs.
#'   \item Counts how often each SNP appears in at least one OR.
#'   \item Divides these counts by the total number of draws for that method, yielding
#'         a consistency score \eqn{C \in [0,1]}.
#' }
#'
#' SNPs that never occur in any OR for a given method are assigned \eqn{C = 0} for that
#' method. The output is then aligned to the order of \code{markers}.
#'
#' @return A \code{data.table} with one row per SNP (in the order of \code{markers}) and
#' one column per method. Column names are prefixed by \code{"C_"}, for example
#' \code{"C_EMX"}, \code{"C_LFMM"}, etc.
#'
#' @examples
#' \dontrun{
#' # Suppose 'draws' is a data.table with columns 'method' and 'ORs'
#' markers  <- ids$snp_id   # all SNP IDs in the dataset
#' C_scores <- get_C(draws, markers = markers)
#'
#' head(C_scores)
#' }
#'
#' @export

get_C <- function(draws,markers){
  n_total_draws <- draws[,.N,by=method][1,N]
  C_scores <- draws[,data.table(table(unlist(ORs))),by=method]
  setnames(C_scores,"V1","marker")
  C_scores[,C:=N/n_total_draws]
  
  
  meth <- "emx_F_prime" 
  C_scores <- do.call(cbind,lapply(unique(draws$method),function(meth){
    
    out <- C_scores[method==meth][match(markers,marker),C]
    out[is.na(out)] <- 0
    out
    
  }))
  C_scores <- data.table(C_scores)
  setnames(C_scores,colnames(C_scores),paste0("C_",unique(draws$method)))
  return(C_scores)
}
#C_scores[,plot(C_emx_F_prime)]
