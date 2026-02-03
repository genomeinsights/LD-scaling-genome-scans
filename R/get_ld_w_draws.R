#' Precompute LD-weighted windows (ld_w) across random LD-decay quantiles
#'
#' For each chromosome in a GDS object, this function generates multiple draws
#' of LD window sizes based on random quantile values \eqn{\rho_w} of the
#' LD-decay curve, and then computes a local LD measure \eqn{ld_w} for each SNP.
#'
#' For each draw:
#' \enumerate{
#'   \item A quantile \eqn{\rho_w \sim \mathrm{Uniform}(\mathrm{min}, \mathrm{max})} is sampled.
#'   \item The corresponding physical distance threshold \eqn{d_{\rho_w}} is computed from
#'         the LD-decay rate \code{a} (using \code{\link{d_from_rho}}).
#'   \item For each SNP, \eqn{ld_w} is computed as the median \eqn{r^2} between that SNP and
#'         all other SNPs on the same chromosome within distance \eqn{d_{\rho_w}}.
#' }
#'
#' The result is a long-format table with one row per chromosome × draw combination,
#' and a list-column containing per-SNP \eqn{ld_w} values for that combination.
#'
#' @param gds An open GDS object containing genotype data (see
#'   \code{\link[SNPRelate]{snpgdsCreateGeno}} and \code{\link[SNPRelate]{snpgdsOpen}}).
#' @param slide_win_ld Integer, sliding window size (in SNPs) passed to
#'   \code{\link[SNPRelate]{snpgdsLDMat}} when computing LD. Defaults to \code{1000}.
#' @param n_draws Integer, number of LD-window draws (i.e. distinct \eqn{\rho_w} values)
#'   to generate. Defaults to \code{100}.
#' @param min Numeric, lower bound of the uniform distribution for \eqn{\rho_w}. Defaults
#'   to \code{0.9}.
#' @param max Numeric, upper bound of the uniform distribution for \eqn{\rho_w}. Defaults
#'   to \code{1.0}.
#' @param n_cores Integer, number of CPU cores to use for internal LD calculations. Currently
#'   passed to \code{\link{get_el}}. Defaults to \code{1}.
#'
#' @details
#' This function relies on:
#' \itemize{
#'   \item \code{\link{read_gds_ids}} to obtain SNP IDs, chromosomes, and positions.
#'   \item \code{\link{get_el}} to compute pairwise LD (\eqn{r^2}) and distances per chromosome.
#'   \item \code{\link{d_from_rho}} to convert LD-decay parameters and \eqn{\rho_w} to a
#'         physical distance threshold \eqn{d_{\rho_w}}.
#' }
#'
#' In the current implementation, the LD-decay rate \code{a} is taken from a global object
#' \code{LD_decay_3sp$summary}, which should contain per-chromosome estimates of the decay
#' parameter \code{a}. For use in a package, this would typically be refactored to pass the
#' decay summary explicitly as an argument.
#'
#' @return A \code{data.table} with one row per chromosome × draw combination, containing:
#' \itemize{
#'   \item \code{Chr}: chromosome identifier.
#'   \item \code{draw}: integer index of the draw (\code{1:n_draws}).
#'   \item \code{rho_w}: sampled quantile \eqn{\rho_w} used to define the LD window.
#'   \item \code{d_th}: physical distance threshold (in bp) corresponding to \eqn{\rho_w}.
#'   \item \code{ld_w}: list-column, each element a \code{data.table} with columns
#'         \code{SNP1} and \code{median_r2}, giving the local LD measure for each SNP
#'         on that chromosome.
#' }
#'
#' @examples
#' \dontrun{
#' ld_w_draws <- get_ld_w_draws(
#'   gds          = gds,
#'   slide_win_ld = 1000,
#'   n_draws      = 100,
#'   min          = 0.9,
#'   max          = 1.0,
#'   n_cores      = 4
#' )
#'
#' head(ld_w_draws)
#' }
#'
#' @export
#' 

get_ld_w_draws <- function(gds, decay_tbl,slide_win_ld=1000, n_draws = 100, min=0.9, max=1, n_cores=1){
  t1 <- Sys.time()
  
  ids <- read_gds_ids(gds)
  
  ld_w <- list()
  
  chrs <- unique(ids$snp_chr)
  rho_ws <- runif(n_draws,min,max)
  #ch <- chrs[1]
  for(ch in chrs){
    cat(ch,"..")
    idx <- which(ids$snp_chr==ch)
    
    el <- get_el(gds,idx,slide_win_ld=slide_win_ld, n_cores=n_cores)
    
    #rho_w <- rho_ws[1]
    ld_w[[ch]] <- rbindlist(lapply(rho_ws,function(rho_w){
      
      
      d_th <- as.integer(d_from_rho(decay_tbl$summary[Chr==ch,a],rho_w))
      
      out <- el[d<d_th,.(median_r2=median(r2,na.rm = TRUE)),by=SNP1]
      out <- out[match(ids$snp_id[idx],SNP1)]
      data.table(Chr=ch,
                 rho_w=rho_w,
                 d_th=d_th,
                 ld_w=list(out))
    }))
    #gc()
    ld_w[[ch]][,draw:=.I]
  }
  
  ld_w <- rbindlist(ld_w)
  
  
  t2 <- Sys.time()
  print(difftime(t2,t1))
  cat("\n\n")
  
  return(ld_w)
}
