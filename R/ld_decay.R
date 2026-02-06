
#' Compute LD-decay summaries and background LD for a dataset
#'
#' High-level convenience function that:
#' \enumerate{
#'   \item Creates a GDS object if one is not supplied.
#'   \item Estimates background LD (\eqn{b}) as the upper quantile of inter-chromosomal \eqn{r^2}.
#'   \item Estimates LD-decay parameters (\code{c} and \code{a}) in sliding windows along chromosomes.
#' }
#'
#' @param gds Optional open GDS object. If \code{NULL}, \code{geno} and \code{map} must be
#'   provided, and a temporary GDS will be created internally.
#' @param geno Optional genotype matrix (\eqn{n_\mathrm{ind} \times n_\mathrm{SNP}}). Used
#'   only if \code{gds} is \code{NULL}.
#' @param map Optional data.frame or \code{data.table} with SNP annotations; must be provided
#'   if \code{geno} is used.
#' @param q Numeric, quantile used both for background LD and LD-decay summarization.
#' @param n_sub Integer, number of SNPs to use when estimating background LD.
#' @param slide_win_ld Integer, sliding window size passed to \code{snpgdsLDMat}.
#' @param window_size Numeric, genomic window size (bp) for LD-decay estimation.
#' @param step_size Numeric, genomic step size (bp) between windows.
#' @param n_cores_ld Integer, number of cores used for LD-decay estimation.
#'
#' @return A list with components:
#'   \itemize{
#'     \item \code{data}: LD-decay parameters per window (see \code{\link{ld_decay_by_chr_win}}).
#'     \item \code{summary}: Per-chromosome summary with columns \code{Chr_size}, \code{c}, and \code{a}.
#'     \item \code{b}: Background LD level (upper quantile of inter-chromosomal \eqn{r^2}).
#'   }
#'
#' @examples
#' \dontrun{
#' LD_decay <- ld_decay(gds, q = 0.95, n_sub = 5000)
#' LD_decay$summary
#' }
#'
#' @export

# 
ld_decay <- function(gds,
                     idx=NULL,
                     b,
                     q = 0.95,                               ## Quantile for LD-decay and background LD
                     n_sub = 5000,                           ## How many random SNP pairs to estimate background LD
                     slide_win_ld = 1000,                    ## Sliding window for LD-estimation 
                     window_size = 1e6,                      ## Window size to estimate LD-decay within chromosomes in bps
                     step_size = 5e5,                        ## Step size for window in bps
                     n_cores_ld =  8,                         ## Number of cores to estimate LD
                     dist_unit = 5000
) {
  t1 <- Sys.time()
  ## generate gds file
  
  cat("q95 of between chromosome LD (background LD) is b=",sprintf("%.4f", b),"\n")
  ids <- read_gds_ids(gds)
  if(is.null(idx)) idx <- seq_along(ids$snp_id)
  
  ## LD-decay data in windows within chromosomes
  
  decay_data <- ld_decay_by_chr_win(gds, idx, q = q, b = b, slide_win_ld = slide_win_ld,window_size = window_size,step_size = step_size,n_cores_ld =  n_cores_ld)
  
  decay_summary <- decay_data[,.(Chr_size=max(end),c=median(c,na.rm=TRUE),a=median(a,na.rm=TRUE),b=b),by=Chr]
  t2 <- Sys.time()
  print(difftime(t2,t1))
  cat("\n")
  return(list(data=decay_data,summary=decay_summary))
}


create_gds_from_geno <- function(geno, map, gds_path) {
  stopifnot(ncol(geno) == nrow(map))
  snpgdsCreateGeno(
    gds_path,
    genmat         = t(round(geno)),     # SNPRelate expects SNP × sample if snpfirstdim=TRUE
    sample.id      = paste0("ind_", seq_len(nrow(geno))), ## sample name not important
    snp.id         = map$marker,
    snp.chromosome = map$Chr,
    snp.position   = map$Pos,
    snpfirstdim    = TRUE
  )
  snpgdsOpen(gds_path)
}



read_gds_ids <- function(gds) {
  list(
    snp_id  = read.gdsn(index.gdsn(gds, "snp.id")),
    snp_chr = read.gdsn(index.gdsn(gds, "snp.chromosome")),
    snp_pos = read.gdsn(index.gdsn(gds, "snp.position"))
  )
}


#' Estimate LD-decay parameters in sliding windows along chromosomes
#'
#' Computes LD-decay parameters (\code{c} and \code{a}) in genomic windows along
#' each chromosome, based on pairwise \eqn{r^2} estimates from a GDS object.
#'
#' @param gds An open GDS object.
#' @param slide_win_ld Integer, sliding window size for LD calculation passed to
#'   \code{\link[SNPRelate]{snpgdsLDMat}}.
#' @param q Numeric, quantile for summarizing LD decay within bins (default \code{0.95}).
#' @param dist_unit Numeric, distance bin width for LD-decay estimation (bp).
#' @param window_size Numeric, genomic window size (bp) for fitting LD-decay within chromosomes.
#' @param step_size Numeric, step size (bp) between consecutive windows.
#' @param b Numeric, background LD level used as the lower asymptote in the LD-decay model.
#' @param n_cores_ld Integer, number of cores to use for LD calculations.
#'
#' @return A \code{data.table} with one row per window and columns:
#'   \itemize{
#'     \item \code{Chr}: chromosome identifier.
#'     \item \code{start}, \code{end}: genomic coordinates of the window.
#'     \item \code{c}, \code{a}: LD-decay parameters for that window (may be \code{NA} if fit fails).
#'   }
#'
#' @examples
#' \dontrun{
#' decay_data <- ld_decay_by_chr_win(gds, slide_win_ld = 1000,
#'                                   window_size = 1e6, step_size = 5e5)
#' }
#'
#' @export

ld_decay_by_chr_win <- function(gds, idx, slide_win_ld = 10000, q = 0.95, dist_unit = 5000, window_size=1e7,
                                step_size =5e+05,b = 0.05, n_cores_ld = 1) {
  ids <- read_gds_ids(gds)
  chrs <- unique(ids$snp_chr[idx])
  
  out <- list()
  #ch="Chr1"
  for(ch in chrs){
    cat(ch,"..")
    
    
    el <- get_el(gds,slide_win_ld = slide_win_ld,idx = which(ids$snp_chr[idx] == ch),n_cores = n_cores_ld)
    
    # get windows for ld-decay analyses
    min_pos <- min(el$pos1, el$pos2)
    max_pos <- max(el$pos1, el$pos2)
    starts <- seq(min_pos, max_pos - window_size, by = step_size)
    ends <- starts + window_size
    
    #i <- 1
    # el_ld <- sub
    out[[ch]] <- rbindlist(lapply(seq_along(starts),function(i) {
      
      sub <- el[pos1 >= starts[i] & pos1 < ends[i] & pos2 >= starts[i] & pos2 < ends[i]]
      
      coefs <- tryCatch({
        coef_ld_dec(sub,q = 0.95, dist_unit = 5000, b = b)
      }, error = function(e) NULL)
      
      
      if(is.null(coefs)){
        data.table(Chr = ch, start=starts[i],end=ends[i],c = NA, a = NA)  
      }else{
        data.table(Chr = ch, start=starts[i],end=ends[i],t(coefs))
      }
    }))
    
  }
  
  rbindlist(out, use.names = TRUE, fill = TRUE)
}


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

