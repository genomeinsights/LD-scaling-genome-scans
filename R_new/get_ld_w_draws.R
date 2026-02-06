#' Precompute LD-weighted windows (ld_w) across random LD-decay quantiles
#'
#' For each chromosome in a GDS object, this function generates multiple draws
#' of LD window sizes based on random quantile values \eqn{\rho_w} of the
#' LD-decay curve, and then computes a local LD measure \eqn{ld_w} for each SNP.
#'
#' For each draw:
#' \enumerate{
#'   \item A quantile \eqn{\rho_w \sim \mathrm{Uniform}(\rho_{\min}, \rho_{\max})} is sampled.
#'   \item The corresponding physical distance threshold \eqn{d_{\rho_w}} is computed from
#'         the LD-decay rate \code{a} (using \code{\link{d_from_rho}}).
#'   \item For each SNP, \eqn{ld_w} is computed as the median \eqn{r^2} between that SNP and
#'         all other SNPs on the same chromosome within distance \eqn{d_{\rho_w}}.
#' }
#'
#' @param gds An open GDS object.
#' @param decay_tbl Output from \code{\link{ld_decay}} (must contain \code{$summary} with per-chromosome \code{a}).
#' @param slide_win_ld Integer, sliding window size (in SNPs) passed to \code{\link[SNPRelate]{snpgdsLDMat}} via \code{\link{get_el}}.
#' @param n_draws Integer, number of LD-window draws (distinct \code{rho_w}) to generate.
#' @param rho_min Numeric in (0,1). Lower bound of Uniform(\code{rho_min}, \code{rho_max}) for \code{rho_w}.
#' @param rho_max Numeric in (0,1). Upper bound of Uniform(\code{rho_min}, \code{rho_max}) for \code{rho_w}.
#' @param n_cores Integer, number of cores used for LD calculations (passed to \code{get_el}).
#'
#' @return A \code{data.table} with one row per chromosome × draw, containing \code{Chr}, \code{draw},
#'   \code{rho_w}, \code{d_th}, and a list-column \code{ld_w}. Each \code{ld_w} element is a
#'   \code{data.table} with columns \code{SNP} and \code{median_r2} (one row per SNP on that chromosome).
#'
#' @export
get_ld_w_draws <- function(gds,
                           decay_tbl,
                           slide_win_ld = 1000,
                           n_draws = 100,
                           rho_min = 0.9,
                           rho_max = 1.0) {
  
  t1 <- Sys.time()
  
  if(rho_max==1) rho_max <- rho_max-1e-9
  
  n_draws <- as.integer(n_draws)
  if (!is.finite(n_draws) || n_draws <= 0L) stop("`n_draws` must be a positive integer.")
  if (!is.finite(rho_min) || !is.finite(rho_max) || rho_min <= 0 || rho_max >= 1 || rho_min >= rho_max) {
    stop("`rho_min` and `rho_max` must satisfy 0 < rho_min < rho_max < 1.")
  }
  if (is.null(decay_tbl) || is.null(decay_tbl$summary)) {
    stop("`decay_tbl` must be supplied and contain a `$summary` table with per-chromosome decay parameter `a`.")
  }
  
  ids  <- read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)
  
  ## Sample rho_w values once and reuse across chromosomes
  rho_ws <- stats::runif(n_draws, rho_min, rho_max)
  
  out_by_chr <- vector("list", length(chrs))
  names(out_by_chr) <- as.character(chrs)
  
  
  for (ch in chrs) {
    cat(ch, "..")
    
    chr_idx <- which(ids$snp_chr == ch)
    snp_ids_chr <- ids$snp_id[chr_idx]
    
    ## Edge list for chromosome (unidirectional pairs)
    ## Expected columns at least: SNP1, SNP2, r2, d
    el <- get_el(gds, idx = chr_idx, slide_win_ld = slide_win_ld)
    
    ## Decay parameter for this chromosome
    summ <- decay_tbl$summary[decay_tbl$summary$Chr == ch, ]
    if (nrow(summ) != 1L) {
      stop("decay_tbl$summary must contain exactly one row per chromosome (missing/duplicate for Chr = ", ch, ").")
    }
    a_chr <- summ$a
    
    out_by_chr[[as.character(ch)]] <- data.table::rbindlist(
      lapply(seq_len(n_draws), function(draw_i) {
        
        rho_w <- rho_ws[draw_i]
        d_th  <- as.integer(d_from_rho(a_chr, rho = rho_w))
        
        ## Filter edges to those within distance threshold
        sub <- el[d < d_th, .(SNP1, SNP2, r2)]
        
        ## True per-SNP local LD:
        ## treat each edge as contributing to BOTH endpoints
        med_dt <- data.table::rbindlist(list(
          sub[, .(SNP = SNP1, r2)],
          sub[, .(SNP = SNP2, r2)]
        ))[, .(median_r2 = stats::median(r2, na.rm = TRUE)), by = SNP]
        
        ## Ensure every SNP appears (even if it has no neighbors within d_th)
        res <- data.table::data.table(SNP = snp_ids_chr)
        res <- med_dt[res, on = "SNP"]  # left join, keeps SNP order
        
        data.table::data.table(
          Chr   = ch,
          draw  = draw_i,
          rho_w = rho_w,
          d_th  = d_th,
          ld_w  = list(res)
        )
      }),
      use.names = TRUE
    )
  }
  
  ld_w_draws <- data.table::rbindlist(out_by_chr, use.names = TRUE)
  
  t2 <- Sys.time()
  print(difftime(t2, t1))
  cat("\n\n")
  
  ld_w_draws
}

