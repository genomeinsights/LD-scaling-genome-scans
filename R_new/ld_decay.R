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
#'   if \code{geno} is used. Must contain \code{marker}, \code{Chr}, \code{Pos}.
#' @param q Numeric in (0,1). Quantile used both for background LD and LD-decay summarization.
#' @param n_sub Integer. Number of SNPs (or SNP pairs, depending on implementation of \code{get_bg_ld})
#'   to use when estimating background LD.
#' @param slide_win_ld Integer. Sliding window size passed to \code{\link[SNPRelate]{snpgdsLDMat}}.
#' @param window_size Numeric. Genomic window size (bp) for LD-decay estimation.
#' @param step_size Numeric. Genomic step size (bp) between windows.
#' @param dist_unit Numeric. Distance bin width (bp) for LD-decay fitting within windows.
#' @param n_cores_ld Integer. Number of cores used for LD calculations.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{data}: LD-decay parameters per window (see \code{\link{ld_decay_by_chr_win}}).
#'   \item \code{summary}: Per-chromosome summary with columns \code{Chr_size}, \code{c}, \code{a}, \code{b}.
#'   \item \code{b}: Background LD level (upper quantile of inter-chromosomal \eqn{r^2}).
#' }
#'
#' @examples
#' \dontrun{
#' LD_decay <- ld_decay(gds, q = 0.95, n_sub = 5000)
#' LD_decay$summary
#' }
#'
#' @export
#' 
#' 


ld_decay <- function(gds,
                     geno = NULL,
                     map = NULL,
                     q = 0.95,
                     n_sub = 5000,
                     slide_win_ld = 1000,
                     window_size = 1e6,
                     step_size = 5e5,
                     dist_unit = 5000,
                     n_cores = 8) {
  
  t1 <- Sys.time()
  
  if (!is.finite(q) || q <= 0 || q >= 1) stop("`q` must be in (0,1).")
  if (!is.finite(window_size) || window_size <= 0) stop("`window_size` must be > 0.")
  if (!is.finite(step_size) || step_size <= 0) stop("`step_size` must be > 0.")
  if (!is.finite(dist_unit) || dist_unit <= 0) stop("`dist_unit` must be > 0.")
  
  ## Create temporary GDS if needed
  if (is.null(gds)) {
    if (is.null(geno) || is.null(map)) stop("If `gds` is NULL, both `geno` and `map` must be provided.")
    gds_path <- tempfile(fileext = ".gds")
    message("Generating temporary GDS object")
    
    gds <- create_gds_from_geno(geno, map, gds_path = gds_path)
    
    ## ensure we close and delete temporary file
    on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) }, add = TRUE)
  }
  
  ## Background LD (upper quantile of inter-chromosomal r^2)
  b <- get_bg_ld(gds, n_sub = n_sub, q = q)
  message(sprintf("q%.0f of between-chromosome LD (background LD) is b = %.4f", 100*q, b))
  
  ## LD-decay in sliding windows along chromosomes
  decay_data <- ld_decay_by_chr_win(
    gds,
    slide_win_ld = slide_win_ld,
    q = q,
    dist_unit = dist_unit,
    window_size = window_size,
    step_size = step_size,
    b = b,
    n_cores_ld = n_cores_ld
  )
  
  decay_summary <- decay_data[
    ,
    .(
      Chr_size = max(end, na.rm = TRUE),
      c = median(c, na.rm = TRUE),
      a = median(a, na.rm = TRUE),
      b = b
    ),
    by = Chr
  ]
  
  t2 <- Sys.time()
  print(difftime(t2, t1))
  cat("\n")
  
  list(data = decay_data, summary = decay_summary, b = b)
}


#' Create a SNPRelate GDS from a genotype matrix and map
#' @keywords internal
create_gds_from_geno <- function(geno, map, gds_path) {
  
  map <- data.table::as.data.table(map)
  
  req <- c("marker", "Chr", "Pos")
  if (!all(req %in% names(map))) stop("`map` must contain columns: marker, Chr, Pos.")
  stopifnot(ncol(geno) == nrow(map))
  
  SNPRelate::snpgdsCreateGeno(
    gds.fn          = gds_path,
    genmat         = t(round(geno)), # SNPRelate expects SNP x sample when snpfirstdim=TRUE
    sample.id      = paste0("ind_", seq_len(nrow(geno))),
    snp.id         = map$marker,
    snp.chromosome = map$Chr,
    snp.position   = map$Pos,
    snpfirstdim    = TRUE
  )
  
  SNPRelate::snpgdsOpen(gds_path)
}


#' Read SNP ids/positions/chromosomes from a GDS
#' @keywords internal
read_gds_ids <- function(gds) {
  list(
    snp_id  = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.id")),
    snp_chr = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.chromosome")),
    snp_pos = gdsfmt::read.gdsn(gdsfmt::index.gdsn(gds, "snp.position"))
  )
}


#' Estimate LD-decay parameters in sliding windows along chromosomes
#'
#' @param gds An open GDS object.
#' @param slide_win_ld Integer, sliding window size for LD calculation passed to
#'   \code{\link[SNPRelate]{snpgdsLDMat}}.
#' @param q Numeric, quantile for summarizing LD decay within bins.
#' @param dist_unit Numeric, distance bin width for LD-decay estimation (bp).
#' @param window_size Numeric, genomic window size (bp) for fitting LD-decay within chromosomes.
#' @param step_size Numeric, step size (bp) between consecutive windows.
#' @param b Numeric, background LD level (lower asymptote) in the LD-decay model.
#' @param n_cores_ld Integer, number of cores to use for LD calculations.
#'
#' @return A \code{data.table} with one row per window and columns:
#' \itemize{
#'   \item \code{Chr}: chromosome identifier.
#'   \item \code{start}, \code{end}: genomic coordinates of the window.
#'   \item \code{c}, \code{a}: LD-decay parameters for that window.
#' }
#'
#' @export
ld_decay_by_chr_win <- function(gds,
                                slide_win_ld = 10000,
                                q = 0.95,
                                dist_unit = 5000,
                                window_size = 1e7,
                                step_size = 5e5,
                                b = 0.05,
                                n_cores = 1) {
  
  ids <- read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)
  ch <- "Chr1"
  out <- mclapply(chrs,function(ch){
    cat(ch, "..")
    
    idx <- which(ids$snp_chr == ch)
    
    ## get_el() is assumed to return a data.table with pos1, pos2, r2
    el <- get_el(gds, slide_win_ld = slide_win_ld, idx = idx, n_cores = 1)
    
    el <- el[abs(pos2 - pos1) <= window_size]
    
    if (nrow(el) == 0) {
      out[[as.character(ch)]] <- data.table::data.table(
        Chr = ch, start = numeric(0), end = numeric(0), c = numeric(0), a = numeric(0)
      )
      next
    }
    
    ## Define windows based on observed positions in el
    min_pos <- min(el$pos1, el$pos2, na.rm = TRUE)
    max_pos <- max(el$pos1, el$pos2, na.rm = TRUE)
    
    if (!is.finite(min_pos) || !is.finite(max_pos) || (max_pos - min_pos) < window_size) {
      out[[as.character(ch)]] <- data.table::data.table(
        Chr = ch, start = min_pos, end = min_pos + window_size, c = NA_real_, a = NA_real_
      )
      next
    }
    
    starts <- seq(min_pos, max_pos - window_size, by = step_size)
    ends <- starts + window_size
    
    el[, dist := abs(pos2 - pos1)]
    el[, dist_bin := dist %/% dist_unit]
    
    data.table::rbindlist(lapply(seq_along(starts), function(i) {
      
      sub <- el[pos1 >= starts[i] & pos1 < ends[i] & pos2 >= starts[i] & pos2 < ends[i]]
      
      coefs <- tryCatch(
        coef_ld_dec(sub, q = q, dist_unit = dist_unit, b = b),
        error = function(e) c(c = NA_real_, a = NA_real_)
      )
      
      data.table::data.table(Chr = ch, start = starts[i], end = ends[i], c = coefs["c"], a = coefs["a"])
    }))
    
  },mc.cores=n_cores)
 
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}


#' Fit LD-decay coefficients within a chromosome window
#'
#' Fits a simple LD-decay model of the form
#' \deqn{r^2(d) = b + \frac{c - b}{1 + a d}}
#' to binned LD–distance data, using non-linear least squares.
#'
#' @param el_ld A \code{data.table} with columns \code{pos1}, \code{pos2}, and \code{r2}
#'   (pairwise LD) for SNP pairs.
#' @param q Numeric in (0,1). Quantile of the \eqn{r^2} distribution within distance bins.
#' @param dist_unit Positive numeric. Distance bin width (bp).
#' @param b Numeric in [0,1]. Fixed background LD level (lower asymptote).
#'
#' @return A named numeric vector \code{c(c = ..., a = ...)}. Returns \code{NA} values on failure.
#'
#' @export

coef_ld_dec <- function(el_ld, q = 0.95, dist_unit = 5000, b = 0.05) {
  
  if (!is.finite(q) || q <= 0 || q >= 1) stop("`q` must be in (0,1).")
  if (!is.finite(dist_unit) || dist_unit <= 0) stop("`dist_unit` must be > 0.")
  if (!is.finite(b) || b < 0 || b > 1) stop("`b` must be in [0,1].")
  
  
  req <- c("pos1", "pos2", "r2")
  if (!all(req %in% names(el_ld))) stop("`el_ld` must contain columns: pos1, pos2, r2.")
  
  bin_dt <- el_ld[
    ,
    .(
      r2_q      = stats::quantile(r2, probs = q, na.rm = TRUE, names = FALSE),
      w         = .N,
      dist_mean = mean(dist, na.rm = TRUE)
    ),
    by = dist_bin
  ][order(dist_bin)]
  
  if (length(unique(sub$dist_bin)) < 5) return(NA)
  
  d0 <- stats::median(bin_dt$dist_mean[bin_dt$dist_mean > 0], na.rm = TRUE)
  if (!is.finite(d0)) return(c(c = NA_real_, a = NA_real_))
  
  starts <- list(
    c = max(bin_dt$r2_q, na.rm = TRUE),
    a = 1 / (d0 + 1e-6)
  )
  
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
  c(c = unname(co["c"]), a = unname(co["a"]))
}



