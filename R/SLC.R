#' Single-linkage clustering of SNPs based on LD and distance thresholds
#'
#' Partitions SNPs into clusters using single-linkage clustering based on LD and
#' distance thresholds. For each chromosome, SNPs connected by edges with \eqn{r^2}
#' above a threshold and distance below a threshold are grouped into connected
#' components (i.e., LD clusters / outlier regions). Singleton SNPs are returned
#' as clusters of size 1.
#'
#' @param gds Optional open GDS object. If \code{NULL}, \code{geno} and \code{map} must
#'   be provided, and a temporary GDS will be created.
#' @param geno Optional genotype matrix (\eqn{n_\mathrm{ind} \times n_\mathrm{SNP}}).
#' @param map Optional SNP map; required if \code{geno} is used.
#' @param idx Optional integer vector of SNP indices (into the GDS SNP arrays) to cluster.
#' @param decay_tbl Optional LD-decay output from \code{\link{ld_decay}} (must contain
#'   \code{$summary} with per-chromosome \code{a} and \code{b}). If supplied, thresholds
#'   are derived per chromosome.
#' @param b Numeric in [0,1]. Fallback background LD level used if \code{decay_tbl} is NULL.
#' @param slide_win_ld Integer. Sliding window size for LD estimation (passed to \code{get_el}).
#' @param q Numeric in (0,1). Quantile position along LD-decay curve used to define thresholds.
#' @param n_cores_ld Integer. Number of cores used for LD estimation (passed to \code{get_el}).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{CLS}: a list of clusters; each element is a character vector of SNP IDs.
#'   \item \code{ths}: a list with per-chromosome thresholds: \code{ld_ths} and \code{d_ths}.
#' }
#'
#' @examples
#' \dontrun{
#' cls <- SLC(gds, decay_tbl = LD_decay, b = LD_decay$b, slide_win_ld = 1000, q = 0.75)
#' sapply(cls$CLS, length)
#' }
#'
#' @export
SLC <- function(gds,
                geno = NULL,
                map = NULL,
                idx = NULL,
                decay_tbl = NULL,
                slide_win_ld = 1000,
                q = 0.75,
                n_cores_ld = 1) {
  
  t1 <- Sys.time()
  
  if (!is.finite(q) || q <= 0 || q >= 1) stop("`q` must be in (0,1).")
  
  ## Create temporary GDS if needed
  if (is.null(gds)) {
    if (is.null(geno) || is.null(map)) stop("If `gds` is NULL, both `geno` and `map` must be provided.")
    gds_path <- tempfile(fileext = ".gds")
    message("Generating temporary GDS object")
    
    gds <- create_gds_from_geno(geno, map, gds_path)
    
    on.exit({ SNPRelate::snpgdsClose(gds); unlink(gds_path) }, add = TRUE)
  }
  b <- decay_tbl$b
  
  ids <- read_gds_ids(gds)
  
  ## Determine which SNPs to cluster
  if (!is.null(idx)) {
    keep_ids <- ids$snp_id[idx]
    chrs <- unique(ids$snp_chr[idx])
  } else {
    keep_ids <- ids$snp_id
    chrs <- unique(ids$snp_chr)
  }
  
  ## Threshold containers (per chromosome)
  ld_ths_by_chr <- setNames(rep(NA_real_, length(chrs)), as.character(chrs))
  d_ths_by_chr  <- setNames(rep(NA_real_, length(chrs)), as.character(chrs))
  
  CLS_by_chr <- vector("list", length(chrs))
  names(CLS_by_chr) <- as.character(chrs)
  #ch <- chrs[1]
  for (ch in chrs) {
    cat(ch, "..")
    
    chr_idx <- which(ids$snp_chr == ch & ids$snp_id %in% keep_ids)
    chr_snp_ids <- ids$snp_id[chr_idx]
    
    ## Edge list with LD info (assumed columns: SNP1, SNP2, pos1, pos2, r2, dist)
    el <- get_el(gds, slide_win_ld = slide_win_ld, idx = chr_idx, n_cores = n_cores_ld)
    
    ## Derive thresholds
    if (is.null(decay_tbl)) {
      d_ths  <- Inf
      ld_ths <- b + (1 - b) * (1 - q)
    } else {
      summ <- decay_tbl$summary[decay_tbl$summary$Chr == ch, ]
      if (nrow(summ) != 1) stop("decay_tbl$summary must contain exactly one row per chromosome (missing/duplicate for Chr = ", ch, ").")
      
      d_ths  <- as.integer(d_from_rho(summ$a, rho = q))
      ld_ths <- summ$b + (1 - summ$b) * (1 - q)
    }
    
    ld_ths_by_chr[as.character(ch)] <- ld_ths
    d_ths_by_chr[as.character(ch)]  <- d_ths
    
    ## If no LD pairs available, everything is singleton
    if (is.null(el) || nrow(el) == 0L) {
      CLS_by_chr[[as.character(ch)]] <- as.list(chr_snp_ids)
      next
    }
    
    ## Filter edges by thresholds
    
    edges_dt <- el[r2 > ld_ths & d < d_ths, .(SNP1, SNP2)]
    
    ## If no edges pass thresholds -> all singleton clusters
    if (nrow(edges_dt) == 0L) {
      CLS_by_chr[[as.character(ch)]] <- as.list(chr_snp_ids)
      next
    }
    
    ## Build graph and connected components (single-linkage clusters)
    el_g <- as.matrix(edges_dt)
    g <- igraph::graph_from_edgelist(el_g, directed = FALSE)
    
    comps <- igraph::components(g)
    comp_members <- split(names(comps$membership), comps$membership)
    
    ## Add singleton SNPs not present in any edge
    in_graph <- unique(unlist(comp_members, use.names = FALSE))
    singletons <- chr_snp_ids[!chr_snp_ids %in% in_graph]
    
    cls <- c(comp_members, as.list(singletons))
    names(cls) <- seq_along(cls)
    
    CLS_by_chr[[as.character(ch)]] <- cls
  }
  
  ## Flatten across chromosomes
  n_loci <- length(unlist(CLS_by_chr, use.names = FALSE))
  CLS <- unlist(CLS_by_chr, recursive = FALSE, use.names = FALSE)
  
  red <- (n_loci - length(CLS)) / n_loci
  cat("\n", round(red * 100), "% reduction in number of loci (including singleton clusters)\n")
  t2 <- Sys.time()
  cat("Time elapsed ", difftime(t2, t1), "\n")
  
  list(
    CLS = CLS,
    ths = list(
      ld_ths_by_chr = ld_ths_by_chr,
      d_ths_by_chr  = d_ths_by_chr,
      q = q
    )
  )
}
