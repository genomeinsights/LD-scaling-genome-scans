get_bg_ld <- function(gds, idx, n_sub = 5000, q = 0.95, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ids <- read_gds_ids(gds)
  n_snps <- length(idx)
  
  if (n_snps < 2L) stop("Not enough SNPs to estimate background LD.")
  
  ## sample SNPs, but ensure multiple chromosomes
  
  snp_pool <- unlist(lapply(split(seq_len(n_snps), ids$snp_chr[idx]), function(ix) {
    sample(ix, min(length(ix), ceiling(n_sub / length(unique(ids$snp_chr[idx])))))
  }))
  
  chr_pool <- ids$snp_chr[snp_pool]
  
  if (length(unique(chr_pool)) < 2L) {
    stop("Sampled SNPs fall on a single chromosome; increase n_sub.")
  }
  
  ## compute LD only once
  ld <- snpgdsLDMat(
    gds,
    snp.id = ids$snp_id[snp_pool],
    method = "r",
    slide = -1,
    verbose = FALSE
  )
  
  ## extract inter-chromosomal r^2 without melting full matrix
  r2 <- ld$LD^2
  chr <- chr_pool
  
  inter_idx <- outer(chr, chr, FUN = "!=")
  r2_inter <- r2[inter_idx]
  
  if (length(r2_inter) == 0L) {
    stop("No inter-chromosomal SNP pairs found.")
  }
  
  stats::quantile(r2_inter, probs = q, na.rm = TRUE)
}
