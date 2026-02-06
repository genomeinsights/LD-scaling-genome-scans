#' Compute SNP-level consistency scores (C) from OR draws
#'
#' @export
get_C <- function(draws, markers) {
  stopifnot("method" %in% names(draws), "ORs" %in% names(draws))
  
  ## number of draws per method
  n_draws <- draws[, .N, by = method]
  setnames(n_draws, "N", "n_draws")
  
  ## for each draw × method, get unique SNPs appearing in any OR
  snp_hits <- draws[
    lengths(ORs) > 0,
    .(marker = unique(unlist(ORs))),
    by = method
  ]
  
  ## count how many draws each SNP appears in (per method)
  snp_counts <- snp_hits[, .N, by = .(method, marker)]
  setnames(snp_counts, "N", "n_hits")
  
  ## compute C = hits / draws
  snp_counts <- merge(snp_counts, n_draws, by = "method")
  snp_counts[, C := n_hits / n_draws]
  
  ## cast to wide SNP × method matrix
  C_wide <- dcast(
    snp_counts,
    marker ~ method,
    value.var = "C",
    fill = 0
  )
  
  ## align to reference SNP order
  C_wide <- C_wide[match(markers, marker)]
  C_wide[is.na(C_wide)] <- 0
  
  ## drop marker column and prefix names
  C_mat <- as.data.table(C_wide[, -1])
  setnames(C_mat, paste0("C_", colnames(C_mat)))
  
  C_mat
}
