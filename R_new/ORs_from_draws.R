#' Generate parameter draws and outlier regions for multiple methods
#'
#' @export
ORs_from_draws <- function(gds,
                           decay_tbl,
                           ld_w_draws,
                           F_vals,
                           n_draws = 25,
                           n_inds = 125,
                           rho_OR_lim = list(min = 0.5, max = 1.0),
                           rho_ld_lim = NULL,
                           rho_d_lim  = NULL,
                           alpha_lim  = list(min = 1.31, max = 4.0),
                           lmin_lim   = list(min = 2, max = 20),
                           n_cores = 1) {
  
  t1 <- Sys.time()
  cat("Preparing data\n")
  
  F_vals <- as.matrix(F_vals)
  if (!is.numeric(F_vals)) stop("`F_vals` must be numeric.")
  if (nrow(F_vals) != length(read_gds_ids(gds)$snp_id)) {
    stop("`F_vals` must have one row per SNP in the GDS (in the same order).")
  }
  
  df2 <- n_inds - 2
  if (df2 <= 0) stop("`n_inds` must be >= 3.")
  
  ## q-values (FDR-adjusted) for original F
  q_orgs <- apply(F_vals, 2, function(Fval) {
    p <- stats::pf(Fval, df1 = 1, df2 = df2, lower.tail = FALSE)
    stats::p.adjust(p, "fdr")
  })
  
  cat("Starting draws\n")
  
  ORs <- get_ORs_for_draw(
    gds        = gds,
    F_vals     = F_vals,
    q_orgs     = q_orgs,
    decay_tbl  = decay_tbl,
    ld_w_draws = ld_w_draws,
    n_inds     = n_inds,
    n_draws    = n_draws,
    rho_OR_lim = rho_OR_lim,
    rho_ld_lim = rho_ld_lim,
    rho_d_lim  = rho_d_lim,
    alpha_lim  = alpha_lim,
    lmin_lim   = lmin_lim,
    n_cores    = n_cores
  )
  
  t2 <- Sys.time()
  print(difftime(t2, t1))
  ORs
}

#' @keywords internal
get_ORs_for_draw <- function(gds,
                             F_vals,
                             q_orgs,
                             decay_tbl,
                             ld_w_draws,
                             n_inds,
                             n_draws,
                             rho_OR_lim,
                             rho_ld_lim = NULL,
                             rho_d_lim  = NULL,
                             alpha_lim,
                             lmin_lim,
                             n_cores = 1) {
  
  ids <- read_gds_ids(gds)
  
  ## sanity: ld_w_draws must have draw, rho_w, d_th, ld_w(list)
  req <- c("draw", "rho_w", "d_th", "ld_w")
  if (!all(req %in% names(ld_w_draws))) stop("`ld_w_draws` must contain: draw, rho_w, d_th, ld_w.")
  
  ## draws per LD-window draw
  out <- list()
  # dr  <- sort(unique(ld_w_draws$draw))[1]
  for (dr in sort(unique(ld_w_draws$draw))) {
    cat(dr, "..")
    
    ## ld_w table for this draw: one row per SNP, columns SNP + median_r2
    
    ld_w_tbl <- data.table::rbindlist(ld_w_draws[draw == dr, ld_w])
    if (!all(c("SNP", "median_r2") %in% names(ld_w_tbl))) {
      stop("Each ld_w element must contain columns: SNP, median_r2.")
    }
    
    ## align ld_w to SNP order in GDS
    ld_w <- ld_w_tbl[match(ids$snp_id, SNP), median_r2]
    ## LD-scaled q-values (your function)
    q_primes <- scale_F_with_ldw(
      F_vals = F_vals,
      ld_w   = ld_w,
      n_rep  = 10,
      n_inds = n_inds,
      full   = FALSE
    )
    colnames(q_primes) <- paste0(colnames(q_primes), "_prime")
    
    #cor(map$max_LD_with_QTN,q_primes[,2])^2
    
    q_vals <- cbind(q_orgs, q_primes)
    
    ## “Joint” = per-SNP min q across methods (as in your code)
    q_min <- apply(q_vals, 1, min)
    q_vals <- as.matrix(cbind(q_vals, joint = q_min))
    
    ## prefilter outliers using the *least stringent* alpha threshold you might draw
    alpha_prefilter <- 1 / (10^alpha_lim$min)
    all_outliers <- ids$snp_id[q_min < alpha_prefilter]
    
    ## precompute LD edges among candidate SNPs (within each chr), once per dr
    if (length(all_outliers) > 0) {
      outl_chr <- ids$snp_chr[match(all_outliers, ids$snp_id)]
      chrs <- unique(outl_chr)
      
      el_list <- vector("list", length(chrs))
      names(el_list) <- as.character(chrs)
      
      for (ch in chrs) {
        outl_ch <- all_outliers[outl_chr == ch]
        chr_idx <- which(ids$snp_id %in% outl_ch)
        
        ## all-pairs LD within candidate set (be careful: can be large)
        el_list[[as.character(ch)]] <- get_el(
          gds = gds,
          idx = chr_idx,
          slide_win_ld = -1,
          n_cores = n_cores
        )
      }
      
      el_chr <- data.table::rbindlist(el_list, use.names = TRUE, fill = TRUE)
    } else {
      el_chr <- NULL
    }
    
    ## generate n_draws OR parameter draws for this ld_w draw
    draw_ORs <- data.table::rbindlist(
      #x <- 1
      parallel::mclapply(seq_len(n_draws), function(x) {
        cbind(
          draw = x,
          get_ORs_for_draw_w(
            ids        = ids,
            q_vals     = q_vals,
            decay_tbl  = decay_tbl,
            el_chr     = el_chr,
            rho_OR_lim = rho_OR_lim,
            rho_ld_lim = rho_ld_lim,
            rho_d_lim  = rho_d_lim,
            alpha_lim  = alpha_lim,
            lmin_lim   = lmin_lim
          )
        )
      }, mc.cores = n_cores),
      use.names = TRUE,
      fill = TRUE
    )
    
    out[[as.character(dr)]] <- cbind(
      draw_w = dr,
      ld_w_draws[draw == dr, .(rho_w = rho_w[1], mean_d_th = mean(d_th))],
      draw_ORs
    )
  }
  
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

#' @keywords internal
get_ORs_for_draw_w <- function(ids,
                               q_vals,
                               decay_tbl,
                               el_chr,
                               rho_OR_lim,
                               rho_ld_lim = NULL,
                               rho_d_lim  = NULL,
                               alpha_lim,
                               lmin_lim) {
  
  ## ---- draw thresholds ----
  if (!is.null(rho_OR_lim)) {
    rho_d  <- runif(1, rho_OR_lim$min, rho_OR_lim$max)
    rho_ld <- rho_d
  } else {
    stopifnot(!is.null(rho_ld_lim), !is.null(rho_d_lim))
    rho_ld <- runif(1, rho_ld_lim$min, rho_ld_lim$max)
    rho_d  <- runif(1, rho_d_lim$min,  rho_d_lim$max)
  }
  
  l_min <- sample(seq.int(lmin_lim$min, lmin_lim$max), 1)
  alpha <- 1 / 10^(runif(1, alpha_lim$min, alpha_lim$max))
  
  methods <- colnames(q_vals)
  
  ## ---- no LD edges at all ----
  if (is.null(el_chr) || nrow(el_chr) == 0) {
    return(data.table::data.table(
      method = methods,
      alpha  = alpha,
      l_min  = l_min,
      rho_d  = rho_d,
      rho_ld = rho_ld,
      mean_d_th_OR  = NA_real_,
      mean_ld_th_OR = NA_real_,
      ORs = replicate(length(methods), list(), simplify = FALSE)
    ))
  }
  
  ## ---- attach per-chromosome thresholds ----
  d_th_by_chr <- setNames(
    d_from_rho(decay_tbl$summary$a, rho = rho_d),
    decay_tbl$summary$Chr
  )
  
  ld_th_by_chr <- setNames(
    decay_tbl$summary$b + (1 - decay_tbl$summary$b) * (1 - rho_ld),
    decay_tbl$summary$Chr
  )
  
  el_chr[, d_th  := d_th_by_chr[Chr1]]
  el_chr[, ld_th := ld_th_by_chr[Chr1]]
  
  mean_d_th  <- mean(el_chr$d_th,  na.rm = TRUE)
  mean_ld_th <- mean(el_chr$ld_th, na.rm = TRUE)
  
  ## ---- generate ORs per method ----
  out <- vector("list", length(methods))
  
  # i <- 1
  for (i in seq_along(methods)) {
    me <- methods[i]
    
    outl_me <- ids$snp_id[q_vals[, i] < alpha]
    
    ed <- el_chr[
      r2 > ld_th & d < d_th &
        SNP1 %in% outl_me & SNP2 %in% outl_me,
      .(SNP1, SNP2)
    ]
    
    if (nrow(ed) == 0L) {
      ors <- list()
    } else {
      g <- igraph::graph_from_data_frame(ed, directed = FALSE)
      comps <- igraph::components(g)
      
      ors <- split(names(comps$membership), comps$membership)
      ors <- ors[vapply(ors, length, integer(1)) >= l_min]
      names(ors) <- seq_along(ors)
    }
    
    out[[i]] <- data.table::data.table(
      method = me,
      alpha  = alpha,
      l_min  = l_min,
      rho_d  = rho_d,
      rho_ld = rho_ld,
      mean_d_th_OR  = mean_d_th,
      mean_ld_th_OR = mean_ld_th,
      ORs = list(ors)
    )
  }
  
  data.table::rbindlist(out)
}
