find_ORs <- function(gds, decay_tbl,el=NULL,outliers,rho_ld=0.999,rho_d=0.999){
  
  ids <- read_gds_ids(gds)
  
  idx <- which(ids$snp_id %in% outliers)
  
  if(is.null(el)){
    
    outl_chromosomes <- ids$snp_chr[idx]
    el_chr <- list()
    
    chrs <- unique(outl_chromosomes)
    
    ## only estimate LD for outliers within chromosomes
    for(ch in chrs){
      
      outl_ch <- outliers[outl_chromosomes==ch]
      
      el[[ch]] <- get_el(gds = gds,
                         idx =  which(ids$snp_id %in% outl_ch),
                         slide_win_ld = -1,
                         n_cores = n_cores)
    }
    
    el <- rbindlist(el)
  }
  
  
  d_th_by_chr <- setNames(as.integer(d_from_rho(decay_tbl$summary$a, rho_d)), decay_tbl$summary$Chr)
  ld_th_by_chr <- setNames(decay_tbl$summary$b + (1 - decay_tbl$summary$b) * (1 - rho_ld), decay_tbl$summary$Chr)
  el[,d_th:=d_th_by_chr[Chr1]]
  el[,ld_th:=ld_th_by_chr[Chr1]]
  
  ed <- el[r2 > ld_th & d < d_th & SNP1 %in% outliers & SNP2 %in% outliers ,.(SNP1, SNP2)]
  
  g <- graph_from_data_frame(ed, directed = FALSE)
  comps <- decompose(g)
  ORs <- lapply(comps, function(cc) V(cc)$name)
  
  
  ORs <- data.table(marker=unlist(ORs),OR=rep(seq_len(length(ORs)),sapply(ORs,length)),n_loci=rep(sapply(ORs,length),sapply(ORs,length)))
  return(ORs)
}
