
#' Single-linkage clustering of SNPs based on LD and distance thresholds
#'
#' Partitions SNPs into clusters (outlier regions or LD clusters) using single-linkage
#' clustering based on LD and distance thresholds. For each chromosome, SNPs connected
#' by edges with \code{r2} above a threshold and distance below a threshold are grouped
#' into connected components.
#'
#' @param gds Optional open GDS object. If \code{NULL}, \code{geno} and \code{map} must
#'   be provided, and a temporary GDS will be created.
#' @param geno Optional genotype matrix (\eqn{n_\mathrm{ind} \times n_\mathrm{SNP}}).
#' @param map Optional SNP map; required if \code{geno} is used.
#' @param decay_tbl Optional list or object containing LD-decay summary and background LD.
#'   If provided, LD and distance thresholds are derived from it.
#' @param b Numeric, background LD level used as a fallback LD threshold if \code{decay_tbl}
#'   is not supplied.
#' @param slide_win_ld Integer, sliding window size for LD estimation.
#' @param q Numeric, quantile position along the LD-decay curve used to define thresholds.
#'
#' @return A list of clusters, where each element is a character vector of SNP IDs belonging
#'   to the same single-linkage cluster. Singleton clusters (isolated SNPs) are included.
#'
#' @examples
#' \dontrun{
#' cls <- SLC(gds, decay_tbl = LD_decay, b = LD_decay$b, slide_win_ld = 1000, q = 0.75)
#' lengths(cls)  # cluster sizes
#' }
#'
#' @export

SLC <- function(gds,geno=NULL,map=NULL,idx=NULL,decay_tbl=NULL,b=0.05,slide_win_ld=1000,q=0.75,n_cores_ld=1){
  
  t1 <- Sys.time()
  ## generate gds file
  if(is.null(gds)){
    gds_path = tempfile(fileext = ".gds")
    
    cat("Generating GDS object\n")
    gds <- create_gds_from_geno(geno, map,  gds_path)
    
    ## close at exit
    on.exit({ snpgdsClose(gds); unlink(gds_path) }, add = TRUE) 
  }
  
  ids <- read_gds_ids(gds)
  chrs <- unique(ids$snp_chr)
  
  if(!is.null(idx)){
    keep <- ids$snp_id[idx]  
    chrs <- unique(ids$snp_chr[idx])
  }else{
    keep <- ids$snp_id
  }
  
  
  #q=0.999
  #decay_tbl = LD_decay_3sp
  
  
  CLS <- list()
  #ch="Chr19"
  for(ch in chrs){
    cat(ch,"..")
    idx <- which(ids$snp_chr == ch & ids$snp_id %in% keep)
    
    el <- get_el(gds,slide_win_ld = slide_win_ld,idx = idx,n_cores = n_cores_ld)
    
    d_ths = ifelse(is.null(decay_tbl), Inf, as.integer(d_from_rho(decay_tbl$summary[Chr==ch,a],rho=q)))
    ld_ths = ifelse(is.null(decay_tbl), b + (1 - b)*(1-q),decay_tbl$summary[Chr==ch,b] + (1 - decay_tbl$summary[Chr==ch,b])*(1-q))
    
    el_g <- as.matrix(el[r2>ld_ths & d<d_ths,.(SNP1, SNP2)])
    el_g <- graph_from_edgelist(el_g, directed = FALSE)
    
    d_g <- decompose(el_g,min.vertices = 0)
    if(length(d_g)==1){
      cls <- list(unique(unlist(el[,.(SNP1,SNP2)])))
    }else{
      cls <- sapply(d_g,function(cl) V(cl)$name)  
    }
    
    
    
    in_cls <- unlist(cls)
    
    cls <- c(cls,ids$snp_id[idx][!ids$snp_id[idx] %in% in_cls])
    names(cls) <- 1:length(cls)
    CLS[[ch]] <-  cls
  }
  
  n_loci <- length(unlist(CLS))
  CLS <- unlist(CLS,recursive = FALSE)
  
  red <- (n_loci-length(CLS))/n_loci
  cat("\n",round(red*100),"% reduction in number of loci (including singleton clusters)\n")
  t2 <- Sys.time()
  cat("Time elapsed ",difftime(t2,t1))
  cat("\n")
  d_by_chr <- setNames(d_from_rho(decay_tbl$summary[,a], rho=q), decay_tbl$summary$Chr) # updated per draw below
  
  
  
  return(list(CLS=CLS,ths=list(ld_ths=ld_ths,d_by_chr=d_by_chr)))
  
}
