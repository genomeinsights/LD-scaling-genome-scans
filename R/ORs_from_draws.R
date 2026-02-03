
#' Generate parameter draws and outlier regions for multiple methods
#'
#' High-level driver that:
#' \enumerate{
#'   \item Converts F-statistics (\code{F_vals}) to q-values for each method.
#'   \item Combines original and LD-scaled statistics via \code{\link{scale_F_with_ldw}}.
#'   \item For each LD-window draw, generates multiple random draws of clustering
#'         parameters and constructs outlier regions (ORs).
#' }
#'
#' @param gds An open GDS object.
#' @param decay_tbl LD-decay summary list from \code{\link{get_LD_decay_data}}.
#' @param ld_w_draws Object with precomputed \code{ld_w} values and \code{rho_w} per draw.
#' @param F_vals Numeric vector or matrix of F-statistics (one column per method).
#' @param n_draws Integer, number of OR draws per LD-window setting (default \code{25}).
#' @param rho_OR_lim List with components \code{min} and \code{max} for \eqn{\rho} in OR thresholds.
#' @param alpha_lim List with components \code{min} and \code{max} for the -log10(\eqn{\alpha})
#'   range sampled.
#' @param max_lmin Integer, maximum minimum OR size (\code{l_min}) considered in draws.
#' @param n_cores Integer, number of cores for parallel computation.
#'
#' @return A \code{data.table} where each row corresponds to a particular random draw of
#'   parameter values and method, with columns describing the thresholds used and a list
#'   of ORs (\code{ORs}) identified under those settings.
#'
#' @examples
#' \dontrun{
#' draws <- get_draws(gds, decay_tbl, ld_w_draws, F_vals = F_mat,
#'                    n_draws = 25,
#'                    rho_OR_lim = list(min = 0.9, max = 1.0),
#'                    alpha_lim  = list(min = 0.602, max = 4.0),
#'                    max_lmin   = 20,
#'                    n_cores    = 8)
#' }
#'
#' @export

ORs_from_draws <- function(gds, 
                      decay_tbl, 
                      ld_w_draws,
                      F_vals,
                      n_draws=25,
                      n_inds=125,
                      rho_OR_lim = list(min=0.5, max=1.0),
                      alpha_lim  = list(min=1.31, max=4.0),
                      lmin_lim   = list(min=2, max=20),
                      n_cores = 1) {
  t1 <- Sys.time()
  cat("Preparing data\n")
  F_vals <- as.matrix(F_vals)
  #ids <- read_gds_ids(gds)
  
  
  #chrs <- unique(ids$snp_chr)  
  
  #b <- decay_tbl$b
  
  ## get q-values from original F-values
  df2 = n_inds-2
  q_orgs <- apply(F_vals,2,function(Fval){
    p <- pf(Fval, df1=1, df2, lower.tail = FALSE)
    q <- p.adjust(p,"fdr")
  })
  
  ## get F´ from F_vals
  
  cat("Starting draws\n")
  
  ORs <- get_ORs_for_draw(gds,
                          decay_tbl,
                          ld_w_draws=ld_w_draws,
                          F_vals=F_vals,
                          q_orgs=q_orgs,
                          n_draws=n_draws,
                          n_inds=n_inds,
                          rho_OR_lim,
                          alpha_lim=alpha_lim,
                          lmin_lim=lmin_lim,
                          n_cores=n_cores)
  t2 <- Sys.time()
  print(difftime(t2,t1))
  return(ORs)
  
}


#' Generate outlier regions across multiple LD window draws
#'
#' Internal function that, for each LD-window draw (i.e. a given \code{rho_w} and
#' corresponding \code{ld_w} vector), computes LD-scaled q-values, identifies outliers,
#' and constructs ORs for multiple random combinations of clustering parameters.
#'
#' @param gds An open GDS object.
#' @param decay_tbl LD-decay summary object as returned by \code{\link{get_LD_decay_data}}.
#' @param ld_w_draws A \code{data.table} or similar structure containing per-draw LD-window
#'   information, including per-SNP \code{ld_w} values and the corresponding \code{rho_w}
#'   or window sizes.
#' @param q_orgs Matrix of original q-values (FDR-adjusted), one column per method.
#' @param n_draws Integer, number of OR draws to generate per LD-window draw.
#' @param rho_OR_lim List specifying the minimum and maximum \eqn{\rho} values for deriving
#'   LD and distance thresholds when forming ORs.
#' @param max_lmin Integer, maximum minimum cluster size (\code{l_min}) for ORs.
#' @param alpha_lim List specifying the range for significance thresholds on the -log10 scale.
#' @param n_cores Integer, number of cores for parallel computation.
#'
#' @return A \code{data.table} with one row per parameter combination and method, including
#'   the drawn thresholds and the list-column \code{ORs} containing SNP clusters.
#'
#' @keywords internal


get_ORs_for_draw <- function(gds,
                             F_vals,
                             q_orgs,
                             decay_tbl,
                             ld_w_draws,
                             n_inds,
                             n_draws,
                             rho_OR_lim,
                             rho_ld_lim,
                             rho_d_lim,
                             alpha_lim,
                             lmin_lim,
                             n_cores=1){
  
  #rm(ids)
  
  draws <- list()
  ids <- read_gds_ids(gds)
  
 
  
 #dr <- 1
  for(dr in unique(ld_w_draws$draw)){
    cat(dr,"..")
    ld_w <- rbindlist(ld_w_draws[draw==dr,ld_w])
    ld_w <- ld_w[match(ids$snp_id,SNP1),median_r2]
    q_primes <- scale_F_with_ldw(F_vals = F_vals,  ld_w = ld_w,n_rep = 10,n_inds = n_inds, full = FALSE)
    colnames(q_primes) <- paste(colnames(q_primes),"prime",sep="_")
    q_vals <- cbind(q_orgs,q_primes)
    
    #plot(-log10(q_vals[1:10000,1]),-log10(q_vals[1:10000,2]))
    ## prepare edge lists for all chromosomes in andvance
    q_min <- apply(q_vals,1,min)
    #table(q_min<0.05)
    q_vals <- cbind(q_vals,joint=q_min)
    
    all_outliers <- ids$snp_id[q_min<1/10^alpha_lim$min]
    
    
    if(length(all_outliers)>0){
      outl_chromosomes <- ids$snp_chr[ids$snp_id %in% all_outliers]

      chrs <- unique(outl_chromosomes)
      
      el_chr <- list()
      
      #all_outliers
      #ch = "Chr7"
      for(ch in chrs){
        outl_ch <- all_outliers[outl_chromosomes==ch]
        idx <- which(ids$snp_id %in% outl_ch)
        el_chr[[ch]] <- get_el(gds = gds,idx = idx,slide_win_ld = -1,n_cores = n_cores)
      }
      
      el_chr <- rbindlist(el_chr)
      
      
     
    }else{
      el_chr <- NULL
    }
    
    ## generate n_draws per ld_w vector
    
    draw_ORs <- rbindlist(mclapply(1:n_draws, function(x) {
      
      cbind(draw=x,get_ORs_for_draw_w(ids,
                                      q_vals,
                                      decay_tbl,
                                      el_chr,
                                      rho_OR_lim,
                                      rho_ld_lim,
                                      rho_d_lim,
                                      alpha_lim,
                                      lmin_lim
      ))
    },mc.cores=n_cores))
    
    draws[[dr]] <- cbind(draw_w=dr,ld_w_draws[draw==dr,.(rho_w=rho_w[1],mean_d_th=mean(d_th))],draw_ORs)   
  }
  rbindlist(draws)
}



#' Generate outlier regions (ORs) for a single parameter draw
#'
#' Internal helper used to generate outlier regions (ORs) for a given set of threshold
#' parameters, given q-values, LD information, and precomputed edge lists per chromosome.
#'
#' @param q_vals Matrix of q-values (FDR-adjusted p-values), with one row per SNP and
#'   one column per method.
#' @param decay_tbl LD-decay table with summary parameters (including distance and LD thresholds).
#' @param el_chr Named list of per-chromosome edge lists as produced by \code{\link{get_el}}.
#' @param chrs Vector of chromosome identifiers corresponding to SNPs considered.
#' @param rho_OR_lim List with components \code{min} and \code{max} specifying the range
#'   of \eqn{\rho} values used to draw LD and distance thresholds for OR definition.
#' @param max_lmin Integer, maximum value of the minimum OR size (\code{l_min}) sampled
#'   at random.
#' @param alpha_lim List with components \code{min} and \code{max} specifying the range
#'   for the significance threshold on the -log10 scale (for \code{alpha} draws).
#'
#' @return A \code{data.table} with one row per method, containing the drawn parameters
#'   and the corresponding list of ORs (clusters of SNPs).
#'
#' @keywords internal

get_ORs_for_draw_w <- function(ids,q_vals,decay_tbl,el_chr,rho_OR_lim,rho_ld_lim,rho_d_lim,alpha_lim,lmin_lim){
  
  ## draw values from uniform distr.
  
  
    if(is.null(rho_OR_lim)){
      rho_ld <- runif(1,rho_ld_lim$min,rho_ld_lim$max)
      rho_d <- runif(1,rho_d_lim$min,rho_d_lim$max)  
    }else{
      rho_d <- rho_ld <- runif(1,rho_OR_lim$min,rho_OR_lim$max)  
    }
    
    #l_min = 1
    
    l_min <- sample(lmin_lim$min:lmin_lim$max,1)
    alpha <- 1/10^(runif(1,alpha_lim$min,alpha_lim$max))
    
    if(length(el_chr)>0){
      ## get thresholds based on drawn rho values
      d_th_by_chr <- setNames(as.integer(d_from_rho(decay_tbl$summary$a, rho_d)), decay_tbl$summary$Chr)
      ld_th_by_chr <- setNames(decay_tbl$summary$b + (1 - decay_tbl$summary$b) * (1 - rho_ld), decay_tbl$summary$Chr) 
      el_chr[,d_th:=d_th_by_chr[Chr1]]
      el_chr[,ld_th:=ld_th_by_chr[Chr1]]
   }
    
    
    
    #for all columns corresponding to F-values from different methods
    # i <- 1
    
    if (length(el_chr) == 0){
      ORs <- data.table(method=colnames(q_vals),mean_d_th_OR=NA,mean_ld_th_OR=NA,ORs=list())
      
    } else {
      ORs <- list()
      for(i in 1:ncol(q_vals)){
        
        
        outl_me <- ids$snp_id[q_vals[,i]<alpha]
        ed <- el_chr[r2 > ld_th & d < d_th & SNP1 %in% outl_me & SNP2 %in% outl_me ,.(SNP1, SNP2)]
        
        if (nrow(ed) < 1){
          ORs[[colnames(q_vals)[i]]] <- data.table(method=colnames(q_vals)[i],mean_d_th_OR=mean(el_chr$d_th),mean_ld_th_OR=mean(el_chr$d_th),ORs=list())
          next
        } 
        
        
        g <- graph_from_data_frame(ed, directed = FALSE)
        comps <- decompose(g)
        ors <- lapply(comps, function(cc) V(cc)$name)
        
        ## remove small clusters
        ors <- ors[sapply(ors, length) >= l_min] 
        
        names(ors) <- seq_len(length(ors))
        if(length(ors)>0){
          ORs[[colnames(q_vals)[i]]] <- data.table(method=colnames(q_vals)[i],mean_d_th_OR=mean(el_chr$d_th),mean_ld_th_OR=mean(el_chr$ld_th),ORs=list(ors)) 
        }else{
          ORs[[colnames(q_vals)[i]]] <- data.table(method=colnames(q_vals)[i],mean_d_th_OR=mean(el_chr$d_th),mean_ld_th_OR=mean(el_chr$ld_th),ORs=list())  
        }
        
      }
      ORs <- rbindlist(ORs)
    }
    
    
    
    return(data.table(alpha,l_min,rho_d=rho_d,rho_ld=rho_ld,ORs))
 
  
}
