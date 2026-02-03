
#' Build an LD edge list for a subset of SNPs
#'
#' Computes pairwise LD (\eqn{r^2}) for a specified subset of SNPs in a GDS file
#' and returns an edge list with genomic coordinates and distances.
#'
#' @param gds An open GDS object.
#' @param idx Integer vector of SNP indices (1-based, referring to the GDS SNP index).
#' @param slide_win_ld Integer. If positive, LD is estimated using a sliding window of
#'   this size; if \code{<= 0}, all pairwise LD values are computed.
#' @param n_cores Integer, number of threads to use in LD computation.
#'
#' @return A \code{data.table} with columns:
#'   \itemize{
#'     \item \code{Var1}, \code{Var2}: internal indices of SNP pairs.
#'     \item \code{r2}: pairwise LD (\eqn{r^2}).
#'     \item \code{Chr1}, \code{Chr2}: chromosomes of SNP1 and SNP2.
#'     \item \code{pos1}, \code{pos2}: positions of SNP1 and SNP2.
#'     \item \code{SNP1}, \code{SNP2}: SNP IDs.
#'     \item \code{d}: physical distance between SNPs in base pairs.
#'   }
#'
#' @examples
#' \dontrun{
#' idx <- 1:1000
#' el  <- get_el(gds, idx, slide_win_ld = 1000)
#' head(el)
#' }
#'
#' @export

get_el <- function(gds, idx, slide_win_ld=1000, n_cores=1){
  ids <- read_gds_ids(gds)
  if(slide_win_ld>0){
    
    ldmat <- snpgdsLDMat(gds, snp.id = ids$snp_id[idx], method = "r", slide = slide_win_ld,verbose=FALSE,num.thread = n_cores)
    el <- as.data.table(reshape2::melt(ldmat$LD^2,value.name="r2"))
    el[,Var1:=Var1+Var2]
    el <- rbind(el,el[,.(Var1=Var2,Var2=Var1,r2)]) ## this is necessary when sliding window is used!!
    el <- el[!is.nan(r2)]
    el[,Chr1:=ids$snp_chr[idx][Var1]]
    el[,Chr2:=ids$snp_chr[idx][Var2]]
    el[,pos1:=ids$snp_pos[idx][Var1]]
    el[,pos2:=ids$snp_pos[idx][Var2]]
    el[,SNP1:=ids$snp_id[idx][Var1]]
    el[,SNP2:=ids$snp_id[idx][Var2]]
    el[,d:=abs(pos1-pos2)]
    return(el)
    
  }else{
    ldmat <- snpgdsLDMat(gds, snp.id = ids$snp_id[idx], method = "r", slide = -1,verbose=FALSE,num.thread = n_cores)
    el <- as.data.table(reshape2::melt(ldmat$LD^2,value.name="r2"))
    el <- el[!is.nan(r2)]
    el[,Chr1:=ids$snp_chr[idx][Var1]]
    el[,Chr2:=ids$snp_chr[idx][Var2]]
    el[,pos1:=ids$snp_pos[idx][Var1]]
    el[,pos2:=ids$snp_pos[idx][Var2]]
    el[,SNP1:=ids$snp_id[idx][Var1]]
    el[,SNP2:=ids$snp_id[idx][Var2]]
    el[,d:=abs(pos1-pos2)]
    
    return(el)
  }
}
