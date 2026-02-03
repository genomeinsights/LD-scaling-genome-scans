######################################################
##  Generates the necessary command line arguments used to run Nemo
##   simulations in parallel through sbatch or similar on a cluster.
##
##  Relies on template/base file sim_base.ini where raw code for Nemo simulations
##   are saved. The code below modifies and creates a copy of sim_base.ini for
##   a given parameter combination
######################################################


#----------------------------------------------------------
#  Load packages and set wd
#----------------------------------------------------------

library(data.table,quietly = TRUE, verbose = FALSE,warn.conflicts=FALSE)
library(parallel,quietly = TRUE, verbose = FALSE,warn.conflicts=FALSE)
library(ggplot2,quietly = TRUE, verbose = FALSE,warn.conflicts=FALSE) ## only for plotting
library(cowplot) ## only for plotting
library(sp)
library(gstat)
setwd("./Nemo")
#----------------------------------------------------------
#  Define functions to generate migration rates
#----------------------------------------------------------


kernel = function(d, a, c) { (c/(2*a*gamma(1/c)))*exp(-(d/a)^(c)) }

## helper functions for euclidean distance
c_ij = function(p1, p2,ncol) {
  
  c1 = (p1-1) %% ncol
  c2 = (p2-1) %% ncol
  
  return(c2 - c1)
}

r_ij = function(p1, p2,ncol) {
  
  r1 = floor( (p1-1) / ncol)
  r2 = floor( (p2-1) / ncol)
  
  return(r2 - r1) 
}

## Builds 2-dimensional distance matrix of dimensions nrow x ncol ####
buildDistanceMatrix = function(nrow, ncol) {
  
  num_patch = nrow * ncol
  
  distmat = matrix(0, nrow=num_patch, ncol=num_patch)
  
  for(i in 1:num_patch) {
    
    for(j in i:num_patch) {
      
      distmat[i,j] = sqrt( c_ij(i,j,ncol)^2 + r_ij(i,j,ncol)^2)
    }
  }
  
  lowtri = t(distmat)
  
  diag(lowtri) <- 0
  
  distmat = distmat + lowtri
  
  return(distmat)
}

# ------------------------------- #
# Basic nemo matrix-writing function
# It writes a matrix (mat) in the nemo matrix format to a text file (outfile)
# ------------------------------- #
write.matrix.nemo = function(mat, outfile) 
{
  rows = dim(mat)[1]
  cols = dim(mat)[2]
  
  cat("{",file=outfile)
  for(i in 1:rows) {
    
    cat("{",file=outfile, append=TRUE)
    cat(mat[i,],sep=",", file=outfile, append=TRUE)
    cat("}\n",file=outfile, append=TRUE)
  }
  
  cat("}\n",file=outfile, append=TRUE)
}

# ------------------------------- #
# builds distance and dispersal matrix and writes dispersal matrix to file
# The colSums are adjusted to 1 by changing the diagonal. This means that if gene flow is too high, 
# the diagonal becomes negative and a warning is given (in that case try decreasing a and/or increasing c)
# if return=TRUE outputs also the distance and dispersal matrix to check outcome
# side = number of patches for each side of the matrix (matrix always symmetrical)
# a = mean of the distance dispersal kernel
# c = exponent of the distance dispersal kernel
# d_thresh = minimum dispersal rate under which dispersal rate is set to 0
# ------------------------------- #


buildDispersalMatrix <- function(side,a=0.5,c=2,d_thresh=1e-7,file_name="./disp_mat.txt",return=FALSE){
  distance_matrix <- buildDistanceMatrix(nrow=side,ncol=side)
  
  disp_matrix <- kernel(distance_matrix,a,c) 
  
  disp_matrix[disp_matrix<d_thresh] <- 0
  
  for(j in 1:ncol(disp_matrix)){
    
    vals <- disp_matrix[,j]
    disp_matrix[j,j] <- 1-sum(vals[-j])
    
  }
  
  if(any(disp_matrix<0)){
    cat("Warning: diagonals<0") ## with this implementation, if gene flow is high, there diagonals can be negative
  }else{
    write.matrix.nemo(disp_matrix,outfile=file_name)
    if(return) return(list(dist=distance_matrix,disp=disp_matrix))
  } 
  
}

#----------------------------------------------------------
# Define parameters
#----------------------------------------------------------

if(file.exists("cmds.sh")) system("rm cmds.sh")

params <- expand.grid(chr=1:10,
                      V=c(0.5,1,2),
                      r=10,
                      a=0.5,
                      c=c(1,1.5,2),
                      rep=1:10)

side = 48

params <- split(as.matrix(params),row(params))

#----------------------------------------------------------
# Write dispersal matrices (the same three are used for all simulations)
#----------------------------------------------------------


if(!file.exists("disp_mat_1.txt")){
  matrices <-   buildDispersalMatrix(side, 0.5, 1,return=TRUE,file_name = "disp_mat_1.txt") 
  saveRDS(matrices,"matrices_1.txt")
}
if(!file.exists("disp_mat_1.5.txt")){
  matrices <-   buildDispersalMatrix(side, 0.5, 1.5,return=TRUE,file_name = "disp_mat_1.5.txt") 
  saveRDS(matrices,"matrices_1.5.txt")
}
if(!file.exists("disp_mat_2.txt")){
  matrices <-   buildDispersalMatrix(side, 0.5, 2,return=TRUE,file_name = "disp_mat_2.txt")
  saveRDS(matrices,"matrices_2.txt")
}

#----------------------------------------------------------
# Generate environments
#----------------------------------------------------------

scale_01 <- function(x) (x-min(x))/max((x-min(x)))

side = 48

# 1. Create a grid of points in 2D space
x.range <- 1:side
y.range <- side:1 # needs to be reversed because counting starts from the top in nemo
grid <- expand.grid(x = x.range, y = y.range)
coordinates(grid) <- ~x + y

vmodel <- vgm(psill = 10, model = "Exp", range = 50,nugget=0.01)

## examples
plots <- lapply(1:9,function(rep){
  
  sim <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 5, model = vmodel, nmax = 20)
  simulated_data <- predict(sim, newdata = grid, nsim = 1)
  
  
  dt <- data.table(simulated_data["sim1"]@coords,Z=unlist(simulated_data["sim1"]@data))
  
  p1 <- ggplot(dt, aes(x,y,fill=scale_01(Z))) +
    geom_tile() +
    scale_fill_viridis_c(option="turbo") +
    theme_bw() +
    theme(aspect.ratio = 1,
          legend.position = "none",
          strip.background = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()
    )
  return(p1)
})

plot_grid(plotlist = plots)

## now generate one for nemo  
lapply(1:10,function(rep){
  sim <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = 5, model = vmodel, nmax = 20)
  simulated_data <- predict(sim, newdata = grid, nsim = 1)
  dt <- data.table(simulated_data["sim1"]@coords,env=scale_01(unlist(simulated_data["sim1"]@data)))
  
  ## check results, redo if not happy
  pdf(paste("env_",rep,".pdf"),width = 5,height = 5)
  print(ggplot(dt, aes(x,y,fill=env)) +
          geom_tile() +
          scale_fill_viridis_c(option="turbo") +
          theme_bw() +
          theme(aspect.ratio = 1,
                legend.position = "none",
                strip.background = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()
          ))
  dev.off()
  
  
  
  #write infile for nemo
  if(!file.exists(paste0("env_",rep,".txt"))){
    write.table(paste0("{",paste(paste0("{",dt$env,"}"),collapse=""),"}"),paste0("env_",rep,".txt"),col.names = FALSE,quote = FALSE,row.names = FALSE)
  }
  
})


#----------------------------------------------------------
# Generate recombination maps and allelic values and save data/command line arguments for running on cluster
#----------------------------------------------------------

lapply(params,function(par){
  
  # ------------------------------- #
  #  Set up parameters, folders and files
  # ------------------------------- #
  chr = as.numeric(par[1])
  V=as.numeric(par[2])
  L=as.numeric(par[3])
  nQTN=200 ## total QTN 100/chr
  nSNPs=2e5 ## total nSNPs
  a = as.numeric(par[4])
  c = as.numeric(par[5])
  rep =  as.numeric(par[6])
  
  outfolder <- paste("chr",chr, "_V",V, "_c",c,"_rep",rep, sep="")
  
  message("Working on",outfolder)
  
  system(paste0("cp ","run_sim.sh ",paste("run_sim",outfolder,sep="_"),".sh"))
  system(paste0("sed -i '' \"s/","name","/",outfolder,"/g\" ", paste0(paste("run_sim",outfolder,sep="_"),".sh")))
  
  if(!dir.exists(outfolder)) dir.create(outfolder)
  
  system(paste("cp sim_base.ini ", paste0(outfolder,".ini")))
  system(paste(paste0("cp env_",rep,".txt "), paste0(outfolder,"/env.txt")))
  
  system(paste(paste0("cp disp_mat_",c,".txt "), paste0(outfolder,"/disp_mat.txt")))
  
  write(paste0(paste("sbatch run_sim",outfolder,sep="_"),".sh"),file="cmds.sh",append=TRUE)
  
  system(paste0("sed -i '' \"s/","root_directory","/",outfolder,"/g\" ", paste0(outfolder,".ini")))
  
  # ------------------------------- #
  #  Generate recombination maps
  # ------------------------------- #
  rec_map <- fread("./recombinatin_maps.txt")[CHR==chr] ## recombination map comes from external file
  
  ## use mean position between male and female
  rec_map[,MEAN_POS := apply(rec_map[,.(MALE_POS, FEMALE_POS)], 1, mean)]
  rec_map[,rec_rate := c((rec_map$MEAN_POS[-1]-rec_map$MEAN_POS[-nrow(rec_map)])/(rec_map$POS[-1]-rec_map$POS[-nrow(rec_map)]),0)*10^6]
  rec_map[,cM:=MEAN_POS] 
  rec_map[,POS:=as.integer(POS)]
  
  
  system(paste0("sed -i '' \"s/m_r/",signif(L/max(rec_map$POS),3),"/g\" ", paste0(outfolder,".ini")))
  system(paste0("sed -i '' \"s/sel_var/",V,"/g\" ", paste0(outfolder,".ini")))
  
  # sample positions
  QTN_pos <- sort(sample(1:max(rec_map$POS),nQTN/2))
  
  # estimate density of qtn and draw neutral SNPs based on the inverse CDF based on kernel density estimation
  d <- density(QTN_pos,bw=100,from = 0,to=max(rec_map$POS))  # KDE estimate
  
  dx <- diff(d$x)[1]  # uniform spacing
  cdf_vals <- cumsum(d$y) * dx
  cdf_vals <- cdf_vals / max(cdf_vals)  # Normalize
  
  inverse_cdf <- suppressWarnings(approxfun(cdf_vals, d$x, rule = 1)) # inverse CDF interpolator
  
  n_samples <- nSNPs/2-nQTN/2   # results in 1/3 of uniformly distributed SNPs 
  
  u <- runif(n_samples)         # uniform samples in [0, 1]
  bp <- inverse_cdf(u)          # transform via inverse CDF
  bp <- round(sort(c(na.omit(bp),runif(length(which(is.na(bp))),min=0,max=max(rec_map$POS))))) ## sort and replace NAs with random positions
  bp[which(duplicated(bp))] <- round(runif(length(which(duplicated(bp))),min=0,max=max(rec_map$POS)))
  bp <- sort(bp)
  
  cat("sampling positions\n")
  while(any(duplicated(bp))){
    u <- runif(length(which(duplicated(bp)))) 
    bp[which(duplicated(bp))] <- round(inverse_cdf(u))
    bp <- round(sort(c(na.omit(bp),runif(nSNPs/2-length(bp)-nQTN/2,min=0,max=max(rec_map$POS)))))
    bp <- sort(bp)
  }
  ## make sure no bp's are duplicaed and that the correlation between the observed density and the drawn density is at least 60%
  while(cor(density(bp,bw=100)$y,density(QTN_pos,bw=100)$y)^2<0.6 ){
    QTN_pos <- sort(sample(1:max(rec_map$POS),nQTN/2))
    
    
    d <- density(QTN_pos,bw=100,from = 0,to=max(rec_map$POS))  # KDE estimate
    dx <- diff(d$x)[1] 
    cdf_vals <- cumsum(d$y) * dx
    cdf_vals <- cdf_vals / max(cdf_vals)  
    
    inverse_cdf <- suppressWarnings(approxfun(cdf_vals, d$x, rule = 1)) # inverse CDF interpolator
    
    u <- runif(n_samples)            
    bp <- inverse_cdf(u)          
    bp <- as.vector(na.omit(bp))
    
    bp <- round(sort(c(na.omit(bp),runif(nSNPs/2-length(bp)-nQTN/2,min=0,max=max(rec_map$POS)))))
    
    while(any(duplicated(bp))){
      u <- runif(length(which(duplicated(bp)))) 
      bp[which(duplicated(bp))] <- round(inverse_cdf(u))
      bp <- round(sort(c(na.omit(bp),runif(nSNPs/2-length(bp)-nQTN/2,min=0,max=max(rec_map$POS)))))
      bp <- sort(bp)
    }
  }
  
  ## concatenate neutral and QTN SNPs and order them
  map_bp <- rbind(data.table(type="ntrl",bp),data.table(type="QTN",bp=QTN_pos))
  map_bp[,bp:=trunc(bp)]
  setorder(map_bp,bp)
  

  ## interpolate map
  ## this takes some time, probably a better way to do this
  cat("building map\n")
 
  
  cM <- unlist(mclapply(map_bp$bp,function(pos){
    ps <- rec_map[,which(POS>pos)][1]
    rec_map[ps,MEAN_POS]
  },mc.cores = 8)) 
  
  ## generate new map
  new_map <- data.table(Chr=1,map_bp,cM=cM)
  new_map[,rec_rate := c((new_map$cM[-1]-new_map$cM[-nrow(new_map)])/(new_map$bp[-1]-new_map$bp[-nrow(new_map)]),0)*10^6] ## estimate recombination rate
  
  new_map[,pos_nemo := trunc(cM*max(bp)/100)] ## generate positions for Nemo
  
  if(TRUE){
    
    setnames(rec_map,c("POS"),c("bp"))
    new_map[,map:="Extrapolated"]
    rec_map[,map:="Original"]
    dt <- rbind(new_map[,.(bp,cM,map)],rec_map[,.(bp,cM,map)])
    
    
    p1 <- ggplot(new_map, aes(bp,rec_rate)) +
      geom_line(col="steelblue") +
      theme_bw()+
      theme(aspect.ratio = 0.3)+
      ggtitle("Extrapolated map")
    
    p2 <- ggplot(rec_map,aes(bp,rec_rate)) +
      geom_line(col="salmon") +
      theme_bw()+
      theme(aspect.ratio = 0.3) +
      ggtitle("Original map")
    
    p3 <- ggplot(dt, aes(bp,cM,col=map)) +
      scale_color_manual(values=c("steelblue","salmon")) +
      geom_point(data=dt[map=="Extrapolated"]) +
      geom_line(data=dt[map=="Original"]) +
      theme_bw()+
      theme(aspect.ratio = 0.3,
            legend.position = "inside",
            legend.position.inside = c(0.1,0.7)) +
      ggtitle("Extrapolated vs. original")
    
    ## option to print to file
    jpeg(paste(outfolder,"recombination_map.jpeg",sep="/"),width = 10,height = 10,units = "in",res = 600)
    print(plot_grid(p1,p2,p3,ncol = 1))
    dev.off()
    
  }
  
  # ------------------------------- #
  #  Generate final maps 
  # ------------------------------- #
  
  message("writing files\n")
  ## duplicate map (second chromosomes has the same QTN but all with effect size=0)
  new_map2 <- copy(new_map)
  new_map2[,Chr:=2]
  ## merge them
  new_map <- rbind(new_map,new_map2)
  
  ## generate neutral map positions for nemo
  map_ntrl <- new_map[type=="ntrl"]
  
  ## write to file
  map_ntrl <- paste0("{",paste(sapply(split(map_ntrl,map_ntrl$Chr),function(chr){
    paste0("{",paste(chr$pos_nemo,collapse = ", "),"}")
  }),collapse = ""),"}")
  
  write.table(map_ntrl, file=paste(outfolder,"map_ntrl.txt",sep="/"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  ## QTN map positions for nemo
  map_QTN <- new_map[type=="QTN"]
  
  ## write to file
  map_QTN <- paste0("{",paste(sapply(split(map_QTN,map_QTN$Chr),function(chr){
    paste0("{",paste(chr$pos_nemo,collapse = ", "),"}")
  }),collapse = ""),"}")
  
  write.table(map_QTN, file=paste(outfolder,"map_QTN.txt",sep="/"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  
  # ------------------------------- #
  #  Allelic values
  # ------------------------------- #
  
  ## draw allelic values from  exponential distribution with mean=0.1 (rate=10). Second chromosome is neutral (all allelic values=0)
  ## make sure the starting phenotype is close to 0 (between -0.1 and 0.1) otherwise the simulations might not run.
  
  # draw allelic values
  allelic_values <- c(rexp(nQTN/2,rate = 10)+0.05,rep(0,nQTN/2))
  
  while(abs(sum(allelic_values*rep(c(1,-1),length(allelic_values)/2)))>0.1){
    allelic_values <- c(rexp(nQTN/2,rate = 10),rep(0,nQTN/2))      
  }
  

  ## write to file to be used by nemo
  write.table(paste0("{{",paste(allelic_values,collapse = ", "),"}}"), file=paste(outfolder,"allelic_values.txt",sep="/"), col.names = FALSE, row.names = FALSE, quote = FALSE)
  
  ## add allelic values to map we will save and use for downstream analyses
  new_map[type=="QTN",allelic_values:=allelic_values]
  
  ## save file for downstream analyses
  saveRDS(new_map,paste(outfolder,"rec_map.rds",sep="/"))

  
})
