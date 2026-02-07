library(data.table)
library(parallel)
library(egg)
library(ggridges)
library(igraph)
library(purrr)
library(dplyr)
library(SNPRelate)
library(stats)
library(ggplot2)
library(SNPRelate)
library(data.table)
library(ggplot2)
invisible(lapply(c("./R/emmax.R",
                   "./R/scale_F_with_ldw.R",
                   "./R/coef_ld_dec.R",
                   "./R/emmax.R",
                   "./R/get_bg_ld.R",
                   "./R/get_C.R",
                   "./R/find_ORs.R",
                   "./R/get_el.R",
                   "./R/get_ld_w_draws.R",
                   "./R/ld_decay.R",
                   "./R/ORs_from_draws.R",
                   "./R/prep_manhattan.R",
                   "./R/SLC.R",
                   "./R/d_from_rho.R"),source))

###################################################
#### LD-decay, LD-scaling and C-scores for 9sp ####
###################################################
#### ==== prepare data and do initial EMMAX and LFMM analyses ==== ####

if(!file.exists("./empirical_data")) message("Please download data form Zenondo")

data_9sp <-readRDS("./empirical_data/9sp/9sp_data.rds") ## contains SNP_res_9sp, GTs_9sp and pheno_9sp
GTs_9sp <- data_9sp$GT
map_9sp <- data_9sp$map
pheno_9sp <- data_9sp$pheno
eco_9sp <- data_9sp$ecotype_bin
##filter by maf
keep <- map_9sp$maf>0.1
GTs_9sp <- GTs_9sp[,map_9sp$maf>0.1]
map_9sp <- map_9sp[maf>0.1]
rm(data_9sp)
gc()
if(any(ls() %in% "gds_9sp")) snpgdsClose(gds_9sp); unlink(gds_9sp)
gds_9sp <- create_gds_from_geno(geno = GTs_9sp, map=map_9sp,"gds_9sp.gds")

# get id's for later use
ids <- read_gds_ids(gds_9sp)



## Analysis pipeline for nine spined sticklebacks.
if(!file.exists("./empirical_data/9sp/SNP_res_9sp.rds")){
  
  #### ==== Estimate LD-decay for chromosomes ==== ####
  if(!file.exists("./empirical_data/9sp/LD_decay_9sp.rds")){
    b <- get_bg_ld(gds_9sp, n_sub = 5000, q = 0.95)
    
    LD_decay_9sp <- ld_decay(gds=gds_9sp,
                             q = 0.95,
                             b = b, 
                             n_sub = 5000,
                             slide_win_ld = 1000,
                             window_size = 1e6,
                             step_size = 5e5,
                             dist_unit = 5000,
                             n_cores_ld =  8)
    
    LD_decay_9sp$summary[,plot(Chr_size,1/a)]
    saveRDS(LD_decay_9sp,"./empirical_data/9sp/LD_decay_9sp.rds")
    
  }else{
    LD_decay_9sp <- readRDS("./empirical_data/9sp/LD_decay_9sp.rds")  
  }
  
  
  
  #### ==== Get LD-pruned data using single linkage clustering with threshold based on 75% quantile of LD-decay rate ==== ####
  if(!file.exists("./empirical_data/9sp/pruning.rds")){
    
    
    cls_05_9sp <- SLC(gds_9sp,decay_tbl = LD_decay_9sp,idx = NULL,slide_win_ld = 1000,q = 0.5)$CLS ## based on both ld and d thresholds, as for ORs, but singleton clusters are kept
    pruned_SNPs <- sapply(cls_05_9sp,sample,1) ## from clusters, pick one randomly
    
    saveRDS(list(cls=cls_05_9sp,pruned_SNPs=pruned_SNPs),"./empirical_data/9sp/pruning.rds")
    
  }else{
    pruned_SNPs <- readRDS("./empirical_data/9sp/pruning.rds")$pruned_SNPs
  }
  
  #### ==== EMMAX analyses based on LD-pruned data ==== ####
  
  
  if(!file.exists("./empirical_data/9sp/emx_9sp.rds")){
    ## GRM based on pruned SNPs
    GRM <- snpgdsGRM(gds_9sp,method = "GCTA",snp.id = pruned_SNPs,verbose = FALSE,autosome.only = FALSE)$grm
    
    ## the binary phenotype
    eco_bin  <- as.numeric(as.factor(data_9sp$ecotype_bin))
    
    
    ## EMMAX does not expect a file in 012 format so the maximum likelihood genotypes can be used (without rounding)
    emx <- emmax(eco_9sp,GTs_9sp,K = GRM,Covar =pheno_9sp$lineage) ## function lives in ./R/emmax.R
    
    map_9sp[,emx_F:=emx$F] ## add to map
    emx_gif = map_9sp[,median(emx_F)/qf(0.5,1,nrow(GTs_9sp),lower.tail = FALSE)] ## inflation factor
    
    #map_9sp[,emx_F_GC:=emx_F/emx_gif]  ## genomic control
    #emx_gif<1; no genomic control 
    
    map_9sp[,emx_p:=pf(emx_F,1,nrow(GTs_9sp),lower.tail = FALSE)] ## p-value
    map_9sp[,emx_q:=p.adjust(emx_p,"fdr")] ## fdr correction
    saveRDS(emx,"./empirical_data/9sp/emx_9sp.rds")
    
  }else{
    map_9sp[,emx_F:=readRDS("./empirical_data/9sp/emx_9sp.rds")$F] ## add to map
    map_9sp[,emx_p:=pf(emx_F,1,nrow(GTs_9sp),lower.tail = FALSE)] ## p-value
    map_9sp[,emx_q:=p.adjust(emx_p,"fdr")] ## fdr correction
  }
  
  
  #### ==== Latent factor mixed model analyses (LFMM) ==== #### 
  #lfmm takes very long time
  if(!file.exists("./empirical_data/9sp/lfmm_F.rds")){
    library(LEA)
    
    eco_bin  <- as.numeric(as.factor(data_9sp$ecotype_bin))
    write.lfmm(GTs_9sp, "./tmp/genotypes.lfmm")
    write.env(eco_bin, "./tmp/gradients.env")
    project = NULL
    
    project = lfmm2("./tmp/genotypes.lfmm", "./tmp/gradients.env", K=data_9sp$pheno[,length(unique(pop_locality))]) # K is number of populations
    pv = lfmm2.test(project, "./tmp/genotypes.lfmm", "./tmp/gradients.env",genomic.control = TRUE,full = TRUE)
    
    saveRDS(pv$f,"./empirical_data/9sp/lfmm_F.rds") # only F-value is needed
    q("no")
  }else{
    map_9sp[,lfmm_F:=readRDS("./empirical_data/9sp/lfmm_F.rds")[keep]]
    lfmm_gif = map_9sp[,median(lfmm_F)/qf(0.5,1,nrow(GTs_9sp),lower.tail = FALSE)] ## inflation factor
    map_9sp[,lfmm_F_GC:=lfmm_F/lfmm_gif]  ## genomic control
    map_9sp[,lfmm_P:=pf(lfmm_F_GC,1,nrow(GTs_9sp),lower.tail = FALSE)]
    map_9sp[,lfmm_q:=p.adjust(lfmm_P,"fdr")]
  }
  
  #### ==== Draw 100 rho_w values and get ld_w vector for each of them ==== ####
  
  if(!file.exists("./empirical_data/9sp/ld_w_draws_9sp_200.rds")){
    ld_w_draws_9sp <- get_ld_w_draws(gds_9sp,decay_tbl = LD_decay_9sp,slide_win_ld = 1000,n_draws = 200,rho_min = 0.75,rho_max = 1,n_cores = 8)
    saveRDS(ld_w_draws_9sp,"./empirical_data/9sp/ld_w_draws_9sp_200.rds")
    cat("Done")
    
  }
  
  #### ==== For each of the 100 draws, draw 25 values for alpha, rho_OR and l_min and get C-scores ==== ####
  if(!file.exists("./empirical_data/9sp/OR_draws_9sp_lmin1_45_2.rds")){
    
    df2 = nrow(GTs_9sp)-2
    F_vals <- map_9sp[,.(emx=emx_F,lfmm=lfmm_F_GC)]
    
    q_orgs <- apply(F_vals,2,function(Fval){
      p <- pf(Fval, df1=1, df2, lower.tail = FALSE)
      q <- p.adjust(p,"fdr")
    })
    
    #ld_w_draws_9sp <- readRDS("./empirical_data/9sp/ld_w_draws_9sp_200.rds")
    OR_draws_9sp <- cbind(sub=i,get_ORs_for_draw(gds=gds_9sp,
                                                 F_vals = F_vals,
                                                 q_orgs = q_orgs,
                                                 decay_tbl = LD_decay_9sp,
                                                 ld_w_draws = ld_w_draws_9sp,
                                                 n_inds=nrow(GTs_9sp),
                                                 n_draws = 25,
                                                 rho_OR_lim = NULL,
                                                 rho_ld_lim = list(min=0.5, max=1.0),
                                                 rho_d_lim = list(min=0.90, max=1.0),
                                                 alpha_lim = list(min=-log10(0.05), max=4),
                                                 lmin_lim = list(min=1, max=45),
                                                 n_cores = 1))
    
    
    
    saveRDS(OR_draws_9sp,"./empirical_data/9sp/OR_draws_9sp_lmin1_45.rds")
    rm(ld_w_draws_9sp)
    gc()
  }
  saveRDS(map_9sp,"./empirical_data/9sp/SNP_res_9sp.rds")

}

## --------------------------
## Manhattan for 9sp stickleback data
## --------------------------

## data created above and/or available from Zenondo
OR_draws_9sp <- readRDS("./empirical_data/9sp/OR_draws_9sp_lmin1_45.rds")
SNP_res_9sp <- readRDS("./empirical_data/9sp/SNP_res_9sp.rds")
SNP_res_9sp <- cbind(SNP_res_9sp,get_C(OR_draws_9sp,markers=SNP_res_9sp$marker))
rm(OR_draws_9sp)
gc()
LD_decay_9sp <- readRDS("./empirical_data/9sp/LD_decay_9sp.rds")

data_9sp <-readRDS("./empirical_data/9sp/9sp_data.rds") ## contains SNP_res_9sp, GTs_9sp and pheno_9sp
GTs_9sp <- data_9sp$GT
map_9sp <- data_9sp$map
keep <- map_9sp$maf>0.1
GTs_9sp <- GTs_9sp[,map_9sp$maf>0.1]
map_9sp <- map_9sp[maf>0.1]
rm(data_9sp)
gc()
gds_9sp <- create_gds_from_geno(geno = GTs_9sp, map=map_9sp,"gds_9sp.gds")

## specify tau_C
tau_C <- 0.2

## get ORs with rho_ld=0.999 and rho_d=0.999
ORs <- find_ORs(gds_9sp,LD_decay_9sp,outliers = SNP_res_9sp[C_joint>tau_C,marker],rho_ld=0.999,rho_d=0.999)

## add ORs to data
SNP_res_9sp[,OR:=ORs$OR[match(marker, ORs$marker)]]
SNP_res_9sp[,n_loci:=ORs$n_loci[match(marker, ORs$marker)]]

plot_data_manh <- prep_manhattan(SNP_res_9sp[,.(bp=Pos,Chr,marker,C_joint, OR = as.character(OR),     # <- Joint_C ORs
                                                LFMM_log = -log10(lfmm_q),
                                                lfmm_q)],chr_cols = c("white","grey80"),spacer =0)



## manhattan plot for lfmm
p1 <- ggplot(plot_data_manh$data, aes(BPcum, LFMM_log, col = OR)) +
  geom_rect(
    data = plot_data_manh$rect,
    aes(xmin = x1, xmax = x2, ymin = y1, ymax = Inf),
    fill = plot_data_manh$rect$col,
    alpha = 0.5,
    linewidth = 0.25,
    inherit.aes = FALSE
  ) +
  geom_point(size = 1) +
  geom_point(data = plot_data_manh$data[!is.na(OR)],size = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = rep(col_vector, 4), guide = "none") +
  scale_x_continuous(
    label = plot_data_manh$axis$Chr,
    breaks = plot_data_manh$axis$center
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 0.25,
    axis.text.x = element_text(angle = 90),
    panel.grid = element_blank()
  ) +
  xlab(NULL) +
  ylab(expression(-log[10](q)))

## manhattan plot for Joint-C
p2 <- ggplot(plot_data_manh$data, aes(BPcum, C_joint, col = OR)) +
  geom_rect(
    data = plot_data_manh$rect,
    aes(xmin = x1, xmax = x2, ymin = y1, ymax = Inf),
    fill = plot_data_manh$rect$col,
    alpha = 0.5,
    linewidth = 0.25,
    inherit.aes = FALSE
  ) +
  geom_point(size = 1) +
  geom_point(data = plot_data_manh$data[!is.na(OR)],size = 1) +
  geom_hline(yintercept = tau_C, linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = rep(col_vector, 4), guide = "none") +
  scale_x_continuous(
    label = plot_data_manh$axis$Chr,
    breaks = plot_data_manh$axis$center
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 0.25,
    axis.text.x = element_text(angle = 90),
    panel.grid = element_blank()
  ) +
  xlab(NULL) +
  ylab(expression(Joint[C]))


p_9sp <- grid.arrange(p1+ggtitle(expression("g) Nine-spined sticklebacks | LFMM" )),
                      p2+ggtitle(expression("h) Nine-spined sticklebacks | Joint" ))
                      ,ncol=1)

saveRDS(p_9sp,"./figures/p_9sp.rds")

## stats 
nrow(SNP_res_9sp)
SNP_res_9sp[!is.na(OR),length(unique(OR))]
SNP_res_9sp[!is.na(OR),length(unique(OR)),by=Chr]

all_chr <- SNP_res_9sp[, .(Chr = unique(Chr))]

# count ORs per chromosome
cnt <- SNP_res_9sp[!is.na(OR),
                   .(n_OR = uniqueN(OR)),
                   by = Chr]

# merge and fill zeros
res <- merge(all_chr, cnt, by = "Chr", all.x = TRUE)
res[is.na(n_OR), n_OR := 0]
res[,mean(n_OR)]
res[,table(n_OR==0)]
res[,table(n_OR>1)]

33/17

26/2.9

17/9
33/26
