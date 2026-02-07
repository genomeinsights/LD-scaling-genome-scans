library(data.table)
library(egg)
library(dplyr)
library(SNPRelate)
library(stats)
library(ggplot2)


if(!file.exists("./results_sim_new/")) message("Please downlad files form Zenondo")

invisible(lapply(c("./R/get_ld_w_draws.R",         ## Generates w for ld_w for drawn w-values
                   "./R/prep_manhattan.R"),         ## Function for manhattan plots
                   source))


p_Va_th <- 0.05

res_file <- "./results_sim/c1_V0.5_rep4.rds"

## extract data
data <- readRDS(res_file)
PR_data <- readRDS(gsub(".rds","_PR.rds",res_file))

## genotype data set used to set up gds object below
GTs <- data$GTs

## this is where all data will be added
SNP_res <- data$SNP_res

## decay information
decay_tbl <- data$decay_tbl

## performance data
PR_data$PR_data[method=="lfmm_F"][which.max(PR),paste(paste0("TP=",TP),paste0("FP=",FP),paste0("FN=",FN),paste0("PR=",round(PR,2)),sep=" | ")]
PR_data$PR_data_C[method=="C_joint"][which.max(PR),paste(paste0("TP=",TP),paste0("FP=",FP),paste0("FN=",FN),paste0("PR=",round(PR,2)),sep=" | ")]

## threshols for best draws
tau_C <- PR_data$PR_data_C[method=="C_joint"][which.max(PR),tau_C]
alpha_lfmm <- PR_data$PR_data[method=="lfmm_F"][which.max(PR),alpha]


## get true and false posities for joint analyses
ORs_max_PR <- PR_data$PR_data_C[method=="C_joint"][which.max(PR),ORs][[1]]
ORs_max_PR <- data.table(unlist(ORs_max_PR),rep(names(ORs_max_PR),sapply(ORs_max_PR,length)))
SNP_res[,OR_joint := ORs_max_PR$V2[match(marker, ORs_max_PR$V1)]]

## get true and false posities for LFMM
ORs_max_PR <- PR_data$PR_data[method=="lfmm_F"][which.max(PR),ORs][[1]]
ORs_max_PR <- data.table(unlist(ORs_max_PR),rep(names(ORs_max_PR),sapply(ORs_max_PR,length)))
SNP_res[,OR_lfmm := ORs_max_PR$V2[match(marker, ORs_max_PR$V1)]]


## create gds object
# close if necessary before
#snpgdsClose(gds); unlink(gds_path)
gds_path = tempfile(fileext = ".gds")

gds <- create_gds_from_geno(geno = GTs, map=SNP_res,gds_path)
ids <- read_gds_ids(gds)

## get ld_w for rho_w=0.9 and add to SNP_res
ld_w_draws <- get_ld_w_draws(gds,decay_tbl = decay_tbl,slide_win_ld = 1000,n_draws = 1,rho_min = 0.9,rho_max = 0.91)
ld_w <- rbindlist(ld_w_draws[draw==1,ld_w])
SNP_res[,ld_w := ld_w[match(ids$snp_id,SNP),median_r2]]

SNP_res[,cor(max_LD_with_QTN,C_joint)^2]
SNP_res[,cor(max_LD_with_QTN,lfmm_q)^2]
SNP_res[,cor(max_LD_with_QTN,ld_w,use = "pair")^2]


## prepare data for manhattan plot
plot_data_manh <- prep_manhattan(SNP_res[,.(bp=Pos,Chr,marker,C=C_joint,LFMM=-log10(lfmm_q),OR_joint,OR_lfmm,type,ld_w,max_LD_with_QTN)],chr_cols = c("white","grey80"),spacer =0)


p1 <-  ggplot(plot_data_manh$data, aes(BPcum,C,col=OR_joint,shape=type)) +
  geom_rect(data=plot_data_manh$rect, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=Inf),fill=plot_data_manh$rect$col,alpha=0.5, linewidth = 0.25,inherit.aes = FALSE) +
  geom_point(data=plot_data_manh$data,size=1,alpha=1) +
  geom_point(data=plot_data_manh$data[type=="QTN"],size=3,alpha=1,shape=3) +
  scale_shape_manual(values = c(20,3),name=NULL)+
  scale_color_manual(values =rep(col_vector,3),guide="none")+
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  geom_hline(yintercept = tau_C,linetype=2,linewidth=0.5)+
  scale_x_continuous(label = plot_data_manh$axis$Chr, breaks= plot_data_manh$axis$center ) +
  theme(
    aspect.ratio = 0.25,
    axis.text.x = element_text(angle=90),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -7),
    strip.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  xlab(NULL)+
  ylab(expression(Joint[C]))
p1

p2 <-  ggplot(plot_data_manh$data, aes(BPcum,LFMM,col=OR_joint,shape=type)) +
  geom_rect(data=plot_data_manh$rect, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=Inf),fill=plot_data_manh$rect$col,alpha=0.5, linewidth = 0.25,inherit.aes = FALSE) +
  geom_point(data=plot_data_manh$data[],size=1,alpha=1) +
  geom_point(data=plot_data_manh$data[LFMM<1.32 & type=="QTN"],size=1,alpha=1,aes(col=OR_joint)) +
  geom_point(data=plot_data_manh$data[type=="QTN" & LFMM>1.32],size=3,alpha=1,shape=3,aes(col=OR_joint)) +
  scale_shape_manual(values = c(20,3),name=NULL)+
  scale_color_manual(values =rep(col_vector,3),guide="none")+
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = 3)))+
  geom_hline(yintercept =1.31,linetype=2,linewidth=0.5)+
  scale_x_continuous(label = plot_data_manh$axis$Chr, breaks= plot_data_manh$axis$center ) +
  theme(
    aspect.ratio = 0.25,
    axis.text.x = element_text(angle=90),
    strip.background = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    legend.box.margin = margin(t = 0, r = 0, b = 0, l = -7),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  xlab(NULL) +
  ylab(expression(-log10(q)))
p2

p3 <-  ggplot(plot_data_manh$data, aes(BPcum,C,col=max_LD_with_QTN)) +
  geom_rect(data=plot_data_manh$rect, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=Inf),fill=plot_data_manh$rect$col,alpha=0.5, linewidth = 0.25,inherit.aes = FALSE) +
  geom_point(data=plot_data_manh$data,size=1,alpha=1) +
  geom_point(data=plot_data_manh$data[type=="QTN" & !is.na(OR_joint)],size=3,alpha=1,shape=3) +
  scale_shape_manual(values = c(2,20,3))+
  scale_color_gradientn(colours = pal,name=expression(r[f]^2),labels = number_format(accuracy = 0.1)) +
  theme_bw() +
  geom_hline(yintercept = tau_C,linetype=2,linewidth=0.5)+
  scale_x_continuous(label = plot_data_manh$axis$Chr, breaks= plot_data_manh$axis$center ) +
  theme(
    aspect.ratio = 0.25,
    legend.position.inside =  c(0.3,0.9),
    axis.text.x = element_text(angle=90),
    strip.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  xlab(NULL)+
  ylab(expression(Joint[C]))



p4 <-  ggplot(plot_data_manh$data, aes(BPcum,ld_w,col=max_LD_with_QTN)) +
  geom_rect(data=plot_data_manh$rect, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=Inf),fill=plot_data_manh$rect$col,alpha=0.5, linewidth = 0.25,inherit.aes = FALSE) +
  geom_point(data=plot_data_manh$data,size=1,alpha=0.5) +
  geom_point(data=plot_data_manh$data[type=="QTN" & !is.na(OR_joint)],size=3,alpha=1,shape=3) +
  scale_shape_manual(values = c(2,20,3))+
  scale_color_gradientn(colours = pal,name=expression(r[f]^2)) +
  theme_bw() +
  scale_x_continuous(label = plot_data_manh$axis$Chr, breaks= plot_data_manh$axis$center ) +
  theme(
    aspect.ratio = 0.25,
    legend.position.inside =  c(0.3,0.9),
    axis.text.x = element_text(angle=90),
    strip.background = element_blank(),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  xlab(NULL)+
  ylab(expression(ld[w]))


## get AUC data
data_AUC <- readRDS("./data/data_AUC.rds")
AUC_lfmm <- data_AUC[c==SNP_res[1,c] & V==SNP_res[1,V] & rep==SNP_res[1,rep]  & method=="lfmm_q",round(AUC_norm,2)]
AUC_Joint <- data_AUC[c==SNP_res[1,c] & V==SNP_res[1,V] & rep==SNP_res[1,rep]  & method=="C_joint",round(AUC_norm,2)]
AUC_lfmm
AUC_Joint

## load manhattanf figures for empirical data
p_9sp <- readRDS("./figures/p_9sp.rds")
p_3sp <- readRDS("./figures/p_3sp.rds")
p_sticklebakcs <-grid.arrange(p_3sp,p_9sp,ncol=1)

p_mahh_sim <- grid.arrange(p2+ggtitle(expression("a) Simulation example | LFMM" )),
                           p1+ggtitle(expression("b) Simulation example | Joint" )),
                           p3+ggtitle(expression("c) Simulation example | Joint" )),
                           p4+ggtitle(expression("d) Simulation example | "*ld["w"*","*rho*"="*0.9])),ncol=1)


jpeg("./figures/Fig5.jpeg",height = 8,width = 12,units = "in",res = 300)
grid.arrange(p_mahh_sim,p_sticklebakcs,nrow=1)
dev.off()

