
#----------------------------------------------------------
#  Supplementary Figure S4
#----------------------------------------------------------

if(length(list.files("./data/"))==0) message("Please download data form Zenond")
## read data generated in ./R_sim/outlier_analyses_sim.R
top_10 <- readRDS("./data/top_10.rds")
data_AUC <- readRDS("./data/data_AUC.rds")
## false positives ##

dt <- melt(top_10, measure.vars = c("FN","TP","PR","FP_ntrl","precision","recall"))
dt <- dt[,.(value=mean(value)),.(sim, gene_flow, selection_intensity,base_method,AUC_meth, rep,variable)]

dt[variable=="FP_ntrl" ,variable:="FP[ntrl]"]
dt[,base_method:=factor(base_method,levels=c("EMX","LFMM","EMX´","LFMM´","Joint"))]

dt <- rbind(dt,melt(data_AUC, measure.vars = c("AUC_norm"))[,.(sim, gene_flow, selection_intensity,base_method, AUC_meth,    rep, variable,value)],fill=TRUE,use.names=TRUE)

dt[variable=="AUC_norm",variable:="AUC"]
dt[,variable:=factor(variable,levels=c("TP","FP[QTN]","FP[ntrl]","FN","precision","recall","PR","AUC"))]


p_sep <- ggplot(dt,aes(base_method,value,col=AUC_meth))+
  geom_boxplot(outlier.size = 0.5) +
  scale_color_manual(values = alpha(pal_jco()(10)[c(4,1)],1),labels = function(x) parse(text = x),name="AUC-")+
  facet_grid(variable~sim,labeller = labeller(variable=label_parsed),scales="free_y") +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()) 

p_sep

pdf("./figures/FigS4.pdf",width = 12,height = 10)
p_sep
dev.off()

