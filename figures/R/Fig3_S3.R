#----------------------------------------------------------
#  Effect sizes (Fig. 3 and Fig. S3)
#----------------------------------------------------------
### === Fig S3 === ###
pop_stats <- readRDS("./data/pop_stats.rds")

eta_data <- do.call(rbind,lapply(c("adlt.fst","adlt.q1.Va","adlt.q1.Qst","Qst/Fst","allelic_values", "Va","a","b","n_QTN","cor_g_env","a_QTN/a_ntr","fitness.mean"),function(x){
  eta <- pop_stats[variable==x,eta_squared(aov(value~gene_flow*selection_intensity),alternative = "two.sided")]  
  eta$variable <- x
  eta
}))

eta_dt <- as.data.table(eta_data)
eta_dt[,Par:=c("sigma[S]^2","k","sigma[S]^2%*%k")[match(Parameter,c("selection_intensity","gene_flow","gene_flow:selection_intensity"))]]

eta_dt[,Var:=c("F[ST]","V[A]","Q[ST]","Q[ST]/F[ST]","bar(beta)","bar(n)[QTN]","r[g%~%e]^2","a[QTN]/a[ntrl]","bar(a)","b","bar(p)[V[A]*','*l]","bar(w)")[match(variable,c("adlt.fst","adlt.q1.Va","adlt.q1.Qst","Qst/Fst","allelic_values", "n_QTN","cor_g_env","a_QTN/a_ntr","a","b","Va","fitness.mean"))]]

eta_dt[, total_eta := sum(Eta2_partial), by = Var]
eta_dt[, Var := factor(Var, levels = unique(Var[order(-total_eta)]))]

pd <- position_dodge(width = 0.9, preserve = "single")

eta_dt[,tmp:=""]
p1 <- ggplot(eta_dt, aes(Var, Eta2_partial, fill = Par)) +
  geom_col(position = pd) +
  facet_grid(.~tmp)+
  geom_errorbar(
    aes(ymin = CI_low, ymax = CI_high, group = Par),
    position = pd,
    width = 0.2
  ) +
  scale_fill_manual(
    values = alpha(pal_jco()(10)[4:6], 0.75),
    name   = "Term",
    labels = function(x) parse(text = x)
  ) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_bw(base_size = 15) +
  theme(
    
    aspect.ratio = 1,
    axis.text.x  = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.8,0.8),
    plot.title = element_text(size = 13),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ylab(expression(italic(eta["partial"]^2))) 



pop_stats[,Var:=c("F[ST]","V[A]","Q[ST]","Q[ST]/F[ST]","bar(beta)","bar(n)[QTN]","r[g%~%e]^2","a[QTN]/a[ntrl]","bar(a)","b","bar(p)[V[A]*','*l]","bar(w)")[match(variable,c("adlt.fst","adlt.q1.Va","adlt.q1.Qst","Qst/Fst","allelic_values", "n_QTN","cor_g_env","a_QTN/a_ntr","a","b","Va","fitness.mean"))]]

pop_stats[, Var := factor(Var, levels = eta_dt[,levels(Var)])]

p2<- ggplot(pop_stats[!is.na(Var),.(value=mean(value)),by=.(Var,gene_flow,selection_intensity,rep,Chr_9sp,chr)], aes(selection_intensity,value,fill=gene_flow))+
  geom_violin(scale = "area",,draw_quantiles = c(0.5),adjust = 1.5) +
  facet_wrap(Var~.,scales="free",nrow=2,labeller = label_parsed) +
  scale_fill_manual(values = alpha(pal_jco()(10)[1:3],0.75),name="Gene flow",labels = function(x) parse(text = x))+
  scale_color_manual(values = alpha(pal_jco()(10)[1:3],0.75),name="Gene flow",labels = function(x) parse(text = x))+
  scale_x_discrete(labels = function(x) parse(text = x))+
  theme_bw()+
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle=90),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )
p2


pdf("./Figures/FigS3.pdf",width = 14,height = 5.5)
grid.arrange(p2+ggtitle("a) Effects of simulation parameters on population demography and local adaptation"),
             p1+ggtitle("b) Effect sizes"),nrow=1,widths=c(1,0.48))
dev.off()

## ---------------
## Fig S6
## ---------------

dt <- rbind(dt0[!(base_method %in% c("EMX","LFMM")) ,.(TP=mean(TP),FP=mean(FP),FN=mean(FN),PR=mean(PR),parameter="rho_w"),by=.(value=round(rho_w*20)/20,Method,AUC_meth,base_method)],
            dt0[!(base_method %in% c("EMX","LFMM")) ,.(TP=mean(TP),FP=mean(FP),FN=mean(FN),PR=mean(PR),parameter="alpha"),by=.(value=round(-log10(alpha)*2)/2,Method,AUC_meth,base_method)],
            dt0[!(base_method %in% c("EMX","LFMM")) ,.(TP=mean(TP),FP=mean(FP),FN=mean(FN),PR=mean(PR),parameter="l_min"),by=.(value=l_min,Method,AUC_meth,base_method)],
            dt0[,.(TP=mean(TP),FP=mean(FP),FN=mean(FN),PR=mean(PR),parameter="tau_C"),by=.(value=round(tau_C*20)/20,Method,AUC_meth,base_method)],
            dt0[,.(TP=mean(TP),FP=mean(FP),FN=mean(FN),PR=mean(PR),parameter="rho_d"),by=.(value=round(rho_d*20)/20,Method,AUC_meth,base_method)],
            dt0[,.(TP=mean(TP),FP=mean(FP),FN=mean(FN),PR=mean(PR),parameter="rho_ld"),by=.(value=round(rho_d*20)/20,Method,AUC_meth,base_method)],
            
            fill=TRUE,use.names=TRUE)
# 


par_map <- data.table(
  parameter = c("tau_C","rho_w","alpha","l_min","rho_d","rho_ld"),
  Par = c("tau[C]","w[rho]","alpha","l[min]","tau[d*','*rho]","tau[LD*','*rho]")
)

dt <- merge(dt, par_map, by="parameter", all.x=TRUE)


dt[,base_method:=factor(base_method,levels=c("LFMM",
                                             "LFMM´",
                                             "EMX",
                                             "EMX´",
                                             "Joint"))]

dt[,"n[OR]":=TP+FP]
dt2 <- melt(dt, measure.vars = c("n[OR]","TP","FP","FN","PR"),value.name = "N")

dt2[,Par:=factor(Par, levels=c("w[rho]","alpha","l[min]","tau[d*','*rho]","tau[LD*','*rho]","tau[C]"))]

dt2[AUC_meth=="PR^*",AUC_meth:="PR^'*'"]

p1 <- ggplot(dt2,aes(value,N,col = base_method,group=paste(Method,AUC_meth,base_method),linetype=AUC_meth))+
  geom_line(linewidth=0.5) +
  facet_wrap(variable~Par,scales="free",labeller = label_parsed,ncol=6)+
  scale_linetype_manual(values = c(2,1),name=NULL,labels = function(x) parse(text = x))+
  scale_color_manual(values = pal_jco()(10),name=NULL)+#,labels = function(x) parse(text = x))+
  labs(linetype=NULL)+
  guides(fill = guide_legend(nrow = 1))+
  theme_bw(base_size = 15)+
  scale_x_continuous(labels = number_format(accuracy = 0.1))+
  theme(
    aspect.ratio = 0.75,
    axis.text.x = element_text(angle = 90, hjust = 1),
    #strip.text.x = element_text(size = 12),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank())

 
### === Fig 3 === ###
dt <- pop_stats[Var %in% c("Q[ST]/F[ST]","b","a[QTN]/a[ntrl]"),
                .(value = mean(value)), 
                by = .(Var,gene_flow,selection_intensity,rep,Chr_9sp,chr)]

dt[Var=="b",Var:="'Bakcground LD'"]

dt[,Var := factor(Var, levels = c("Q[ST]/F[ST]",
                                  "'Bakcground LD'",
                                  "a[QTN]/a[ntrl]"))]

dt_hline <- data.table(
  Var = c("Q[ST]/F[ST]", "'Bakcground LD'", "a[QTN]/a[ntrl]"),
  yintercept = c(1, NA, 1)
)

dt_hline[, Var := factor(Var,levels = c("Q[ST]/F[ST]",
                                        "'Bakcground LD'",
                                        "a[QTN]/a[ntrl]"))]



p_sim <- ggplot(dt,aes(selection_intensity, value, fill = gene_flow)) +
  geom_violin(scale = "width", draw_quantiles = c(0.5), adjust = 1) +
  facet_wrap(Var ~ ., scales = "free", nrow = 1, labeller = label_parsed) +
  scale_fill_manual(values = alpha(pal_jco()(10)[1:3],0.75),
                    name="Gene flow",
                    labels=function(x) parse(text=x)) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  geom_hline(
    data = dt_hline,
    aes(yintercept = yintercept),
    linetype=2, colour="black", linewidth=0.5
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank(),
    #legend.position = "bottom",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )
pdf("./Figures/Fig3.pdf",height = 3,width = 9)
p_sim
dev.off()

