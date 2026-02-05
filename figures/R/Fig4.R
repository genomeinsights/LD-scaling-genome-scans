library(data.table)
library(effectsize)
library(egg)
library(dplyr)
library(stats)
library(ggplot2)

if(!file.exists("./data/")) message("Please downlad files form Zenondo")
## from /R_sim/outlier_analyses_sim.R
top_10 <- readRDS("./data/top_10.rds")

dt2 <- melt(top_10[Method %in% c("'EMX'","'LFMM'","'LFMM´'","'EMX´'","'EMX-C'","'LFMM-C'","'LFMM´-C'","'EMX´-C'","'Joint-C'","'Joint'")],
            measure.vars = c("precision","recall","PR"))

dt2[variable=="alpha",value:=-log10(value)]
dt2[variable=="alpha",variable:="-log[10]*'('*alpha*')'"]
dt2[variable=="rho_w",variable:="w"]
dt2[variable=="rho_ld",variable:="tau[LD]"]
dt2[variable=="rho_d",variable:="tau[d]"]
dt2[variable=="l_min",variable:="l[min]"]
dt2[variable=="C_th",variable:="C"]


dt3 <- dt2[,.(median=median(value),sd=sd(value)),by=.(Method,variable)]

# format values
df <- dt3 %>%
  mutate(value = sprintf("%.2f (%.2f)", median, sd)) %>%
  select(Method,variable, value)

# pivot to wide format
df_wide <- df %>%
  pivot_wider(names_from = variable, values_from = value)

# print LaTeX table
kable(df_wide,
      format = "latex",
      booktabs = TRUE,
      caption = "Summary of median (SD) values for each method and variable.")

#### ==== variance explained by parameters ==== ####
library(effectsize)

eta_list <- list()
params <- c("alpha", "l_min", "rho_w")
for (p in params) {
  print(p)
  f <- as.formula(paste(p, "~ gene_flow + selection_intensity + Method"))
  model <- aov(f, data = top_10[!is.na(Method) & !grepl("C",Method)])
  eta <- eta_squared(model, partial = TRUE,alternative = "two.sided")
  eta$parameter <- p
  eta_list[[p]] <- eta
}

params <- c("rho_ld", "rho_d")
for (p in params) {
  print(p)
  f <- as.formula(paste(p, "~ gene_flow + selection_intensity + Method"))
  model <- aov(f, data = top_10[!is.na(Method)])
  eta <- eta_squared(model, partial = TRUE,alternative = "two.sided")
  eta$parameter <- p
  eta_list[[p]] <- eta
}

params <- c("tau_C")
for (p in params) {
  print(p)
  f <- as.formula(paste(p, "~ gene_flow + selection_intensity + Method"))
  model <- aov(f, data = top_10[!is.na(Method) & grepl("C",Method)])
  eta <- eta_squared(model, partial = TRUE,alternative = "two.sided")
  eta$parameter <- p
  eta_list[[p]] <- eta
}

eta_summary <- rbindlist(eta_list)
## supplementary table ##
knitr::kable(eta_summary, format = "latex", booktabs = TRUE, digits = 3)

# add a new group column
eta_summary$group <- "other"
eta_summary$group[eta_summary$Parameter == "Method"] <- "Method"
eta_summary$group[eta_summary$Parameter %in% c("gene_flow", "Selection_intensity",
                                               "gene_flow:Selection_intensity")] <- "Simulation parameters"

# 3) summarize total eta² by group for each parameter
unique_params <- unique(eta_summary$parameter)
groups <- c("Method", "Simulation parameters")

eta_group_summary <- data.frame(parameter = character(),
                                group = character(),
                                total_eta = numeric(),
                                stringsAsFactors = FALSE)


for (p in unique_params) {
  for (g in groups) {
    vals <- eta_summary$Eta2_partial[eta_summary$parameter == p & eta_summary$group == g]
    eta_group_summary <- rbind(eta_group_summary,
                               data.frame(parameter = p,
                                          group = g,
                                          total_eta = sum(vals)))
  }
}

param_labels <- c(
  alpha      = "alpha",
  rho_ld = "tau[LD]",
  rho_d  = "tau[d]",
  l_min      = "l[min]",
  rho_w      = "w",
  tau_C      = "tau[C]"
)

eta_group_summary$group <- factor(eta_group_summary$group,
                                  levels = c("Simulation parameters", "Method"))
eta_group_summary$parameter <- factor(eta_group_summary[,"parameter"],levels=c("tau_C","rho_ld","rho_d","alpha","l_min","rho_w"))
eta_group_summary <- data.table(eta_group_summary)

levs <- eta_group_summary[,sum(total_eta),by=parameter][order(-V1)]$parameter
eta_group_summary[,parameter:=factor(parameter,levels=levs)]


p_PVE <- ggplot(eta_group_summary, aes(x = parameter, y = total_eta, fill = group)) +
  geom_col(position = "stack", color = "black", width = 0.7,linewidth=0.5) +
  scale_fill_manual(values = c("Simulation parameters" = pal_jco()(2)[1],
                               "Method" = pal_jco()(2)[2], name = ""),name=NULL) +
  scale_x_discrete(labels = function(x) parse(text = param_labels[x])) +
  labs(x =" ",
       y = expression("Total partial "*eta^2)) +
  theme_bw(base_size = 20) +
  theme(aspect.ratio = 0.5,
        legend.position = "inside",
        legend.position.inside =  c(0.7,0.8),
        legend.background = element_blank(),
        axis.text.x=element_text(size=20)
  )

### --- Detection power --- ###

  
res_files <- list.files( "./results_sim/",full.names = TRUE) ### available from Zenond
PR_files <- res_files[grepl("_PR.rds",res_files)]
all_qtn <- readRDS("./data/all_qtn.rds")

QTN_map <- list()
for (PR_file in PR_files) {
  cat("Working on", PR_file, "\n")
  
  base_name <- sub("_PR\\.rds$", "", basename(PR_file))
  
  data_path <- res_files[grep(paste0(base_name, "."), res_files, fixed = TRUE)][1]
  data <- readRDS(data_path)
  QTN_map[[PR_file]] <- data$SNP_res[type=="QTN",.(c,V,rep,marker,p_Va)]
}
qtn_map <- rbindlist(QTN_map)

keys <- c("c","V","rep")
brks <- c(0, 0.25, 0.5, 0.75, 1)

# all_qtn: one row per dataset, with list-column p_Va for *all* segregating QTNs
avail_long <- all_qtn[, .(p_Va = unlist(p_Va)), by = keys]
avail_long[, bin := cut(p_Va, breaks = brks, include.lowest = TRUE)]
avail_counts <- avail_long[, .(N_avail = .N), by = c(keys, "bin")]

# Keep only draws with at least one recovered focal QTN
rec_long <- all_draws[AUC_meth=="PR[C]"&
                        !is.na(OR_focals) & lengths(OR_focals) > 0,
                      .(focal_QTN = unlist(OR_focals)),          # expand recovered QTNs
                      by = c(keys,  "base_method")
]



rec_long <- merge(rec_long, qtn_map, by.x = c(keys, "focal_QTN"), by.y = c(keys, "marker"), all.x = TRUE)

rec_long[, bin := cut(p_Va, breaks = brks, include.lowest = TRUE)]

# Count unique QTNs recovered in each bin, per draw
rec_counts <- rec_long[, .(N_rec = uniqueN(focal_QTN)), by = c(keys, "base_method", "bin")]

power_draw <- merge(rec_counts, avail_counts, by = c(keys, "bin"), all.x = TRUE)

power_draw[, power := N_rec / N_avail]


# Summarize across draws first (within dataset), then across datasets
power_ds <- power_draw[, .(power = mean(power, na.rm = TRUE)), by = c(keys,  "base_method", "bin")]
power_sum <- power_ds[, .(
  power_mean = mean(power, na.rm = TRUE),
  power_sd   = sd(power, na.rm = TRUE)
), by = .(base_method, bin)]

power_sum[,base_method:=factor(base_method,levels=c("'LFMM'",
                                                    "'LFMM´'",
                                                    "'EMX'",
                                                    "'EMX´'",
                                                    "'Joint'"))]

p_recovered <- ggplot(power_sum, aes(bin, power_mean, group = base_method,color = base_method)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  ylab("Pr(QTN recovered in a draw | QTN in bin)") +
  xlab("Binned p(V_A) of focal QTN") +
  scale_color_manual(values = alpha(pal_jco()(10),1),labels = function(x) parse(text = x))+
  theme_bw(20) +
  theme(aspect.ratio = 0.5,
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.text = element_text(size = 20),
        legend.position.inside = c(0.13,0.6),
        legend.background = element_blank(),
  ) +
  ylab(expression("QTN recovered")) +
  xlab(expression(V[A]* " explained by QTN")) 

#### ---- robustness ---- ####

p_R <- ggplot(data_AUC[R!=0], aes(x = base_method, y = R, fill = AUC_meth)) +
  geom_boxplot(width = 1, color = "black",outliers = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = alpha(pal_jco()(10),1),name="Method",labels = function(x) parse(text = x))+
  labs(x =" ", y = expression(italic(R))) +
  theme_bw(base_size = 20) +
  theme(aspect.ratio = 0.5,
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.text = element_text(size = 20),
        legend.position.inside = c(0.1,0.18),
        legend.background = element_blank(),
        axis.title.x = element_blank()
  )
#p_R

#### ---- overall performance ---- ####


dt <- data_AUC[!is.na(base_method),.(AUC=mean(AUC_norm),max_PR=mean(max_PR),mean_PR=mean(mean_PR)),by=.(base_method,AUC_meth)]
dt[,method:=paste(base_method,AUC_meth,sep = "*' | '*")]

dt[,method:=factor(method,levels=dt[order(AUC),method])]


p_performance <- ggplot(dt,aes(method,AUC,fill=base_method))+
  geom_col(alpha=1,col="black") +
  scale_fill_manual(values = pal_jco()(10),name=NULL)+
  labs(linetype=NULL)+
  scale_x_discrete(labels = function(x) parse(text = as.character(dt[order(AUC)][!duplicated(method),AUC_meth])))+
  guides(fill = guide_legend(nrow = 1))+
  theme_bw(base_size = 20) +
  theme(
    aspect.ratio = 0.5,
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "inside",
    legend.position.inside =  c(0.4,0.92),
    legend.key.size = unit(1.3, "lines"),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.background = element_blank(),
    legend.text = element_text(size = 16)
  ) +
  ylab(expression(AUC))



pdf("./figures/Fig4.pdf",height = 10,width = 16)
grid.arrange(p_PVE+ggtitle("a) Sensitivity to parameter choice"),
             p_R+ggtitle("b) Robustness"),
             p_performance+ggtitle("c) Performance"),
             p_recovered+ggtitle("d) Detection power"),ncol=2)
dev.off()
