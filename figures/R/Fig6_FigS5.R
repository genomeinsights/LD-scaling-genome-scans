## ======================================================================
## FigS5 only: heatmaps of PR q0.9 surfaces + interaction screening (ΔAIC)
## ======================================================================

library(data.table)
library(ggplot2)
library(scales)
library(viridis)     # for scale_fill_viridis_c()
library(cowplot)     # for plot_grid()
library(gridExtra)   # only used for arranging columns (optional)

if(length(list.files("./data/"))==0) message("Please download data form Zenond")

## -----------------------------
## 1) Load data
## -----------------------------
dt0 <- readRDS("./data/all_draws.rds") ## very large file!! (20gb)
setDT(dt0)

## -----------------------------
## 2) Small utilities
## -----------------------------

# Robust 0–1 scaling (within facet panels), used for tile fills if desired.
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || diff(rng) == 0) return(rep(NA_real_, length(x)))
  (x - rng[1]) / diff(rng)
}

# Discretize/round parameters so the heatmap grid is manageable.
prep_param_cols <- function(dt) {
  dt <- copy(dt)
  
  # alpha is typically tiny: tile nicely on -log10 scale
  if ("alpha" %in% names(dt)) dt[, alpha := round(-log10(alpha), 1)]
  
  # continuous parameters rounded for a coarser grid
  for (nm in c("rho_w", "rho_ld", "tau_C")) {
    if (nm %in% names(dt)) dt[, (nm) := round(get(nm), 2)]
  }
  
  # rho_d appears to be used on a ~0.5 grid in your original code
  if ("rho_d" %in% names(dt)) dt[, rho_d := round(rho_d * 2, 2) / 2]
  
  # l_min should be integer bins
  if ("l_min" %in% names(dt)) dt[, l_min := as.integer(l_min)]
  
  dt
}

## -----------------------------
## 3) Heatmap generator (q0.9 of PR)
## -----------------------------

# Expression labels used in the figure
axis_lab <- c(
  alpha  = expression(-log[10](alpha)),
  l_min  = expression(l[min]),
  rho_w  = expression(w[rho]),
  rho_ld = expression(LD[rho]),
  rho_d  = expression(d[rho]),
  tau_C  = expression(C)
)

make_pair_heat <- function(dt, xvar, yvar,
                           by_vars = c("base_method", "sim"),
                           measure_var = "PR",
                           q = 0.9) {
  
  stopifnot(all(c(xvar, yvar, by_vars, measure_var) %in% names(dt)))
  
  heat <- dt[
    ,
    .(PR_q90 = quantile(get(measure_var), q, na.rm = TRUE)),
    by = c(by_vars, xvar, yvar)
  ]
  
  # Scale within each (base_method, sim) panel for visual comparability
  heat[, PR_q90_01 := scale01(PR_q90), by = by_vars]
  heat[is.na(PR_q90_01), PR_q90_01 := 0]
  
  # Match your base_method ordering (edit levels to match your data exactly)
  if ("base_method" %in% names(heat)) {
    heat[, base_method := factor(
      base_method,
      levels = rev(c("EMX", "LFMM", "EMX´", "LFMM´", "Joint"))
    )]
  }
  
  ggplot(heat, aes(x = .data[[xvar]], y = .data[[yvar]], fill = PR_q90_01)) +
    geom_tile(linewidth = 0.2) +
    scale_fill_viridis_c(
      option = "turbo",
      name = expression(PR[0.9]),
      na.value = "grey90"
    ) +
    facet_grid(reformulate("sim", "base_method")) +
    scale_x_continuous(labels = number_format(accuracy = 0.1)) +
    theme_minimal(base_size = 8) +
    theme(
      aspect.ratio = 1,
      strip.text = element_text(margin = margin(1, 1, 1, 1)),
      plot.margin = margin(1, 1, 1, 1),
      panel.spacing = unit(0.05, "lines"),
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.background = element_blank(),
      panel.grid.major = element_blank(),
      legend.position = "right"
    ) +
    xlab(axis_lab[[xvar]] %||% xvar) +
    ylab(axis_lab[[yvar]] %||% yvar)
}

## -----------------------------
## 4) Build the *specific* four heatmaps used in FigS5
##    (based on the indices you used originally)
## -----------------------------

dt0 <- prep_param_cols(dt0[!is.na(PR)])

dt0[,base_method:=gsub("'","",base_method)]
# STAR set: AUC_meth == "PR^'*'" and you used:
#   [[4]] = (l_min vs rho_w)
#   [[2]] = (alpha vs rho_w)
dt_star <- dt0[AUC_meth == "PR^'*'"]

g_b <- make_pair_heat(dt_star, xvar = "l_min", yvar = "rho_w")   + ggtitle("b)")
g_c <- make_pair_heat(dt_star, xvar = "alpha", yvar = "rho_w")   + ggtitle("c)")

# C set: AUC_meth != "PR^'*'" and method starts with "C_"
# You used:
#   [[2]] = (tau_C vs rho_d)
#   [[1]] = (tau_C vs rho_ld)
dt_C <- dt0[AUC_meth != "PR^'*'" & grepl("^C_", method)]

g_d <- make_pair_heat(dt_C, xvar = "tau_C", yvar = "rho_d")      + ggtitle("d)")
g_e <- make_pair_heat(dt_C, xvar = "tau_C", yvar = "rho_ld")     + ggtitle("e)")

## Arrange the four heatmaps into the same layout as your original:
## left column: b over d, right column: c over e
p_heat <- plot_grid(
  plot_grid(g_b, g_d, ncol = 1, align = "v"),
  plot_grid(g_c, g_e, ncol = 1, align = "v"),
  nrow = 1
)

## -----------------------------
## 5) Interaction screening: ΔAIC for adding p1:p2 to main-effects model
## -----------------------------

# Fast LR / ΔAIC comparison for main-effects vs interaction (binomial glm.fit)
screen_pair_glm_fast <- function(d, p1, p2, y = "highPR") {
  yv <- d[[y]]
  if (is.logical(yv)) yv <- as.integer(yv)
  if (!is.numeric(yv)) yv <- as.integer(as.factor(yv)) - 1L
  
  x1 <- d[[p1]]
  x2 <- d[[p2]]
  
  # Standardize numeric predictors for stability (optional, but helps)
  if (is.numeric(x1)) x1 <- as.numeric(scale(x1))
  if (is.numeric(x2)) x2 <- as.numeric(scale(x2))
  
  keep <- !is.na(x1) & !is.na(x2) & !is.na(yv)
  x1 <- x1[keep]; x2 <- x2[keep]; yv <- yv[keep]
  
  # Reduced: 1 + x1 + x2
  Xr <- cbind(1, x1, x2)
  # Full:    1 + x1 + x2 + x1:x2
  Xf <- cbind(1, x1, x2, x1 * x2)
  
  fr <- suppressWarnings(glm.fit(Xr, yv, family = binomial()))
  ff <- suppressWarnings(glm.fit(Xf, yv, family = binomial()))
  
  LR <- fr$deviance - ff$deviance
  df <- (ff$rank - fr$rank)
  p  <- pchisq(LR, df = df, lower.tail = FALSE)
  
  AIC_red  <- fr$deviance + 2 * fr$rank
  AIC_full <- ff$deviance + 2 * ff$rank
  
  data.table(
    p1 = p1, p2 = p2,
    delta_AIC = AIC_red - AIC_full,
    LR_Chisq  = LR,
    df        = df,
    p_value   = p
  )
}

screen_all_pairs_by_stratum <- function(dt, params, keys,
                                        y = "highPR",
                                        max_n = 50000L) {
  stopifnot(all(c(keys, y, params) %in% names(dt)))
  
  # All unique parameter pairs
  pairs <- t(combn(params, 2))
  pairs <- data.table(p1 = pairs[, 1], p2 = pairs[, 2])
  
  dt[
    ,
    {
      dsub <- .SD
      if (!is.null(max_n) && nrow(dsub) > max_n) {
        dsub <- dsub[sample.int(.N, max_n)]
      }
      rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
        screen_pair_glm_fast(dsub, pairs$p1[i], pairs$p2[i], y = y)
      }))
    },
    by = keys
  ]
}

## Define "highPR" within each scenario × method (same as your original intent)
# For C methods:
dt_C2 <- dt_C[!is.na(PR)]
dt_C2[, highPR := PR >= quantile(PR, 0.9, na.rm = TRUE),
      by = .(gene_flow, selection_intensity, base_method)]

res_C <- screen_all_pairs_by_stratum(
  dt = dt_C2,
  params = c("tau_C", "rho_ld", "rho_d"),
  keys = c("gene_flow", "selection_intensity", "base_method"),
  y = "highPR",
  max_n = 50000L
)
res_C[, AUC_meth := "PR[C]"]

# For STAR methods:
dt_star2 <- dt_star[!is.na(PR)]
dt_star2[, highPR := PR >= quantile(PR, 0.9, na.rm = TRUE),
         by = .(gene_flow, selection_intensity, base_method)]

res_star <- screen_all_pairs_by_stratum(
  dt = dt_star2,
  params = c("alpha", "l_min", "rho_w", "rho_ld", "rho_d"),
  keys = c("gene_flow", "selection_intensity", "base_method"),
  y = "highPR",
  max_n = 50000L
)
res_star[, AUC_meth := "PR^'*'"]

dt_int <- rbind(res_C, res_star, fill = TRUE)

## Pretty pair labels (parsed expressions)
xlab_map <- c(
  alpha  = "alpha",
  l_min  = "l[min]",
  rho_w  = "w[rho]",
  rho_ld = "LD[rho]",
  rho_d  = "d[rho]",
  tau_C  = "C"
)

dt_int[, pair_exp := paste0(xlab_map[p1], " %*% ", xlab_map[p2])]

# Add the same panel lettering you used previously (optional, but matches FigS5)
dt_int[pair_exp == "alpha %*% l[min]",     pair_exp := "alpha %*% l[min]*'      '"]
dt_int[pair_exp == "l[min] %*% w[rho]",    pair_exp := "l[min] %*% w[rho]*' (b)'"]
dt_int[pair_exp == "alpha %*% w[rho]",     pair_exp := "alpha %*% w[rho]*' (c)'"]
dt_int[pair_exp == "C %*% d[rho]",         pair_exp := "C %*% d[rho]*' (d)'"]
dt_int[pair_exp == "C %*% LD[rho]",        pair_exp := "C %*% LD[rho]*' (e)'"]

# Facet label (gene_flow - selection_intensity), keep as in original
dt_int[, facet := interaction(gene_flow, selection_intensity, drop = TRUE, sep = " - ")]

## ΔAIC interaction plot (panel a)
p_a <- ggplot(
  dt_int,
  aes(delta_AIC, reorder(pair_exp, delta_AIC),
      col = AUC_meth, shape = base_method)
) +
  scale_y_discrete(labels = function(x) parse(text = x)) +
  geom_point(size = 2) +
  facet_wrap(~ facet, nrow = 1) +
  scale_shape_manual(values = c(20, 4, 3, 5, 6), name = "Base",labels = function(x) parse(text = x)) +
  scale_color_manual(values = pal_jco()(10)[c(4,1)],name="AUC",labels = function(x) parse(text = x))+
  labs(x = expression(Delta*AIC~"(interaction)"), y = NULL) +
  theme_bw(base_size = 10) +
  theme(
    aspect.ratio = 1.5,
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90)
  ) +
  coord_cartesian(xlim = c(-5, 500)) +
  ggtitle("a)")

## -----------------------------
## 6) Assemble and save FigS5
## -----------------------------

pdf("./figures/FigS5.pdf", width = 15, height = 12)
plot_grid(p_a, p_heat, ncol = 1, rel_heights = c(0.45, 1))
dev.off()

## -------------------
## Figure 6, main document
## ---------------

## AUC, panel a)


dt <- cbind(data_AUC[AUC_meth=="PR[C]",], C=data_AUC[AUC_meth=="PR[C]",AUC_norm],star=data_AUC[AUC_meth!="PR[C]",AUC_norm])

dt <- data_AUC[ Method %in% c("EMX-C","LFMM-C","EMX´-C","LFMM´-C","Joint-C")]


#dt[Method ==  "LFMM" , base_method := "LFMM^'*'"]
dt[AUC_meth=="PR[C]", base_method := paste0(base_method,"[C]")]


dt[,base_method:=factor(base_method,levels=c("EMX[C]",
                                             "LFMM[C]",
                                             "EMX´[C]",
                                             "LFMM´[C]",
                                             "Joint[C]"),labels=c("EMX[C]",
                                                                  "LFMM[C]",
                                                                  "'EMX´'[C]",
                                                                  "'LFMM´'[C]",
                                                                  "Joint[C]"))]

p_AUC <- ggplot(dt, aes(
  base_method, AUC_norm,
  fill = ifelse(grepl("prime", method) | grepl("joint", method),
                "LD-scaled", "Baseline")
)) +
  geom_violin(alpha = 0.75, linewidth = 0.2, adjust = 1.25, scale = "area") + 
  geom_boxplot(width = 0.2, alpha = 0.75, outlier.shape = NA,
               col = "black", linewidth = 0.2) +
  #facet_grid(.~sim) +
  facet_wrap(.~sim,nrow=1)+
  scale_fill_manual(values = alpha(pal_jco()(10), 0.75), name = NULL) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(linetype = NULL) +
  #scale_y_continuous(labels = function(x) parse(text = x)) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_bw(base_size = 12) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 35, hjust = 1),
    axis.title.x = element_blank(),
    legend.position = "inside",
    legend.background = element_blank(),
    legend.position.inside = c(0.95, 0.7),
    legend.direction = "vertical",  
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ylab(expression(AUC))
p_AUC


## heatmap, panel b)
params = c("alpha","l_min")
by_vars = c("base_method","sim")
q = 0.9
fillvar = "PR_q90_"
xvar <- "alpha"
yvar <- "l_min"

p_heat <- make_pair_heat(dt_star, xvar = "alpha", yvar = "l_min")   

#p_heat + theme(plot.title = element_text(size = 12))
## panel c)

dt <- rbind(dt0[,.(PR=mean(PR),parameter="rho_w"),by=.(value=round(rho_w*20)/20,Method,AUC_meth,base_method)],
            dt0[,.(PR=mean(PR),parameter="tau_C"),by=.(value=round(tau_C*20)/20,Method,AUC_meth,base_method)],
            fill=TRUE,use.names=TRUE)
# 


par_map <- data.table(
  parameter = c("tau_C","rho_w"),
  Par = c("tau[C]*(PR[C])","w[rho]*(PR^'*')")
)

dt <- merge(dt, par_map, by="parameter", all.x=TRUE)


dt[,base_method:=factor(base_method,levels=c("LFMM",
                                             "LFMM´",
                                             "EMX",
                                             "EMX´",
                                             "Joint"))]


p2 <- ggplot(dt[Par!="tau[C]*(PR[C])"],aes(value,PR,col = base_method,group=paste(Method,AUC_meth,base_method)))+
  geom_line(linewidth=0.5) +
  scale_linetype_manual(values = c(2,1),name=NULL)+
  scale_color_manual(values = pal_jco()(10),name=NULL)+
  labs(linetype=NULL)+
  guides(fill = guide_legend(nrow = 1))+
  theme_bw(base_size = 10)+
  scale_x_continuous(labels = number_format(accuracy = 0.1))+
  theme(
    aspect.ratio = 1.25,
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()) +
  xlab(expression(w[rho])) +
  ggtitle(expression("c) Optimal values for"~w[rho]))
#p2

p3 <- ggplot(dt[Par=="tau[C]*(PR[C])"],aes(value,PR,col = base_method,group=paste(Method,AUC_meth,base_method)))+
  geom_line(linewidth=0.5) +
  scale_linetype_manual(values = c(2,1),name=NULL)+
  scale_color_manual(values = pal_jco()(10),name=NULL)+
  labs(linetype=NULL)+
  guides(fill = guide_legend(nrow = 1))+
  theme_bw(base_size = 10)+
  scale_x_continuous(labels = number_format(accuracy = 0.1))+
  theme(
    aspect.ratio = 1.25,
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank())+
  xlab(expression(tau[C])) +
  ggtitle(expression("d) Optimal values for"~tau[C]))
#p3

p4 <- grid.arrange(p2,p3,ncol=1)
p3 <- grid.arrange(p_heat+ggtitle(expression("b) Co-dependancy between significance threshold and minimum OG size"))+theme(plot.title = element_text(size = 12)),
                   p4,nrow=1,widths=c(1,0.3))

pdf("./Figures/Fig6.pdf",width = 13,height = 9)
plot_grid(p_AUC+ggtitle(expression("a) Performance"))+theme(title = element_text(size = 11)),p3,ncol=1,rel_heights = c(0.53,1))
dev.off()



