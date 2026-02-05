## ======================================================================
## Figures: Fig S3, Fig S6, Fig 3
## Second pass: minimal libs, annotations, bug fixes, dt0 loading provided
## ======================================================================

## ---- Libraries --------------------------------------------------------
library(data.table)   # core data munging (fast + memory-friendly)
library(ggplot2)      # plotting
library(effectsize)   # eta_squared()
library(scales)       # alpha(), number_format()
library(gridExtra)    # grid.arrange()
library(ggsci)        # pal_jco() palette

## ---- Helpers ----------------------------------------------------------

# Bin/round parameters for tiling / grouping (keeps grids manageable).
prep_param_cols <- function(dt) {
  dt <- copy(dt)
  
  # alpha: show on -log10 scale (often tiny); round for manageable grids
  if ("alpha" %in% names(dt)) dt[, alpha := round(-log10(alpha), 1)]
  
  # continuous params: round for grid size reduction
  for (nm in c("rho_w", "rho_ld", "tau_C")) {
    if (nm %in% names(dt)) dt[, (nm) := round(get(nm), 2)]
  }
  
  # rho_d: your original code used ~0.5 resolution
  if ("rho_d" %in% names(dt)) dt[, rho_d := round(rho_d * 2, 2) / 2]
  
  # l_min should be integer bins
  if ("l_min" %in% names(dt)) dt[, l_min := as.integer(l_min)]
  
  dt
}


## ======================================================================
## Fig S3: Effect sizes (eta^2) + violins of demographic/adaptation stats
## ======================================================================

if (length(list.files("./data/")) == 0) {
  message("Please download data from Zenodo and place it in ./data/")
}

pop_stats <- readRDS("./data/pop_stats.rds")
setDT(pop_stats)

## Variables to analyze in Fig S3
vars_to_test <- c(
  "adlt.fst","adlt.q1.Va","adlt.q1.Qst","Qst/Fst",
  "allelic_values","Va","a","b","n_QTN","cor_g_env",
  "a_QTN/a_ntr","fitness.mean"
)

## Compute partial eta^2 for each response variable using:
## value ~ gene_flow * selection_intensity
eta_dt <- rbindlist(lapply(vars_to_test, function(v) {
  eta <- pop_stats[variable == v,
                   effectsize::eta_squared(
                     aov(value ~ gene_flow * selection_intensity),
                     alternative = "two.sided"
                   )]
  eta$variable <- v
  as.data.table(eta)
}), fill = TRUE)

## Map ANOVA terms to parsed labels (for legend)
eta_dt[, Par := c("sigma[S]^2", "k", "sigma[S]^2%*%k")[
  match(Parameter, c("selection_intensity","gene_flow","gene_flow:selection_intensity"))
]]

## Map response variables to parsed labels (for x-axis + facets)
var_map <- c(
  "adlt.fst"        = "F[ST]",
  "adlt.q1.Va"      = "V[A]",
  "adlt.q1.Qst"     = "Q[ST]",
  "Qst/Fst"         = "Q[ST]/F[ST]",
  "allelic_values"  = "bar(beta)",
  "n_QTN"           = "bar(n)[QTN]",
  "cor_g_env"       = "r[g%~%e]^2",
  "a_QTN/a_ntr"     = "a[QTN]/a[ntrl]",
  "a"               = "bar(a)",
  "b"               = "b",
  "Va"              = "bar(p)[V[A]*','*l]",
  "fitness.mean"    = "bar(w)"
)

eta_dt[, Var := unname(var_map[variable])]

## Order variables by total partial eta^2 across terms (visual importance)
eta_dt[, total_eta := sum(Eta2_partial), by = Var]
eta_dt[, Var := factor(Var, levels = unique(Var[order(-total_eta)]))]

pd <- position_dodge(width = 0.9, preserve = "single")

## (b) Effect sizes barplot
eta_dt[, tmp := ""]
p_eta <- ggplot(eta_dt, aes(Var, Eta2_partial, fill = Par)) +
  geom_col(position = pd) +
  facet_grid(. ~ tmp) +
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
    legend.position.inside = c(0.80, 0.80),
    plot.title = element_text(size = 13),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  ylab(expression(italic(eta["partial"]^2)))

## (a) Violin distributions for the same variables
pop_stats[, Var := unname(var_map[variable])]
pop_stats[, Var := factor(Var, levels = levels(eta_dt$Var))]

## Aggregate to mean(value) per replicate/chr etc. (matches your original intent)
dt_violin <- pop_stats[!is.na(Var),
                       .(value = mean(value)),
                       by = .(Var, gene_flow, selection_intensity, rep, Chr_9sp, chr)]

p_violin <- ggplot(dt_violin, aes(selection_intensity, value, fill = gene_flow)) +
  geom_violin(
    scale = "area",
    draw_quantiles = c(0.5),
    adjust = 1.5
  ) +
  facet_wrap(Var ~ ., scales = "free", nrow = 2, labeller = label_parsed) +
  scale_fill_manual(
    values = alpha(pal_jco()(10)[1:3], 0.75),
    name = "Gene flow",
    labels = function(x) parse(text = x)
  ) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank(),
    legend.position = "bottom",
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

## Save Fig S3
pdf("./Figures/FigS3.pdf", width = 14, height = 5.5)
grid.arrange(
  p_violin + ggtitle("a) Effects of simulation parameters on population demography and local adaptation"),
  p_eta    + ggtitle("b) Effect sizes"),
  nrow = 1, widths = c(1, 0.48)
)
dev.off()

## ======================================================================
## dt0: very large file (20GB). Load once, then reduce columns.
## ======================================================================

## You provided:
dt0 <- readRDS("./data/all_draws.rds")  ## very large file!! (20gb)
setDT(dt0)

## Keep only rows needed for plots and prepare parameter columns
dt0 <- prep_param_cols(dt0[!is.na(PR)])

## Optional but recommended: drop unused columns to save RAM.
## Keep what Fig S6 needs:

need_cols <- c("TP","FP","FN","PR",
               "rho_w","alpha","l_min","tau_C","rho_d","rho_ld",
               "Method","AUC_meth","base_method")
dt0 <- dt0[, ..need_cols]

## ======================================================================
## Fig S6: Performance metrics as functions of threshold parameters
## ======================================================================

## Create binned summaries for each parameter separately.
## NOTE: we keep the same selection rule from your code:
## exclude base_method in c("EMX","LFMM") for alpha/rho_w/l_min panels.

dt_s6 <- rbind(
  dt0[!(base_method %in% c("'EMX'","'LFMM'")),
      .(TP = mean(TP), FP = mean(FP), FN = mean(FN), PR = mean(PR), parameter = "rho_w"),
      by = .(value = round(rho_w * 20) / 20, Method, AUC_meth, base_method)],
  
  dt0[!(base_method %in% c("'EMX'","'LFMM'")),
      .(TP = mean(TP), FP = mean(FP), FN = mean(FN), PR = mean(PR), parameter = "alpha"),
      by = .(value = round(alpha * 2) / 2, Method, AUC_meth, base_method)],
  
  dt0[!(base_method %in% c("'EMX'","'LFMM'")),
      .(TP = mean(TP), FP = mean(FP), FN = mean(FN), PR = mean(PR), parameter = "l_min"),
      by = .(value = l_min, Method, AUC_meth, base_method)],
  
  dt0[, .(TP = mean(TP), FP = mean(FP), FN = mean(FN), PR = mean(PR), parameter = "tau_C"),
      by = .(value = round(tau_C * 20) / 20, Method, AUC_meth, base_method)],
  
  dt0[, .(TP = mean(TP), FP = mean(FP), FN = mean(FN), PR = mean(PR), parameter = "rho_d"),
      by = .(value = round(rho_d * 20) / 20, Method, AUC_meth, base_method)],
  
  ## BUG FIX: this should bin rho_ld (not rho_d)
  dt0[, .(TP = mean(TP), FP = mean(FP), FN = mean(FN), PR = mean(PR), parameter = "rho_ld"),
      by = .(value = round(rho_ld * 20) / 20, Method, AUC_meth, base_method)],
  
  fill = TRUE, use.names = TRUE
)

## Parameter label map (parsed expressions)
par_map <- data.table(
  parameter = c("tau_C","rho_w","alpha","l_min","rho_d","rho_ld"),
  Par = c("tau[C]","w[rho]","alpha","l[min]","tau[d*','*rho]","tau[LD*','*rho]")
)
dt_s6 <- merge(dt_s6, par_map, by = "parameter", all.x = TRUE)

## Order base methods (matches your figure convention)
dt_s6[, base_method := factor(base_method, levels = c("'LFMM'","'LFMM´'","'EMX'","'EMX´'","'Joint'"))]

## Derived: number of outlier regions n[OR] = TP + FP
dt_s6[, "n[OR]" := TP + FP]

## Long format: plot multiple metrics in a grid
dt_s6_long <- melt(
  dt_s6,
  measure.vars = c("n[OR]","TP","FP","FN","PR"),
  value.name = "N"
)

dt_s6_long[, Par := factor(Par, levels = c("w[rho]","alpha","l[min]","tau[d*','*rho]","tau[LD*','*rho]","tau[C]"))]

## Harmonize AUC method label
dt_s6_long[AUC_meth == "PR^*", AUC_meth := "PR^'*'"]



p_s6 <- ggplot(dt_s6_long,
               aes(
                 value, N,
                 col = base_method,
                 group=interaction(Method,AUC_meth,base_method),
                 linetype=AUC_meth)) +
  geom_line(linewidth = 0.5) +
  facet_wrap(variable ~ Par, scales = "free", labeller = label_parsed, ncol = 6) +
  scale_linetype_manual(values = c(2, 1), name = NULL, labels = function(x) parse(text = x)) +
  scale_color_manual(values = pal_jco()(10), name = NULL, labels = function(x) parse(text = x)) +
  scale_x_continuous(labels = number_format(accuracy = 0.1)) +
  theme_bw(base_size = 15) +
  theme(
    aspect.ratio = 0.75,
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

pdf("./Figures/FigS6.pdf", width = 14, height = 14)
p_s6
dev.off()

## ======================================================================
## Fig 3: Violin plots of key summary variables
## ======================================================================

## Reuse pop_stats, but only for the three panels shown in Fig 3.
## Use the same Var mapping scheme, with a corrected label for background LD.

dt_fig3 <- pop_stats[
  Var %in% c("Q[ST]/F[ST]", "b", "a[QTN]/a[ntrl]"),
  .(value = mean(value)),
  by = .(Var, gene_flow, selection_intensity, rep, Chr_9sp, chr)
]

## BUG FIX: correct spelling; keep parsed-friendly quoting if you rely on parse()
dt_fig3[Var == "b", Var := "'Background LD'"]

dt_fig3[, Var := factor(Var, levels = c("Q[ST]/F[ST]", "'Background LD'", "a[QTN]/a[ntrl]"))]

dt_hline <- data.table(
  Var = c("Q[ST]/F[ST]", "'Background LD'", "a[QTN]/a[ntrl]"),
  yintercept = c(1, NA, 1)
)
dt_hline[, Var := factor(Var, levels = levels(dt_fig3$Var))]

p_fig3 <- ggplot(dt_fig3, aes(selection_intensity, value, fill = gene_flow)) +
  geom_violin(scale = "width", draw_quantiles = c(0.5), adjust = 1) +
  facet_wrap(Var ~ ., scales = "free", nrow = 1, labeller = label_parsed) +
  scale_fill_manual(
    values = alpha(pal_jco()(10)[1:3], 0.75),
    name = "Gene flow",
    labels = function(x) parse(text = x)
  ) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  geom_hline(
    data = dt_hline,
    aes(yintercept = yintercept),
    linetype = 2, colour = "black", linewidth = 0.5
  ) +
  theme_bw() +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 90),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank()
  )

pdf("./Figures/Fig3.pdf", height = 3, width = 9)
p_fig3
dev.off()
