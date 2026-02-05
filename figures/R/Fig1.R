######################################################
##  Schematic figure of quantile trasformation based on made up values 
##  Main figure 2
######################################################
library(ggplot2)
library(data.table)

# Quantile grid
q <- seq(0, 1, length.out = 300)

# Theoretical null (1:1)
F0 <- q

# Permuted LD-weighted null (systematic inflation)
Fp <- q + 0.7 * q^2

# Observed LD-weighted statistic (true signal + LD)
Fn <- q + 0.7 * q^2 + 1.2 * q^4

dt <- data.table(
  q  = q,
  F0 = F0,
  Fp = Fp,
  Fn = Fn
)

dt_ld_distortion <- data.table(
  q   = q,
  ymin = F0,
  ymax = Fp
)

dt_signal <- data.table(
  q   = q,
  ymin = Fp,
  ymax = Fn
)

p1 <- ggplot() +
  
  # LD-induced distortion
  geom_ribbon(
    data = dt_ld_distortion,
    aes(x = q, ymin = ymin, ymax = ymax),
    fill = "grey70",
    alpha = 0.5
  ) +
  
  # Selection-consistent excess
  geom_ribbon(
    data = dt_signal,
    aes(x = q, ymin = ymin, ymax = ymax),
    fill = "#D95F02",
    alpha = 0.5
  ) +
  
  # Lines
  geom_line(data = dt, aes(q, F0),
            linetype = "dashed", color = "grey30", linewidth = 0.5) +
  geom_line(data = dt, aes(q, Fp),
            color = "grey50", linewidth = 0.8) +
  geom_line(data = dt, aes(q, Fn),
            color = "black", linewidth = 1) +
  
  # Labels
  labs(
    x = "Expected quantiles under null",
    y = "Observed quantiles"
  ) +
  
  coord_equal() +
  
  theme_bw(base_size = 10) +
  theme(
    aspect.ratio = 1,
    panel.grid = element_blank(),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 9)
  )+ 
  annotate("text", x = 0.9, y = 1.13,
           label = "LD-induced\ndistortion",
           size = 3, color = "grey20")+ 
  annotate("text", x = 0.91, y = 1.75,
           label = "Selection-\nconsistent\nexcess",
           size = 3, color = "#8C2D04")

pdf("./figures/Fig1.pdf", height = 5.5,width = 5.5)
p1
dev.off()
