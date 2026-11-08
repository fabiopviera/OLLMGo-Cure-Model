# install.packages(c("ggplot2", "patchwork"))
library(ggplot2)
library(patchwork)

# Equation (14): compute alpha from the target cure fraction.
alpha_from_eta <- function(eta, beta, gamma) {
  beta * log1p(((1 - eta) / eta)^(1 / gamma))
}

surv_ollmgo <- function(t, eta, beta, gamma) {
  alpha <- alpha_from_eta(eta, beta, gamma)
  q <- -(alpha / beta) * expm1(-beta * t)
  1 / (1 + expm1(q)^gamma)
}

beta <- 0.5
t <- seq(0, 15, length.out = 801)
gamma_values <- c(0.5, 1, 2, 4)
eta_values <- c(0.10, 0.20, 0.40)

left_data <- do.call(rbind, lapply(gamma_values, function(gamma) {
  data.frame(t, survival = surv_ollmgo(t, 0.20, beta, gamma),
             curve = as.character(gamma))
}))
left_data$curve <- factor(left_data$curve, levels = as.character(gamma_values))

right_data <- do.call(rbind, lapply(eta_values, function(eta) {
  data.frame(t, survival = surv_ollmgo(t, eta, beta, 2),
             curve = sprintf("%.2f", eta))
}))
right_data$curve <- factor(right_data$curve, levels = sprintf("%.2f", eta_values))

survival_panel <- function(data, limits, styles, panel, subtitle, legend) {
  ggplot(data, aes(t, survival)) +
    geom_hline(yintercept = limits, colour = "grey70",
               linetype = "dotted", linewidth = 0.35) +
    geom_line(aes(linetype = curve), colour = "black", linewidth = 0.65) +
    scale_linetype_manual(values = styles) +
    scale_x_continuous(breaks = seq(0, 15, 5)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 1)) +
    labs(title = panel, subtitle = subtitle, x = "Time",
         y = "Survival", linetype = legend) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          legend.key.width = grid::unit(1.1, "cm"),
          plot.margin = margin(8, 12, 8, 12))
}

p_left <- survival_panel(
  left_data, 0.20, c("solid", "dashed", "dotted", "dotdash"), "(a)",
  expression(paste(eta == 0.20, ",  ", beta == 0.5)), expression(gamma)
)

p_right <- survival_panel(
  right_data, eta_values, c("solid", "dashed", "dotdash"), "(b)",
  expression(paste(gamma == 2, ",  ", beta == 0.5)), expression(eta)
)

figure <- p_left + p_right + plot_layout(ncol = 2, guides = "keep")
print(figure)

# Optional exports for the manuscript.
# ggsave("OLLMGo_survival_eta.pdf", figure, width = 10, height = 4.3, units = "in")
 ggsave("OLLMGo_survival_eta.eps", figure, width = 10, height = 4.3,
       units = "in", device = "eps")
