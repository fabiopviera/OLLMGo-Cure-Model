# OLLMGo and MGo population hazards in original time units.
# Install once if needed: install.packages(c("ggplot2", "patchwork"))
library(ggplot2)
library(patchwork)

h_ollmgo <- function(t, alpha, beta, gamma) {
  q <- -(alpha / beta) * expm1(-beta * t)
  log_x <- log(expm1(q))
  z <- gamma * log_x
  log_den <- pmax(z, 0) + log1p(exp(-abs(z)))
  exp(log(gamma * alpha) - beta * t + q + (gamma - 1) * log_x - log_den)
}

A <- function(q, gamma) {
  x <- expm1(q)
  (gamma - 1 + gamma * x - x^gamma) / (x * (1 + x^gamma))
}

# Boundary for gamma = 0.5. The label is rounded, the calculation is not.
q0 <- log1p((1 + sqrt(2))^2)
boundary <- optimize(
  function(q) q + 1 / A(q, 0.5),
  interval = c(q0 + 1e-6, 8), tol = 1e-12
)
q_m <- boundary$minimum
lambda_c <- boundary$objective

original_time <- function(q, lambda, beta) -log1p(-q / lambda) / beta

critical_time <- function(lambda, beta, gamma, interval) {
  q <- uniroot(function(q) (lambda - q) * A(q, gamma) - 1,
               interval, tol = 1e-12)$root
  original_time(q, lambda, beta)
}

t_B <- original_time(q_m, lambda_c, 1)
t_C <- c(critical_time(35, 0.5, 0.5, c(q0, q_m)),
         critical_time(35, 0.5, 0.5, c(q_m, 35)))
t_D <- critical_time(2, 1.5, 2, c(1e-6, 2))

panel <- function(title, alpha, beta, gamma, t_max,
                  points = numeric(0), zoom = NULL, alpha_digits = 8) {
  t <- exp(seq(log(0.007 / beta), log(t_max), length.out = 1500))
  t_mgo <- c(0, t)
  curves <- rbind(
    data.frame(t, hazard = h_ollmgo(t, alpha, beta, gamma), model = "OLLMGo"),
    data.frame(t = t_mgo, hazard = alpha * exp(-beta * t_mgo), model = "MGo")
  )
  critical <- data.frame(t = points,
                         hazard = h_ollmgo(points, alpha, beta, gamma))
  
  # Round only the displayed alpha, retaining the full value for calculations.
  alpha_show <- signif(alpha, alpha_digits)
  alpha_term <- if (alpha_show == alpha) {
    bquote(alpha == .(alpha_show))
  } else {
    bquote(alpha %~~% .(alpha_show))
  }
  subtitle <- bquote(.(alpha_term) ~ "," ~
                       beta == .(beta) ~ "," ~ gamma == .(gamma))
  
  p <- ggplot(curves, aes(t, hazard)) +
    geom_line(aes(linetype = model), colour = "black", linewidth = 0.65) +
    geom_point(data = critical, shape = 21, size = 2.3,
               fill = "black", colour = "white", stroke = 0.35) +
    scale_linetype_manual(
      values = c(OLLMGo = "solid", MGo = "dashed"),
      breaks = c("OLLMGo", "MGo"),
      labels = expression("OLLMGo", paste("MGo (", gamma == 1, ")"))
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
    coord_cartesian(xlim = c(0, t_max), ylim = c(0, NA)) +
    labs(title = title, subtitle = subtitle, linetype = NULL,
         x = expression(paste( t)),
         y = expression(paste(h(t)))) +
    theme_classic(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.margin = margin(9, 10, 10, 9)
    )
  
  if (!is.null(zoom)) {
    # The inset isolates the OLLMGo turning behaviour.
    # Both models remain visible in the main panel.
    tz <- seq(zoom[1], zoom[2], length.out = 400)
    detail <- data.frame(t = tz, hazard = h_ollmgo(tz, alpha, beta, gamma))
    inset <- ggplot(detail, aes(t, hazard)) +
      geom_line(colour = "black", linewidth = 0.55) +
      geom_point(data = critical, size = 1.6, colour = "black") +
      scale_x_continuous(n.breaks = 4) +
      scale_y_continuous(n.breaks = 4) +
      labs(title = "OLLMGo detail", x = NULL, y = NULL) +
      theme_classic(base_size = 8) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid       = element_blank()
      )
    
    p <- p + inset_element(inset, 0.51, 0.36, 0.98, 0.96, align_to = "panel")
  }
  p
}

plots <- list(
  panel("(a)", 4, 2, 0.5, 2.5),
  panel("(b)", lambda_c, 1, 0.5, 5,
        points = t_B, zoom = c(0.055, 0.28), alpha_digits = 2),
  panel("(c)", 17.5, 0.5, 0.5, 8,
        points = t_C, zoom = c(0.09, 0.42)),
  panel("(d)", 3, 1.5, 2, 3.4, points = t_D)
)

figure <- wrap_plots(plots, ncol = 2, guides = "collect") +
  plot_annotation(
    title = "",
    subtitle = "",
    caption = paste("",
                    ""),
    theme = theme(plot.title = element_text(face = "bold", size = 18),
                  plot.subtitle = element_text(size = 11, margin = margin(b = 8)),
                  plot.caption = element_text(hjust = 0, size = 9),
                  plot.margin = margin(13, 14, 13, 14))
  )
figure <- figure & theme(legend.position = "bottom")

print(figure)


ggsave("OLLMGo_MGo_hazards.eps", figure, width = 11, height = 8,
       units = "in", dpi = 400, bg = "white")
ggsave("OLLMGo_MGo_hazards.pdf", figure, width = 11, height = 8,
       units = "in", bg = "white")







h_ollmgo <- function(t, alpha, beta, gamma) {
  q <- -(alpha / beta) * expm1(-beta * t)
  log_x <- log(expm1(q))
  z <- gamma * log_x
  log_den <- pmax(z, 0) + log1p(exp(-abs(z)))
  exp(log(gamma * alpha) - beta * t + q + (gamma - 1) * log_x - log_den)
}

curve(h_ollmgo(x,alpha=4,beta=2,gamma=0.5),xlim=c(0,1))



curve(h_ollmgo(x,alpha=25,beta=1,gamma=0.5),xlim=c(0,1))


curve(h_ollmgo(x,alpha=17.5,beta=0.5,gamma=0.5),xlim=c(0.051,1),ylim=c(5.5,8))

curve(h_ollmgo(x,alpha=3,beta=1.5,gamma=2),xlim=c(0.1,3),ylim=c(0,8))



# Zoom no mínimo e no máximo locais
curve(
  h_ollmgo(x, alpha = 17.5, beta = 0.5, gamma = 0.5),
  xlim = c(0.09, 0.42),
  ylim = c(6.78, 7.05),
  n = 1000,
  xaxs = "i", yaxs = "i",
  lwd = 2, col = "black",
  xlab = expression(paste("Time, ", t)),
  ylab = expression(paste("Hazard, ", h(t))),
  main = "Local minimum and maximum"
)

# Tempos dos dois extremos, obtidos numericamente
t_ext <- c(0.1398077, 0.2986140)
h_ext <- h_ollmgo(t_ext, alpha = 17.5, beta = 0.5, gamma = 0.5)

points(t_ext, h_ext, pch = 19)
text(t_ext, h_ext, c("Local minimum", "Local maximum"),
     pos = c(3, 1), cex = 0.85)