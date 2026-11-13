source("R/01_calcular_perfis.R")



REAJUSTAR <- FALSE

APLICACOES <- c(
  "ESCC",
  "COVID_sem_regressao",
  "COVID_com_regressao"
)

N_UNI <- 81
N_BI  <- 35


library(patchwork)
source("R/02_graficos_ggplot.R")

# Carregar os perfis da primeira aplicação
P <- readRDS("perfis/perfis_ESCC.rds")

# Perfil univariado de eta
g1 <- plot_uni(P, "eta", cor = "black",
               xlim = c(0.3, 0.85),
               ylim = c(-0.35, 30),
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )

print(g1)

g2 <- plot_uni(P, "beta", cor = "black",
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )
print(g2)

g3 <- plot_uni(P, "gamma", cor = "black",
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )
print(g3)

figura_uni_p1 <- g1|g2|g3
x11()
figura_uni_p1

ggsave(
  "ini_aplica_1.eps",
  plot = figura_uni_p1,
  width = 6,
  height = 4.5,
  device = "eps"
)



# Perfil bivariado de eta e beta
gg1 <- plot_bi(
  P, "eta", "beta",
  xlim = c(0.35, 0.80),
  ylim = c(0.02, 0.12),
  cor_emv = "black",
  rotulos = TRUE
)

print(gg1)

# Perfil bivariado de eta e beta
gg2 <- plot_bi(
  P, "eta", "gamma",
  cor_emv = "black",
  rotulos = TRUE
)

print(gg2)

# Perfil bivariado de eta e beta
gg3 <- plot_bi(
  P, "beta", "gamma",
  cor_emv = "black",
  rotulos = TRUE
)

print(gg3)


figura_bi_p1  <- gg1|gg2|gg3
figura_bi_p1

ggsave(
  "bi_aplica_1.eps",
  plot = figura_bi_p1,
  width = 6,
  height = 4.5,
  device = "eps"
)



######################33
P2 <- readRDS(
  "perfis/perfis_COVID_com_regressao.rds"
)

P2$univariados

# Perfil univariado de eta
h1 <- plot_uni(P2, "xi0", cor = "black",
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )

print(h1)

h2 <- plot_uni(P2, "xi_sex", cor = "black",
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )

print(h2)

h3 <- plot_uni(P2, "xi_age60", cor = "black",
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )

print(h3)

coeficientes <- h1|h2|h3

h4 <- plot_uni(P2, "beta", cor = "black",
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )

print(h4)

h5 <- plot_uni(P2, "gamma", cor = "black",
               cor_emv = "black")+
  ggplot2::geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "solid",
    linewidth = 0.5
  )

print(h5)

coeficientes <- (h1|h2|h3)/(h4|h5| plot_spacer())
coeficientes 

ggsave(
  "uni_aplica_2.eps",
  plot = coeficientes,
  width = 6,
  height = 4.5,
  device = "eps"
)

####


# Perfil bivariado de eta e beta
hh1 <- plot_bi(
  P2, "xi0", "beta",
  cor_emv = "black",
  rotulos = TRUE
)

print(hh1)


hh2 <- plot_bi(
  P2, "xi_sex", "beta",
  cor_emv = "black",
  rotulos = TRUE
)

print(hh2)

hh3 <- plot_bi(
  P2, "xi_age60", "beta",
  cor_emv = "black",
  rotulos = TRUE
)

print(hh3)




# Perfil bivariado de eta e beta
hh4 <- plot_bi(
  P2, "xi0", "gamma",
  cor_emv = "black",
  rotulos = TRUE
)

print(hh4)


hh5 <- plot_bi(
  P2, "xi_sex", "gamma",
  cor_emv = "black",
  rotulos = TRUE
)

print(hh5)

hh6 <- plot_bi(
  P2, "xi_age60", "gamma",
  cor_emv = "black",
  rotulos = TRUE
)

print(hh6)



# Perfil bivariado de eta e beta
hh7 <- plot_bi(
  P2, "beta", "gamma",
  cor_emv = "black",
  rotulos = TRUE
)

print(hh7)



figura_bi_p2  <- (hh1|hh2|hh3)/(hh4|hh5|hh6)/hh7
figura_bi_p2


figura_bi_p2 <- 
  (hh1 | hh2 | hh3) /
  (hh4 | hh5 | hh6) /
  (hh7 | plot_spacer() | plot_spacer())

figura_bi_p2
ggsave(
  "bi_aplica_2.eps",
  plot = figura_bi_p2,
  width = 6,
  height = 4.5,
  device = "eps"
)


