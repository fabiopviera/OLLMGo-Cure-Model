# Figura do estudo de simulacao OLLMGo (R base, sem pacotes externos).
# Uso: Rscript OLLMGo_figura_simulacao.R pasta_dos_RData pasta_de_saida
# Sem argumentos, le os dados e grava as figuras no diretorio atual.

library(grid)

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[1L] else "."
output_dir <- if (length(args) >= 2L) args[2L] else "."

scenarios <- 1:8
sizes <- c(200, 500, 800, 1000)
metrics <- c(Rbias = "RB", ESE = "MCse", ASE = "Ase", RMSE = "RMSE", CP = "CP")
parameters <- c("gama0", "gama1", "gama2", "b", "g")
parameter_labels <- expression(xi[0], xi[1], xi[2], beta, gamma)
columns <- as.vector(outer(unname(metrics), parameters, paste, sep = "."))

read_table <- function(scenario) {
  path <- file.path(input_dir, sprintf("simulation_results/Res-sc%d.RData", scenario))
  object_name <- paste0("Res", scenario)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object_name, envir = env, inherits = FALSE)) {
    stop("Objeto ausente em ", path, ": ", object_name)
  }
  tab <- get(object_name, envir = env, inherits = FALSE)$Table
  if (!all(c("n", columns) %in% colnames(tab))) {
    stop("Colunas esperadas ausentes em ", path)
  }
  if (nrow(tab) != length(sizes) || !setequal(tab[, "n"], sizes)) {
    stop("Tamanhos amostrais inesperados em ", path)
  }
  tab <- tab[match(sizes, tab[, "n"]), columns, drop = FALSE]
  if (!all(is.finite(tab))) stop("Valores nao finitos em ", path)
  tab
}

tables <- lapply(scenarios, read_table)

# Escalas comuns aos cinco parametros de cada metrica.
ticks <- list(
  Rbias = c(-0.2, -0.1, 0, 0.1, 0.2),
  ESE = c(0, 0.5, 1.0, 1.5),
  ASE = c(0, 0.4, 0.8, 1.2),
  RMSE = c(0, 0.5, 1.0, 1.5),
  CP = c(0.92, 0.95, 0.98)
)
ylimits <- list(
  Rbias = c(-0.28, 0.22),
  ESE = c(-0.07, 1.9),
  ASE = c(-0.05, 1.31),
  RMSE = c(-0.07, 1.9),
  CP = c(0.915, 0.985)
)
references <- c(Rbias = 0, CP = 0.95)

colors <- c("#202020", "#505050", "#272727", "#646464",
            "#303030", "#777777", "#414141", "#555555")
symbols <- c(21, 22, 24, 23, 25, 22, 21, 23)
filled <- c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE)
line_types <- c("solid", "44", "solid", "22", "44", "solid", "13", "44")
line_style <- function(i) gpar(col = colors[i], lwd = 0.90, lty = line_types[i])
point_style <- function(i) {
  gpar(col = colors[i], fill = if (filled[i]) colors[i] else "white", lwd = 0.70)
}

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(gp = gpar(fontfamily = "sans", col = "#202020")))
  # Coordenadas em polegadas para manter a geometria identica nos dois formatos.
  inch <- function(x) unit(x, "inches")
  label <- function(text, x, y, fontsize = 8, ...) {
    grid.text(text, inch(x), inch(y), gp = gpar(fontsize = fontsize), ...)
  }
  strip <- function(x, y, width, height) {
    grid.rect(inch(x), inch(y), inch(width), inch(height),
              just = c("left", "bottom"),
              gp = gpar(fill = "#EFEFEF", col = "#AEAEAE", lwd = 0.55))
  }
  
  x_left <- 0.43
  x_right <- 6.76
  y_top <- 5.20
  y_bottom <- 0.81
  gap_x <- 0.082
  gap_y <- 0.080
  panel_width <- (x_right - x_left - 4 * gap_x) / 5
  panel_height <- (y_top - y_bottom - 4 * gap_y) / 5
  xlim <- c(115, 1085)
  map_x <- function(n) (n - xlim[1L]) / diff(xlim)
  axis_gp <- gpar(col = "#777777", lwd = 0.55)
  
  for (row in seq_along(metrics)) {
    metric <- names(metrics)[row]
    limits <- ylimits[[metric]]
    bottom <- y_top - row * panel_height - (row - 1L) * gap_y
    y_ticks <- ticks[[metric]]
    y_labels <- sprintf(if (metric == "CP") "%.2f" else "%.1f", y_ticks)
    strip(x_right + 0.045, bottom, 0.28, panel_height)
    label(metric, x_right + 0.185, bottom + panel_height / 2, rot = 270)
    
    for (column in seq_along(parameters)) {
      left <- x_left + (column - 1L) * (panel_width + gap_x)
      series_name <- paste(metrics[[metric]], parameters[column], sep = ".")
      values <- vapply(tables, function(tab) tab[, series_name], numeric(length(sizes)))
      if (any(values < limits[1L] | values > limits[2L])) {
        stop("Ajuste os limites do eixo para a metrica ", metric)
      }
      if (row == 1L) {
        strip(left, y_top + 0.035, panel_width, 0.23)
        label(parameter_labels[column], left + panel_width / 2, y_top + 0.15, 10)
      }
      
      pushViewport(viewport(inch(left), inch(bottom), inch(panel_width), inch(panel_height),
                            just = c("left", "bottom"), xscale = xlim, yscale = limits))
      grid.segments(unit(0, "npc"), unit(y_ticks, "native"),
                    unit(1, "npc"), unit(y_ticks, "native"),
                    gp = gpar(col = "#EAEAEA", lwd = 0.50))
      if (metric %in% names(references)) {
        y_ref <- unit(references[[metric]], "native")
        grid.segments(unit(0, "npc"), y_ref, unit(1, "npc"), y_ref,
                      gp = gpar(col = "#6A6A6A", lwd = 0.75, lty = "33"))
      }
      for (i in scenarios) {
        grid.lines(unit(sizes, "native"), unit(values[, i], "native"), gp = line_style(i))
      }
      for (i in scenarios) {
        grid.points(unit(sizes, "native"), unit(values[, i], "native"),
                    pch = symbols[i], size = unit(3.25, "points"), gp = point_style(i))
      }
      grid.rect(gp = gpar(fill = NA, col = "#AEAEAE", lwd = 0.55))
      
      if (column == 1L) {
        grid.segments(unit(-0.035, "inches"), unit(y_ticks, "native"),
                      unit(0, "npc"), unit(y_ticks, "native"), gp = axis_gp)
        grid.text(y_labels, unit(-0.055, "inches"), unit(y_ticks, "native"),
                  just = "right", gp = gpar(fontsize = 7.6))
      }
      if (row == length(metrics)) {
        grid.segments(unit(sizes, "native"), unit(0, "npc"),
                      unit(sizes, "native"), unit(-0.035, "inches"), gp = axis_gp)
      }
      popViewport()
      
      if (row == length(metrics)) {
        label(as.character(sizes), left + panel_width * map_x(sizes), y_bottom - 0.115, 7.3)
      }
    }
  }
  
  label(expression(italic(n)), (x_left + x_right) / 2, 0.49, 10)
  key_y <- 0.22
  label("Scenarios", 1.39, key_y, 8.5, just = "right")
  for (i in scenarios) {
    key_x <- 1.70 + (i - 1L) * 0.57
    grid.segments(inch(key_x - 0.11), inch(key_y), inch(key_x + 0.11), inch(key_y),
                  gp = line_style(i))
    grid.points(inch(key_x), inch(key_y), pch = symbols[i], size = unit(4.0, "points"),
                gp = point_style(i))
    label(as.character(i), key_x + 0.16, key_y, 8.5, just = "left")
  }
  popViewport()
  invisible(NULL)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
width <- 7.2
height <- 5.55
pdf_file <- file.path(output_dir, "OLLMGo_figura_simulacao.pdf")
eps_file <- file.path(output_dir, "OLLMGo_figura_simulacao.eps")
png_file <- file.path(output_dir, "OLLMGo_figura_simulacao.png")

if (capabilities("cairo")) {
  cairo_pdf(pdf_file, width = width, height = height, bg = "white", antialias = "gray")
} else {
  pdf(pdf_file, width = width, height = height, bg = "white", useDingbats = FALSE)
}
draw_figure()
invisible(dev.off())

# Em cairo_ps(), onefile = FALSE gera um arquivo EPS de uma pagina.
if (capabilities("cairo")) {
  cairo_ps(eps_file, width = width, height = height, onefile = FALSE,
           bg = "white", antialias = "gray")
} else {
  postscript(eps_file, width = width, height = height, onefile = FALSE,
             horizontal = FALSE, paper = "special", bg = "white", colormodel = "gray")
}
draw_figure()
invisible(dev.off())

# Quando disponivel, renderiza o PDF em cinza a 600 dpi. No R para Windows,
# o dispositivo PNG nativo tambem permite gerar a figura sem o Poppler.
pdftoppm <- Sys.which("pdftoppm")
if (nzchar(pdftoppm)) {
  png_prefix <- tools::file_path_sans_ext(png_file)
  status <- system2(pdftoppm, c("-r", "600", "-gray", "-png", "-singlefile",
                                shQuote(pdf_file), shQuote(png_prefix)))
  if (status != 0L) stop("Falha ao converter o PDF em PNG")
} else {
  png(png_file, width = width, height = height, units = "in", res = 600, bg = "white")
  draw_figure()
  invisible(dev.off())
}
