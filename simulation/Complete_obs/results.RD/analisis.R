## ============================================================
## analyze_iter_time.R
##
## Calcula el promedio (y desvio estandar) del numero de iteraciones
## de gamlss (archivos giter_p*_*.csv) y del tiempo en segundos
## (archivos time_p*_*.csv), separado por:
##   - escenario (carpeta: nu_1_05, tau_0_5, etc. o "." si no hay subcarpetas)
##   - set de parametros (p1, p2)
##   - tamaño de muestra n (50, 100, 200, 500, o el que corresponda)
##
## Uso:
##   Rscript analyze_iter_time.R [directorio_base]
##
## Si no se indica directorio_base, usa el directorio actual (".").
## Busca de forma RECURSIVA en todas las subcarpetas.
##
## Genera:
##   - summary_iter_time_long.csv   (formato largo: una fila por escenario/p/n)
##   - summary_iter_time_wide.csv   (formato ancho: una fila por escenario/n)
##   - imprime las tablas en consola
## ============================================================

args <- commandArgs(trailingOnly = TRUE)
base_dir <- if (length(args) >= 1) args[1] else "."

if (!dir.exists(base_dir)) {
  stop("El directorio indicado no existe: ", base_dir, call. = FALSE)
}

cat("Buscando resultados en:", normalizePath(base_dir), "\n\n")

## --- 1. Encontrar todos los archivos giter_p<k>_<n>.csv de forma recursiva ---
giter_files <- list.files(base_dir,
                           pattern = "^giter_p[0-9]+_[0-9]+\\.csv$",
                           recursive = TRUE,
                           full.names = TRUE)

if (length(giter_files) == 0) {
  stop("No se encontraron archivos giter_p*_*.csv bajo ", base_dir, call. = FALSE)
}

## --- 2. Funcion auxiliar: extraer escenario, p y n desde la ruta del archivo ---
parse_info <- function(filepath) {
  fname <- basename(filepath)
  # fname tipo: giter_p1_50.csv  ->  p = "p1", n = 50
  m <- regmatches(fname, regexec("^giter_(p[0-9]+)_([0-9]+)\\.csv$", fname))[[1]]
  p_label <- m[2]
  n_val <- as.numeric(m[3])

  # escenario = nombre de la carpeta contenedora (relativo a base_dir)
  scenario_dir <- dirname(filepath)
  scenario <- basename(scenario_dir)
  if (scenario == "." || scenario == basename(normalizePath(base_dir))) {
    scenario <- "(sin subcarpeta)"
  }

  list(scenario = scenario, p = p_label, n = n_val,
       giter_path = filepath,
       time_path = file.path(scenario_dir, sub("^giter_", "time_", fname)))
}

## --- 3. Leer un vector numerico desde un csv de una sola columna ---
read_vec <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  df <- tryCatch(read.csv(path), error = function(e) NULL)
  if (is.null(df) || ncol(df) == 0) return(NA_real_)
  as.numeric(df[[1]])
}

## --- 4. Armar la tabla resumen (formato largo) ---
results_list <- lapply(giter_files, function(f) {
  info <- parse_info(f)

  giter_vals <- read_vec(info$giter_path)
  time_vals  <- read_vec(info$time_path)

  data.frame(
    scenario   = info$scenario,
    p          = info$p,
    n          = info$n,
    reps       = sum(!is.na(giter_vals)),
    mean_iter  = mean(giter_vals, na.rm = TRUE),
    sd_iter    = sd(giter_vals, na.rm = TRUE),
    mean_time  = mean(time_vals, na.rm = TRUE),
    sd_time    = sd(time_vals, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

summary_long <- do.call(rbind, results_list)
summary_long <- summary_long[order(summary_long$scenario, summary_long$p, summary_long$n), ]
rownames(summary_long) <- NULL

cat("===== Resumen (formato largo) =====\n")
print(summary_long, row.names = FALSE)

## --- 5. Version ancha: una fila por escenario/n, columnas para p1 y p2 ---
library(reshape2)  # viene con muchas instalaciones de R; si no esta, usamos reshape base mas abajo

wide_ok <- requireNamespace("reshape2", quietly = TRUE)

if (wide_ok) {
  iter_wide <- dcast(summary_long, scenario + n ~ p, value.var = "mean_iter")
  colnames(iter_wide)[-(1:2)] <- paste0("mean_iter_", colnames(iter_wide)[-(1:2)])

  time_wide <- dcast(summary_long, scenario + n ~ p, value.var = "mean_time")
  colnames(time_wide)[-(1:2)] <- paste0("mean_time_", colnames(time_wide)[-(1:2)])

  summary_wide <- merge(iter_wide, time_wide, by = c("scenario", "n"))
  summary_wide <- summary_wide[order(summary_wide$scenario, summary_wide$n), ]
} else {
  # Fallback sin reshape2: armar manualmente
  scenarios <- unique(summary_long$scenario)
  ns <- sort(unique(summary_long$n))
  rows <- list()
  for (sc in scenarios) {
    for (nn in ns) {
      sub <- summary_long[summary_long$scenario == sc & summary_long$n == nn, ]
      row <- data.frame(scenario = sc, n = nn)
      for (pl in sort(unique(sub$p))) {
        row[[paste0("mean_iter_", pl)]] <- sub$mean_iter[sub$p == pl]
        row[[paste0("mean_time_", pl)]] <- sub$mean_time[sub$p == pl]
      }
      rows[[length(rows) + 1]] <- row
    }
  }
  summary_wide <- do.call(rbind, rows)
}

cat("\n===== Resumen (formato ancho) =====\n")
print(summary_wide, row.names = FALSE)

## --- 6. Guardar CSVs ---
write.csv(summary_long, file = file.path(base_dir, "summary_iter_time_long.csv"), row.names = FALSE)
write.csv(summary_wide, file = file.path(base_dir, "summary_iter_time_wide.csv"), row.names = FALSE)

cat("\nArchivos guardados:\n")
cat(" -", file.path(base_dir, "summary_iter_time_long.csv"), "\n")
cat(" -", file.path(base_dir, "summary_iter_time_wide.csv"), "\n")
