## ============================================================
## analyze_iter_time_tau.R
##
## Calcula el promedio (y desvio estandar) del numero de iteraciones
## de gamlss (archivos giter_p*_*.csv) y del tiempo en segundos
## (archivos time_p*_*.csv) para los modelos RDQ1 y RDQ2, separado por:
##   - modelo (RDQ1, RDQ2)
##   - tau (carpeta: tau_0_1, tau_0_25, tau_0_5, tau_0_75, tau_0_9)
##   - set de parametros (p1, p2)
##   - tamaño de muestra n (50, 100, 200, 500)
##
## Uso:
##   Rscript analyze_iter_time_tau.R <dir_RDQ1> <dir_RDQ2>
##   Rscript analyze_iter_time_tau.R ./results.RDQ1 ./results.RDQ2
##
## Tambien acepta un solo directorio (por si queres correrlo para un
## solo modelo), o mas de dos si tuvieras mas variantes:
##   Rscript analyze_iter_time_tau.R ./results.RDQ1
##
## Genera (en el directorio actual desde donde se corre el script):
##   - summary_iter_time_tau_long.csv
##   - summary_iter_time_tau_wide.csv
##   - imprime las tablas en consola
## ============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Debe indicar al menos un directorio base.\n",
       "Uso: Rscript analyze_iter_time_tau.R <dir_RDQ1> <dir_RDQ2>\n",
       "Ejemplo: Rscript analyze_iter_time_tau.R ./results.RDQ1 ./results.RDQ2",
       call. = FALSE)
}

base_dirs <- args

for (d in base_dirs) {
  if (!dir.exists(d)) {
    stop("El directorio indicado no existe: ", d, call. = FALSE)
  }
}

## --- Funcion auxiliar: extraer info (tau, p, n) desde la ruta del archivo ---
parse_info <- function(filepath, base_dir, model_name) {
  fname <- basename(filepath)
  # fname tipo: giter_p1_50.csv  ->  p = "p1", n = 50
  m <- regmatches(fname, regexec("^giter_(p[0-9]+)_([0-9]+)\\.csv$", fname))[[1]]
  p_label <- m[2]
  n_val <- as.numeric(m[3])

  # carpeta contenedora = escenario de tau (ej. "tau_0_5")
  scenario_dir <- dirname(filepath)
  tau_folder <- basename(scenario_dir)

  # convertir "tau_0_5" -> "0.5" para que quede mas legible
  tau_label <- sub("^tau_", "", tau_folder)
  tau_label <- gsub("_", ".", tau_label)

  if (tau_folder == basename(normalizePath(base_dir))) {
    tau_label <- "(sin subcarpeta)"
  }

  list(model = model_name, tau = tau_label, p = p_label, n = n_val,
       giter_path = filepath,
       time_path = file.path(scenario_dir, sub("^giter_", "time_", fname)))
}

## --- Leer un vector numerico desde un csv de una sola columna ---
read_vec <- function(path) {
  if (!file.exists(path)) return(NA_real_)
  df <- tryCatch(read.csv(path), error = function(e) NULL)
  if (is.null(df) || ncol(df) == 0) return(NA_real_)
  as.numeric(df[[1]])
}

## --- Recorrer cada directorio base (un modelo por directorio) ---
all_results <- list()

for (base_dir in base_dirs) {
  model_name <- basename(normalizePath(base_dir))
  cat("Buscando resultados de", model_name, "en:", normalizePath(base_dir), "\n")

  giter_files <- list.files(base_dir,
                             pattern = "^giter_p[0-9]+_[0-9]+\\.csv$",
                             recursive = TRUE,
                             full.names = TRUE)

  if (length(giter_files) == 0) {
    warning("No se encontraron archivos giter_p*_*.csv bajo ", base_dir)
    next
  }

  model_results <- lapply(giter_files, function(f) {
    info <- parse_info(f, base_dir, model_name)

    giter_vals <- read_vec(info$giter_path)
    time_vals  <- read_vec(info$time_path)

    data.frame(
      model      = info$model,
      tau        = info$tau,
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

  all_results[[length(all_results) + 1]] <- do.call(rbind, model_results)
}

if (length(all_results) == 0) {
  stop("No se encontro ningun archivo giter_p*_*.csv en los directorios indicados.", call. = FALSE)
}

summary_long <- do.call(rbind, all_results)
summary_long <- summary_long[order(summary_long$model, summary_long$tau,
                                    summary_long$p, summary_long$n), ]
rownames(summary_long) <- NULL

cat("\n===== Resumen (formato largo) =====\n")
print(summary_long, row.names = FALSE)

## --- Version ancha: una fila por modelo/tau/n, columnas para p1 y p2 ---
wide_ok <- requireNamespace("reshape2", quietly = TRUE)

if (wide_ok) {
  library(reshape2)

  iter_wide <- dcast(summary_long, model + tau + n ~ p, value.var = "mean_iter")
  key_cols <- c("model", "tau", "n")
  colnames(iter_wide)[!(colnames(iter_wide) %in% key_cols)] <-
    paste0("mean_iter_", colnames(iter_wide)[!(colnames(iter_wide) %in% key_cols)])

  time_wide <- dcast(summary_long, model + tau + n ~ p, value.var = "mean_time")
  colnames(time_wide)[!(colnames(time_wide) %in% key_cols)] <-
    paste0("mean_time_", colnames(time_wide)[!(colnames(time_wide) %in% key_cols)])

  summary_wide <- merge(iter_wide, time_wide, by = key_cols)
  summary_wide <- summary_wide[order(summary_wide$model, summary_wide$tau, summary_wide$n), ]
} else {
  # Fallback sin reshape2: armar manualmente
  keys <- unique(summary_long[, c("model", "tau", "n")])
  rows <- list()
  for (k in seq_len(nrow(keys))) {
    md <- keys$model[k]; tu <- keys$tau[k]; nn <- keys$n[k]
    sub <- summary_long[summary_long$model == md & summary_long$tau == tu & summary_long$n == nn, ]
    row <- data.frame(model = md, tau = tu, n = nn)
    for (pl in sort(unique(sub$p))) {
      row[[paste0("mean_iter_", pl)]] <- sub$mean_iter[sub$p == pl]
      row[[paste0("mean_time_", pl)]] <- sub$mean_time[sub$p == pl]
    }
    rows[[length(rows) + 1]] <- row
  }
  summary_wide <- do.call(rbind, rows)
  summary_wide <- summary_wide[order(summary_wide$model, summary_wide$tau, summary_wide$n), ]
}

cat("\n===== Resumen (formato ancho) =====\n")
print(summary_wide, row.names = FALSE)

## --- Guardar CSVs (en el directorio actual desde donde se corre el script) ---
write.csv(summary_long, file = "summary_iter_time_tau_long.csv", row.names = FALSE)
write.csv(summary_wide, file = "summary_iter_time_tau_wide.csv", row.names = FALSE)

cat("\nArchivos guardados:\n")
cat(" -", normalizePath("summary_iter_time_tau_long.csv"), "\n")
cat(" -", normalizePath("summary_iter_time_tau_wide.csv"), "\n")
