# =============================================================================
# 03_exploracion_inicial.R – QA (NA/duplicados), boxplots y PCA por cohorte
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# **Ejercutar primero 02_curacion_datos.R**
# =============================================================================

# En este script vamos a realizar la exploración inicial para el análisis. Aquí
# evaluamos la calidad de los metadatos; cargamos y emparejamos matrices de expresión
# con los metadatos y generamos diagnósticos por cohorte (gráficos):

# En primer lugar, cargamos las librerías:
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)})

path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")
invisible(lapply(c(path_tabs, path_figs), dir.create, recursive = TRUE, showWarnings = FALSE))


# Ahora, definimos funciones auxiliares para guardar las figuras, exportar las tablas y
# estandarizar los identificadores de las muestras:
save_plot <- function(p, filename, w = 9, h = 6, dpi = 300) {
  ggsave(filename = file.path(path_figs, filename), plot = p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))
}
write_tsv_safe <- function(df, filename) {
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))
}
to_upper_trim <- function(x) { x <- as.character(x); x <- trimws(x); toupper(x) }


# A partir de una tabla de expresión (genes en primera columna, muestras en resto),
# construimos una lista con la matriz numérica, el vector de genes y los IDs de muestra:
prepare_expr <- function(tbl, gene_col = 1) {
  stopifnot(is.data.frame(tbl))
  if (ncol(tbl) < 2) stop("La tabla de expresión no tiene columnas suficientes.")
  genes <- tbl[[gene_col]]
  mat <- as.matrix(tbl[,-gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes, samples = colnames(mat))}


# Emparejamos una matriz de expresión con "meta_master" usando "sample_id:
match_expr_meta <- function(expr_list, meta_df, cohort_tag) {
  X <- expr_list$expr
  cn <- colnames(X)
  meta_ids <- meta_df %>% dplyr::filter(cohort == cohort_tag) %>% dplyr::pull(sample_id) %>% unique()
  cn_up <- to_upper_trim(cn)
  keep <- cn_up %in% meta_ids
  if (!any(keep)) {
    cn2 <- gsub("\\.CEL$|\\.GZ$|\\.TXT$|\\.FASTQ$|\\.BAM$", "", cn_up)
    keep <- cn2 %in% meta_ids
    X <- X[, keep, drop = FALSE]
    colnames(X) <- cn2[keep]
  } else {
    X <- X[, keep, drop = FALSE]
    colnames(X) <- cn_up[keep]
  }
  list(expr = X, genes = expr_list$genes, samples = colnames(X))}


# Boxplot por muestra: visualizamos la distribución global (intensidades/recuentos):
plot_box <- function(X, title) {
  df <- as_tibble(X) %>%
    dplyr::mutate(.gene = dplyr::row_number()) %>%
    tidyr::pivot_longer(-.gene, names_to = "sample", values_to = "value")
  ggplot(df, aes(x = sample, y = value)) +
    geom_boxplot(outlier.size = 0.25) +
    labs(title = title, x = "Muestra", y = "Expresión") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_blank())}


# PCA por cohorte: proyectamos PC1–PC2 para evaluar patrones globales:
plot_pca <- function(X, meta_df, cohort_tag, title) {
  v <- apply(X, 1, stats::var)
  keep <- is.finite(v) & v > 0
  X2 <- X[keep, , drop = FALSE]
  if (ncol(X2) < 3) stop("Menos de 3 muestras tras emparejamiento; no se puede hacer PCA.")
  pcs <- prcomp(t(X2), center = TRUE, scale. = TRUE)
  pdat <- as_tibble(pcs$x[, 1:2]) %>%
    dplyr::mutate(sample_id = rownames(pcs$x)) %>%
    dplyr::left_join(meta_df %>% dplyr::filter(cohort == cohort_tag), by = dplyr::join_by(sample_id))
  ggplot(pdat, aes(x = PC1, y = PC2, color = response, shape = arm)) +
    geom_point(size = 2, alpha = 0.95) +
    labs(title = title, color = "Respuesta", shape = "Brazo") +
    theme_bw(base_size = 12)}

# Cargamos objetos del Script 2:
need <- c("meta_master", "ispy_expr", "medi_raw", "medi_norm", "htg_env", "htg_meta_h")
missing <- need[!need %in% ls(envir = .GlobalEnv)]
if (length(missing)) {
  message("Faltan en memoria: ", paste(missing, collapse = ", "))}
# Cargar meta_master desde TSV:
if (!exists("meta_master", envir = .GlobalEnv)) {
  meta_master_path <- file.path(path_tabs, "meta_master.tsv")
  if (file.exists(meta_master_path)) {
    meta_master <- readr::read_tsv(meta_master_path, show_col_types = FALSE) %>%
      dplyr::mutate(across(c(cohort, technology, sample_id, arm, response, timepoint),
                           ~ ifelse(is.na(.x), NA, as.character(.x))))
    message("meta_master cargado desde: ", meta_master_path)
  } else {stop("No se encontró 'meta_master.tsv'. Ejecutar antes script 2.")}}

# Reconstruir I-SPY2:
if (!exists("ispy_expr", envir = .GlobalEnv)) {
  path_ispy_expr <- file.path("data", "GSE173839",
                              "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz")
  if (file.exists(path_ispy_expr)) {
    ispy_expr <- data.table::fread(path_ispy_expr, sep = "\t", header = TRUE) |> as_tibble()
    message("I-SPY2 (Agilent) cargado desde disco. Dim: ", paste(dim(ispy_expr), collapse = " x "))
  } else {warning("No se encontró la matriz de I-SPY2 en disco.")}}

# Reconstruir MEDI:
if (!exists("medi_norm", envir = .GlobalEnv) && !exists("medi_raw", envir = .GlobalEnv)) {
  path_medi_norm <- file.path("data", "GSE241876", "GSE241876_DeseqNormalizedCount.csv.gz")
  path_medi_raw  <- file.path("data", "GSE241876", "GSE241876_raw_count_matrix.csv.gz")
  if (file.exists(path_medi_norm)) {
    medi_norm <- data.table::fread(path_medi_norm) |> as_tibble()
    message("MEDI4736 (DESeq2) cargado desde disco. Dim: ", paste(dim(medi_norm), collapse = " x "))
  } else if (file.exists(path_medi_raw)) {
    medi_raw <- data.table::fread(path_medi_raw) |> as_tibble()
    message("MEDI4736 (raw) cargado desde disco. Dim: ", paste(dim(medi_raw), collapse = " x "))
  } else {
    warning("No se encontró ninguna matriz de MEDI4736 en disco.")}}

# Cargar htg_env:
if (!exists("htg_env", envir = .GlobalEnv)) {
  path_htg_rdata <- file.path("data","G9_HTG","valid_datasets_HTG.RData")
  if (file.exists(path_htg_rdata)) {
    htg_env <- new.env()
    load(path_htg_rdata, envir = htg_env)
    message("✅ htg_env cargado desde: ", path_htg_rdata)}}

# Ahora, estandarizamos a tibble y calculamos duplicados por cohorte+sample_id
# (podrían indicar merges incorrectos), y número de NA por campo (guía para
# imputación o limpieza posterior):
meta_master <- tibble::as_tibble(meta_master)

# Duplicados por cohorte+sample_id:
dup <- meta_master %>%
  dplyr::group_by(cohort, sample_id) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1)
write_tsv_safe(dup, "H3_duplicados_meta_master.tsv")
print(dup, n = 20)

# Conteo de NA por campo:
na_tbl <- meta_master %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.)))) %>%
  tidyr::pivot_longer(dplyr::everything(), names_to = "campo", values_to = "n_na") %>%
  dplyr::arrange(dplyr::desc(n_na))
write_tsv_safe(na_tbl, "H3_NA_por_campo_meta_master.tsv")
print(na_tbl, n = nrow(na_tbl))

# Exploración I-SPY2 (Agilent). Preparamos la matriz, y luego generamos boxplot y PCA:
if (exists("ispy_expr") && is.data.frame(ispy_expr)) {
  ispy_prep <- prepare_expr(ispy_expr, gene_col = 1)
  ispy_matched <- match_expr_meta(ispy_prep, meta_master, "GSE173839_ISPY2")
  Xispy <- ispy_matched$expr
  message("I-SPY2 – muestras tras emparejar con meta: ", ncol(Xispy))
  
  med_all <- median(Xispy, na.rm = TRUE)
  if (is.finite(med_all) && med_all > 50) {
    Xispy_t <- log1p(Xispy); trans_used_ispy <- "log1p"
  } else {
    Xispy_t <- Xispy; trans_used_ispy <- "as-is"
  }
  
  p1 <- plot_box(Xispy_t, "I-SPY2 (Agilent) – Distribución por muestra")
  save_plot(p1, "H3_ISPY2_boxplot.png")
  
  p2 <- plot_pca(Xispy_t, meta_master, "GSE173839_ISPY2", "I-SPY2 (Agilent) – PCA")
  save_plot(p2, "H3_ISPY2_PCA.png")
  
  tibble(
    cohort = "GSE173839_ISPY2",
    n_genes = nrow(Xispy_t),
    n_samples = ncol(Xispy_t),
    transform = trans_used_ispy
  ) |> write_tsv_safe("H3_ISPY2_resumen.tsv")
} else {
  message("I-SPY2 no disponible en memoria/disco para exploración.")}

# Exploración MEDI4736. Eliminamos columnas de anotación para quedarnos
# únicamente con columnas de muestra:
medi_tbl <- if (exists("medi_norm") && is.data.frame(medi_norm)) medi_norm else {
  if (exists("medi_raw") && is.data.frame(medi_raw)) medi_raw else NULL
}

if (!is.null(medi_tbl)) {
  drop_annot_cols <- function(tbl) {
    bad <- c("ENSEMBLEID","ENTREZID","GENESYMBOL","SYMBOL","GENE","DESCRIPTION","V1")
    keep <- !(toupper(names(tbl)) %in% bad)
    tbl[, keep]
  }
  medi_tbl <- drop_annot_cols(medi_tbl)
  
  medi_prep <- prepare_expr(medi_tbl, gene_col = 1)
  medi_matched <- match_expr_meta(medi_prep, meta_master, "GSE241876_MEDI4736")
  Xmedi <- medi_matched$expr
  message("MEDI4736 – muestras tras emparejar con meta: ", ncol(Xmedi))
  
  Xmedi_t <- log1p(Xmedi)
  
  p3 <- plot_box(Xmedi_t, "MEDI4736 (RNA-seq) – Distribución por muestra")
  save_plot(p3, "H3_MEDI4736_boxplot.png")
  
  p4 <- plot_pca(Xmedi_t, meta_master, "GSE241876_MEDI4736", "MEDI4736 (RNA-seq) – PCA")
  save_plot(p4, "H3_MEDI4736_PCA.png")
  
  tibble(
    cohort = "GSE241876_MEDI4736",
    n_genes = nrow(Xmedi_t),
    n_samples = ncol(Xmedi_t),
    transform = "log1p"
  ) |> write_tsv_safe("H3_MEDI4736_resumen.tsv")
} else {
  message("MEDI4736 no disponible en memoria/disco para exploración.")
}

# Exploración HTG. Recorremos candidatos por nombre, aplicamos log1p y generamos boxplot y P:
if (exists("htg_env")) {
  objs <- ls(htg_env)
  expr_candidates <- objs[grepl("Rseq|expr|mat|counts", objs, ignore.case = TRUE)]
  did_any <- FALSE
  
  for (o in expr_candidates) {
    obj <- get(o, envir = htg_env)
    if (is.matrix(obj) || is.data.frame(obj)) {
      if (is.matrix(obj)) {
        X <- obj
        Xnum <- suppressWarnings(apply(X, 2, as.numeric))
        if (!is.matrix(Xnum)) Xnum <- as.matrix(X)
        colnames(Xnum) <- colnames(X)
        X_t <- log1p(abs(Xnum))
      } else {
        prep <- prepare_expr(as_tibble(obj), gene_col = 1)
        X_t <- log1p(abs(prep$expr))
      }
      
      tag <- gsub("[^A-Za-z0-9_]+", "_", o)
      
      p5 <- plot_box(X_t, paste0("HTG EdgeSeq – Distribución (", o, ")"))
      save_plot(p5, paste0("H3_HTG_", tag, "_boxplot.png"))
      
      v <- apply(X_t, 1, stats::var); keep <- is.finite(v) & v > 0
      if (sum(keep) > 10 && ncol(X_t) >= 3) {
        pcs <- prcomp(t(X_t[keep, , drop = FALSE]), center = TRUE, scale. = TRUE)
        pdat <- as_tibble(pcs$x[, 1:2]) %>% dplyr::mutate(sample = rownames(pcs$x))
        p6 <- ggplot(pdat, aes(PC1, PC2)) + geom_point(size = 2) +
          labs(title = paste0("HTG EdgeSeq – PCA (", o, ")")) + theme_bw(base_size = 12)
        save_plot(p6, paste0("H3_HTG_", tag, "_PCA.png"))
      }
      
      tibble(
        cohort = paste0("HTG_", toupper(tag)),
        n_genes = nrow(X_t),
        n_samples = ncol(X_t),
        transform = "log1p"
      ) |> write_tsv_safe(paste0("H3_", tag, "_resumen.tsv"))
      
      did_any <- TRUE}}
  
  if (!did_any) message("No se identificaron matrices de expresión utilizables dentro de htg_env.")
} else {
  message("htg_env no está disponible; omitimos la parte de HTG.")
}



message("Exploración inicial (QA) completado.")

