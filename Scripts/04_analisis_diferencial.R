# =============================================================================
# 04_analisis_diferencial.R – R vs NR por cohorte
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================

# En este script desarrollamos el análisis diferencial para la cohorte I-SPY2,
# contrastando la expresión génica entre pacientes respondedores (R) y no
# respondedores (NR), y controlando el posible efecto del brazo de tratamiento.

# Cargamos las librerías:
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(limma)
  library(edgeR)})

# Definimos las rutas principales de trabajo:
path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")
invisible(lapply(c(path_tabs, path_figs), dir.create, recursive = TRUE, showWarnings = FALSE))

# Funciones de apoyo:
# Función para limpiar y homogeneizar identificadores:
to_upper_trim <- function(x){ x <- as.character(x); x <- trimws(x); toupper(x) }

# Función para guardar tablas en formato TSV:
write_tsv_safe <- function(df, filename){
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))}

# Función para guardar gráficos en PNG:
save_plot <- function(p, filename, w = 8, h = 6, dpi = 300){
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))}


# Preparamos lista con matriz de expresión y vectores de genes/muestras
# (primera columna = genes. Resto = muestras):
prepare_expr <- function(tbl, gene_col = 1){
  stopifnot(is.data.frame(tbl))
  genes <- tbl[[gene_col]]
  mat <- as.matrix(tbl[, -gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes)}

# Alineamos la matriz de expresión con los metadatos del objeto meta_master:
match_expr_meta <- function(expr_list, meta_df, cohort_tag){
  X  <- expr_list$expr
  cn <- colnames(X)
  meta_ids <- meta_df %>%
    dplyr::filter(cohort == cohort_tag) %>%
    dplyr::pull(sample_id) %>%
    unique()
  cn_up <- to_upper_trim(cn)
  keep  <- cn_up %in% meta_ids
  if (!any(keep)) {
    cn2  <- gsub("\\.CEL$|\\.GZ$|\\.TXT$|\\.FASTQ$|\\.BAM$", "", cn_up)
    keep <- cn2 %in% meta_ids
    X    <- X[, keep, drop = FALSE]
    colnames(X) <- cn2[keep]
  } else {
    X    <- X[, keep, drop = FALSE]
    colnames(X) <- cn_up[keep]
  }
  X
}

# Función para generar un gráfico tipo “volcano plot” --> Evaluar la magnitud (logFC)
# frente a la significación (FDR ajustado):
volcano_plot <- function(tt, title, lfc_thr = log2(1.5), p_thr = 0.05){
  tt2 <- tt %>%
    dplyr::mutate(sig = ifelse(adj.P.Val < p_thr & abs(logFC) >= lfc_thr,
                               "Significativo", "No significativo"))
  ggplot(tt2, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(shape = sig), alpha = 0.85) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
    geom_hline(yintercept = -log10(p_thr), linetype = 2) +
    labs(title = title, x = "log2(FC)", y = "-log10(FDR)") +
    theme_bw(base_size = 12)
}


# Función para generar un gráfico MA (cambio de expresión en función de la media global):
ma_plot <- function(tt, title){
  ggplot(tt, aes(x = AveExpr, y = logFC)) +
    geom_point(alpha = .6) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = title, x = "Media log2 expresión", y = "log2(FC) R vs NR") +
    theme_bw(base_size = 12)
}


# A continuación, cargamos la tabla de metadatos armonizados y la matriz de expresión de la cohorte
# I-SPY2 generados en el script 2:
meta_master <- readr::read_tsv(file.path(path_tabs, "meta_master.tsv"),
                               show_col_types = FALSE) %>%
  dplyr::mutate(across(c(cohort, technology, sample_id, arm, response, timepoint),
                       as.character))

ispy_tbl <- data.table::fread(file.path("data","GSE173839",
                                        "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"),
                              sep = "\t", header = TRUE) %>% tibble::as_tibble()

# Ahora, realizamos el contraste R vs NR ajustando por el brazo de tratamiento:
message("I-SPY2 – preparación del diseño (R vs NR, covariable 'arm').")

# 1. Preparamos expresión y filtramos genes de varianza > 0:
prep   <- prepare_expr(ispy_tbl, gene_col = 1)
X_all  <- prep$expr
genes  <- prep$genes
v_all  <- apply(X_all, 1, stats::var)
keep   <- is.finite(v_all) & v_all > 0
X      <- X_all[keep, , drop = FALSE]
genes  <- genes[keep]

# 2. Emparejamos columnas con metadatos (cohorte I-SPY2):
X <- match_expr_meta(list(expr = X, genes = genes), meta_master, "GSE173839_ISPY2")
samples <- colnames(X)

# 3. Construimos fenodata y ordenamos para que coincida con X (sin usar slice):
pheno <- meta_master %>%
  dplyr::filter(cohort == "GSE173839_ISPY2", sample_id %in% samples) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR","R")),   # coeficiente 'responseR' = R>NR
    arm      = to_upper_trim(arm),
    arm      = forcats::fct_lump_min(arm, min = 5),      # agrupar brazos raros
    arm      = stats::relevel(factor(arm), ref = "CONTROL")
  )
# reordenamos filas para que sample_id == columnas de X
pheno <- pheno[match(samples, pheno$sample_id), , drop = FALSE]
stopifnot(identical(pheno$sample_id, samples))

# 4. Ajustamos limma. Incluimos la variable "response" (efecto principal) y "arm" (covariable):
design <- stats::model.matrix(~ response + arm, data = pheno)
fit <- limma::lmFit(X, design)
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

# Extraemos resultados para el contraste R vs NR:
coef_name <- "responseR"
if (!(coef_name %in% colnames(coef(fit)))) {
  stop("No se encontró el coeficiente 'responseR' en el modelo")}

tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("rowid") %>%
  dplyr::mutate(Gene = genes[as.integer(rowid)]) %>%
  dplyr::select(Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# 5. Guardamos tablas y figuras:
write_tsv_safe(tt, "H4_ISPY2_limma_R_vs_NR.tsv")
write_tsv_safe(tt %>% dplyr::slice_head(n = 50), "H4_ISPY2_limma_R_vs_NR_top50.tsv")

save_plot(volcano_plot(tt, "I-SPY2 – R vs NR (limma, ~ response + arm)"),
          "H4_ISPY2_volcano.png")
save_plot(ma_plot(tt, "I-SPY2 – MA plot (limma)"),
          "H4_ISPY2_MA.png")

# 6. Generamos un resumen con el número de genes diferenciales significativos:
summary_ispy <- tibble::tibble(
  cohort = "GSE173839_ISPY2",
  n_genes = nrow(tt),
  n_sig_fdr_05 = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
  n_sig_absLFC_1.5_FDR05 = sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= log2(1.5), na.rm = TRUE)
)
write_tsv_safe(summary_ispy, "H4_ISPY2_resumen.tsv")
print(summary_ispy)


message("Análisis diferencial para el I-SPY2 completado.")

