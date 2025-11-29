# =============================================================================
# 07_validacion_y_sensibilidad.R – Validación biológica, capacidad
# discriminativa y análisis de sensibilidad/estabilidad.
#
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
#
# Objetivos (H7):
#   1) Validación biológica general en I-SPY2:
#        - Coherencia de señales entre vías Hallmark inmunes y marcadores
#          de respuesta inmune (PDCD1, CD274, CD8A, GZMB, CXCL9, CXCL10).
#   2) Evaluar capacidad discriminativa (separación R/NR) y estabilidad:
#        - Modelo logístico usando el score de HALLMARK_INTERFERON_GAMMA_RESPONSE.
#        - AUC global y distribución de AUC bajo submuestreos balanceados R/NR.
#   3) Priorización final de genes y vías con evidencia convergente
#      entre cohortes (usando tablas del script 06).
# =============================================================================

#### 0) Carga de librerías y rutas ###########################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(data.table)
  library(ggplot2)
  library(pROC)      # para ROC/AUC
})

path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")
invisible(lapply(c(path_tabs, path_figs), dir.create,
                 recursive = TRUE, showWarnings = FALSE))


#### Funciones auxiliares ####################################################

write_tsv_safe <- function(df, filename) {
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))
}

save_plot <- function(p, filename, w = 7, h = 6, dpi = 300) {
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))
}

to_upper_trim <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)
}

# Prepara lista a partir de una tabla de expresión: primera col = gen
prepare_expr <- function(tbl, gene_col = 1) {
  stopifnot(is.data.frame(tbl))
  if (ncol(tbl) < 2) stop("La tabla de expresión no tiene columnas suficientes.")
  genes <- tbl[[gene_col]]
  mat   <- as.matrix(tbl[, -gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes)
}

# Empareja matriz de expresión con meta_master para una cohorte dada
match_expr_meta <- function(expr_list, meta_df, cohort_tag) {
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
    X <- X[, keep, drop = FALSE]
    colnames(X) <- cn_up[keep]
  }
  
  samples <- colnames(X)
  
  pheno <- meta_df %>%
    dplyr::filter(cohort == cohort_tag, sample_id %in% samples) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  pheno <- pheno[match(samples, pheno$sample_id), , drop = FALSE]
  stopifnot(identical(pheno$sample_id, samples))
  
  list(expr = X, genes = expr_list$genes, pheno = pheno)
}

# Colapsa genes duplicados usando la media de expresión
collapse_by_gene <- function(X, genes) {
  ord   <- order(genes)
  genes <- genes[ord]
  X     <- X[ord, , drop = FALSE]
  X2    <- rowsum(X, group = genes) / as.vector(table(genes))
  list(expr = X2, genes = rownames(X2))
}

# Scatter genérico con cuadrantes y color por categoría
plot_scatter_comparison <- function(df,
                                    x_var, y_var,
                                    col_var,
                                    x_lab, y_lab,
                                    title) {
  ggplot(df, aes(x = !!sym(x_var),
                 y = !!sym(y_var),
                 color = !!sym(col_var))) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = title,
      x = x_lab,
      y = y_lab,
      color = "Categoría"
    ) +
    geom_line(linewidth = 1)
  
}


# =============================================================================
# PARTE A – Validación biológica en I-SPY2
#   - Reconstruimos expresión gen a gen y scores ssGSEA Hallmark
#   - Correlaciones entre vías inmunes y marcadores relacionados
# =============================================================================

message("=== PARTE A: Validación biológica en I-SPY2 ===")

## 1) Cargar meta_master e I-SPY2 (igual que en scripts previos)  -------------

meta_master <- readr::read_tsv(
  file.path(path_tabs, "meta_master.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(across(
    c(cohort, technology, sample_id, arm, response, timepoint),
    as.character
  ))

ispy_tbl <- data.table::fread(
  file.path(
    "data", "GSE173839",
    "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"
  ),
  sep = "\t", header = TRUE
) %>% tibble::as_tibble()

message("I-SPY2 cargado. Dim: ", paste(dim(ispy_tbl), collapse = " x "))

## 2) Preparar expresión gen a gen y pheno (R/NR) -----------------------------

prep_ispy <- prepare_expr(ispy_tbl, gene_col = 1)
X_all     <- prep_ispy$expr
genes_all <- prep_ispy$genes

# Filtro mínimo por varianza
v_all <- apply(X_all, 1, stats::var)
keep  <- is.finite(v_all) & v_all > 0
X     <- X_all[keep, , drop = FALSE]
genes <- genes_all[keep]

matched_ispy <- match_expr_meta(
  list(expr = X, genes = genes),
  meta_master,
  "GSE173839_ISPY2"
)

X          <- matched_ispy$expr
genes      <- matched_ispy$genes
pheno_ispy <- matched_ispy$pheno

pheno_ispy <- pheno_ispy %>%
  dplyr::filter(response %in% c("R", "NR")) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR", "R")),
    arm      = to_upper_trim(arm),
    arm      = forcats::fct_lump_min(arm, min = 5),
    arm      = stats::relevel(factor(arm), ref = "CONTROL")
  )

samples_keep_ispy <- pheno_ispy$sample_id
X <- X[, samples_keep_ispy, drop = FALSE]
pheno_ispy <- pheno_ispy[match(colnames(X), pheno_ispy$sample_id), , drop = FALSE]
stopifnot(identical(pheno_ispy$sample_id, colnames(X)))

# Colapsamos genes duplicados → matriz genes x muestras
collapsed_ispy <- collapse_by_gene(X, genes)
X_gene_ispy    <- collapsed_ispy$expr
rownames(X_gene_ispy) <- collapsed_ispy$genes

message("I-SPY2 – genes tras colapso: ", nrow(X_gene_ispy))
message("I-SPY2 – muestras R/NR:"); print(table(pheno_ispy$response))


## 3) Cargar scores ssGSEA Hallmark de I-SPY2 --------------------------------

scores_ispy_df <- readr::read_tsv(
  file.path(path_tabs, "H5_ISPY2_ssGSEA_scores.tsv"),
  show_col_types = FALSE
)

# Primera columna = pathway, resto = muestras
stopifnot("pathway" %in% names(scores_ispy_df))

scores_mat <- scores_ispy_df %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# Reordenamos columnas para casarlas con pheno_ispy
common_samples <- intersect(colnames(scores_mat), pheno_ispy$sample_id)
scores_mat     <- scores_mat[, common_samples, drop = FALSE]
pheno_ispy     <- pheno_ispy[match(colnames(scores_mat), pheno_ispy$sample_id), ]
stopifnot(identical(colnames(scores_mat), pheno_ispy$sample_id))

message("I-SPY2 – dim scores ssGSEA (vías x muestras): ",
        paste(dim(scores_mat), collapse = " x "))


## 4) Correlaciones entre vías inmunes y marcadores de respuesta  -------------

# Seleccionamos vías Hallmark inmunes relevantes (si existen en el objeto):
pathways_interest_all <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"
)
pathways_interest <- intersect(pathways_interest_all, rownames(scores_mat))
message("Vías inmunes disponibles en ssGSEA:"); print(pathways_interest)

# Marcadores inmunes de interés
marker_genes_all <- c("PDCD1", "CD274", "CTLA4", "CD8A", "GZMB", "CXCL9", "CXCL10")
marker_genes     <- intersect(marker_genes_all, rownames(X_gene_ispy))
message("Marcadores presentes en expresión I-SPY2:"); print(marker_genes)

# Construimos tabla larga gen-vía con correlaciones Spearman
corr_list <- list()

for (g in marker_genes) {
  expr_g <- as.numeric(X_gene_ispy[g, colnames(scores_mat)])
  
  for (p in pathways_interest) {
    score_p <- as.numeric(scores_mat[p, ])
    
    cor_val <- suppressWarnings(
      cor(expr_g, score_p, method = "spearman", use = "complete.obs")
    )
    
    corr_list[[length(corr_list) + 1]] <- tibble::tibble(
      gene    = g,
      pathway = p,
      cor_spearman = cor_val
    )
  }
}

corr_table <- dplyr::bind_rows(corr_list)

write_tsv_safe(corr_table, "H7_ISPY2_corr_markers_vs_pathways.tsv")
message("Resumen correlaciones marcadores–vías:")
print(corr_table)

# Figura de ejemplo: scatter para una pareja gene–vía (si existe)
if (length(marker_genes) > 0 && length(pathways_interest) > 0) {
  gene_plot    <- marker_genes[1]
  pathway_plot <- pathways_interest[1]
  
  df_scatter <- tibble::tibble(
    sample_id  = colnames(scores_mat),
    response   = pheno_ispy$response,
    expr_gene  = as.numeric(X_gene_ispy[gene_plot, colnames(scores_mat)]),
    score_path = as.numeric(scores_mat[pathway_plot, ])
  )
  
  p_corr <- ggplot(df_scatter,
                   aes(x = expr_gene, y = score_path, color = response)) +
    geom_point(size = 2, alpha = 0.8) +
    labs(
      title = paste0("I-SPY2 – ", gene_plot, " vs ", pathway_plot),
      x     = paste0("Expresión de ", gene_plot, " (log2 escala Agilent)"),
      y     = paste0("Score ssGSEA: ", pathway_plot),
      color = "Respuesta"
    ) +
    theme_bw(base_size = 11)
  
  save_plot(p_corr, "H7_ISPY2_scatter_ejemplo_gene_vs_pathway.png")
}

message("PARTE A COMPLETADA.\n")


# =============================================================================
# PARTE B – Capacidad discriminativa y análisis de estabilidad
#   - Firma simple basada en HALLMARK_INTERFERON_GAMMA_RESPONSE
#   - Modelo logístico R vs NR y AUC
#   - Submuestreos balanceados R/NR para evaluar estabilidad
# =============================================================================

message("=== PARTE B: Capacidad discriminativa y estabilidad (I-SPY2) ===")

## 1) Construimos data.frame con score IFN-γ y respuesta ----------------------

path_ifng <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"

if (!(path_ifng %in% rownames(scores_mat))) {
  stop("La vía HALLMARK_INTERFERON_GAMMA_RESPONSE no está disponible en scores_mat.")
}

df_model <- tibble::tibble(
  sample_id      = colnames(scores_mat),
  response       = pheno_ispy$response,
  arm            = pheno_ispy$arm,
  IFNG_score     = as.numeric(scores_mat[path_ifng, ])
)

# Comprobación rápida de grupos
message("Distribución respuesta en I-SPY2 para modelo IFNG:")
print(table(df_model$response))

## 2) Modelo logístico global (ajustado por brazo) ---------------------------

fit_ifng <- glm(response ~ IFNG_score + arm,
                data = df_model,
                family = binomial)

# Probabilidades predichas (probabilidad de R)
df_model$prob_R <- predict(fit_ifng, type = "response")

# Curva ROC y AUC
roc_ifng <- pROC::roc(
  response  = df_model$response,
  predictor = df_model$prob_R,
  levels    = c("NR", "R"),
  direction = "<"
)

auc_ifng  <- as.numeric(pROC::auc(roc_ifng))
ci_ifng   <- as.numeric(pROC::ci.auc(roc_ifng))

summary_auc <- tibble::tibble(
  cohort = "GSE173839_ISPY2",
  marker = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  AUC    = auc_ifng,
  CI_low = ci_ifng[1],
  CI_high= ci_ifng[3]
)

write_tsv_safe(summary_auc, "H7_ISPY2_IFNG_AUC_global.tsv")
message("AUC global modelo IFNG (I-SPY2):")
print(summary_auc)

# Boxplot del score por respuesta
p_box_ifng <- ggplot(df_model,
                     aes(x = response, y = IFNG_score, fill = response)) +
  geom_boxplot(alpha = 0.8) +
  labs(
    title = "I-SPY2 – Score ssGSEA IFN-γ vs respuesta",
    x     = "Respuesta clínica",
    y     = "Score ssGSEA HALLMARK_INTERFERON_GAMMA_RESPONSE"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")

save_plot(p_box_ifng, "H7_ISPY2_IFNG_score_boxplot_R_vs_NR.png")

# Curva ROC
roc_df <- tibble::tibble(
  sens = roc_ifng$sensitivities,
  spec = roc_ifng$specificities
)

p_roc_ifng <- ggplot(roc_df, aes(x = 1 - spec, y = sens)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    title = paste0("I-SPY2 – ROC IFN-γ (AUC = ",
                   round(auc_ifng, 3), ")"),
    x = "1 - Especificidad",
    y = "Sensibilidad"
  ) +
  theme_bw(base_size = 11)

save_plot(p_roc_ifng, "H7_ISPY2_IFNG_ROC_curve.png")


## 3) Submuestreos balanceados R/NR para evaluar estabilidad -----------------

set.seed(123)

n_R  <- sum(df_model$response == "R")
n_NR <- sum(df_model$response == "NR")
n_min <- min(n_R, n_NR)

message("Tamaños de grupo en I-SPY2: R = ", n_R, ", NR = ", n_NR,
        " → tamaño balanceado por submuestreo = ", n_min)

B_iter <- 200   # número de submuestreos

get_auc_subsample <- function(b, df) {
  set.seed(1000 + b)
  
  idx_R  <- which(df$response == "R")
  idx_NR <- which(df$response == "NR")
  
  idx_R_sub  <- sample(idx_R,  n_min, replace = FALSE)
  idx_NR_sub <- sample(idx_NR, n_min, replace = FALSE)
  
  idx_sub <- c(idx_R_sub, idx_NR_sub)
  
  df_b <- df[idx_sub, , drop = FALSE]
  
  fit_b <- glm(response ~ IFNG_score,
               data = df_b,
               family = binomial)
  
  prob_R_b <- predict(fit_b, type = "response")
  
  roc_b <- pROC::roc(
    response  = df_b$response,
    predictor = prob_R_b,
    levels    = c("NR", "R"),
    direction = "<"
  )
  
  as.numeric(pROC::auc(roc_b))
}

auc_sub_vec <- sapply(1:B_iter, get_auc_subsample, df = df_model)

auc_sub_df <- tibble::tibble(
  iter = 1:B_iter,
  AUC  = auc_sub_vec
)

write_tsv_safe(auc_sub_df, "H7_ISPY2_IFNG_AUC_submuestreos.tsv")

summary_sub <- tibble::tibble(
  cohort = "GSE173839_ISPY2",
  marker = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  n_iter = B_iter,
  AUC_median = median(auc_sub_vec, na.rm = TRUE),
  AUC_IQR_low = quantile(auc_sub_vec, 0.25, na.rm = TRUE),
  AUC_IQR_high = quantile(auc_sub_vec, 0.75, na.rm = TRUE),
  AUC_min = min(auc_sub_vec, na.rm = TRUE),
  AUC_max = max(auc_sub_vec, na.rm = TRUE)
)

write_tsv_safe(summary_sub, "H7_ISPY2_IFNG_AUC_submuestreos_resumen.tsv")
message("Distribución de AUC en submuestreos balanceados:")
print(summary_sub)

p_hist_auc <- ggplot(auc_sub_df, aes(x = AUC)) +
  geom_histogram(bins = 30, color = "black", fill = "grey80") +
  geom_vline(xintercept = auc_ifng, linetype = "dashed") +
  labs(
    title = "I-SPY2 – Distribución de AUC\nsubmuestreos balanceados R/NR (IFN-γ)",
    x = "AUC (modelo IFN-γ, R vs NR)",
    y = "Frecuencia"
  ) +
  theme_bw(base_size = 11)

save_plot(p_hist_auc, "H7_ISPY2_IFNG_AUC_submuestreos_histograma.png")

message("PARTE B COMPLETADA.\n")


# =============================================================================
# PARTE C – Priorización final de vías y genes con evidencia convergente
#   (basado en tablas generadas en el script 06)
# =============================================================================

message("=== PARTE C: Priorización final genes/vías (evidencia convergente) ===")

## 1) Vías Hallmark: combinamos información de H6 -----------------------------

pathways_merged <- readr::read_tsv(
  file.path(path_tabs, "H6_pathways_ISPY2_vs_GSE241876_merged.tsv"),
  show_col_types = FALSE
)

# Añadimos criterio de concordancia de dirección:
# - logFC_ISPY2 > 0 y NES_GSE > 0  (activadas en R en ambas)
# - logFC_ISPY2 < 0 y NES_GSE < 0  (reprimidas en R en ambas)
pathways_prior <- pathways_merged %>%
  dplyr::mutate(
    same_direction = sign(logFC_ISPY2) == sign(NES_GSE),
    priorizada = (FDR_ISPY2 < 0.05 | FDR_GSE < 0.05) & same_direction
  ) %>%
  dplyr::arrange(FDR_ISPY2, FDR_GSE)

write_tsv_safe(pathways_prior, "H7_pathways_priorizadas.tsv")

message("Resumen vías priorizadas (coherencia entre cohortes):")
print(
  pathways_prior %>%
    dplyr::filter(priorizada) %>%
    dplyr::select(pathway, logFC_ISPY2, NES_GSE, FDR_ISPY2, FDR_GSE) %>%
    dplyr::slice_head(n = 20)
)

## 2) Genes: usamos tabla de H6 y pedimos
##    - FDR_ISPY2 < 0.05 (descubrimiento fuerte)
##    - misma dirección de efecto en ambas cohortes
##      (aunque en GSE241876 no sea significativo)
## ---------------------------------------------------------------------------

genes_merged <- readr::read_tsv(
  file.path(path_tabs, "H6_genes_ISPY2_vs_GSE241876_merged.tsv"),
  show_col_types = FALSE
)

genes_prior <- genes_merged %>%
  dplyr::mutate(
    same_direction = sign(logFC_ISPY2) == sign(logFC_GSE),
    priorizado = (FDR_ISPY2 < 0.05) & same_direction
  ) %>%
  dplyr::arrange(FDR_ISPY2)

write_tsv_safe(genes_prior, "H7_genes_priorizados.tsv")

message("Top genes priorizados (FDR_ISPY2 < 0.05 y misma dirección):")
print(
  genes_prior %>%
    dplyr::filter(priorizado) %>%
    dplyr::select(Gene, logFC_ISPY2, logFC_GSE, FDR_ISPY2, FDR_GSE) %>%
    dplyr::slice_head(n = 30)
)

message("=== Script 07 (validación + sensibilidad + priorización final) COMPLETADO ===")

