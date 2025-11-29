# =============================================================================
# 06_comparacion_cohortes.R – Integración I-SPY2 vs GSE241876 (MEDI)
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# En este script comparamos los resultados obtenidos en:
#
#   - I-SPY2 (GSE173839) – microarrays Agilent, ssGSEA Hallmark + limma
#   - GSE241876 (MEDI)    – RNA-seq, limma-voom a nivel de gen + GSEA (fgsea)
#
# PARTE A: Comparación a nivel de vías Hallmark (ssGSEA vs fgsea)
#   1) Cargar resultados de vías de I-SPY2 (H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv)
#   2) Cargar resultados de GSE241876 (H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv)
#   3) Unir por pathway, calcular estadísticos y correlación
#   4) Figura: diagrama de dispersión logFC_ISPY2 vs NES_GSE
#
# PARTE B: Comparación a nivel de gen (limma R vs NR)
#   5) Cargar limma I-SPY2 (archivo H?_ISPY2_limma_R_vs_NR.tsv)
#   6) Cargar limma GSE241876 (H5_GSE241876_limma_R_vs_NR.tsv)
#   7) Unir por gen, definir categorías de significación
#   8) Calcular correlaciones de logFC y generar figura de dispersión
#
# Todas las tablas se guardan en results/tablas y las figuras en results/figuras.
# =============================================================================


#### 0) Librerías y rutas ####

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
})

path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")

invisible(lapply(c(path_tabs, path_figs), dir.create,
                 recursive = TRUE, showWarnings = FALSE))


#### 1) Funciones auxiliares ####

write_tsv_safe <- function(df, filename) {
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))
}

save_plot <- function(p, filename, w = 7, h = 6, dpi = 300) {
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))
}

# Scatter genérico para comparación entre cohortes
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
    theme_bw(base_size = 11)
}


# =============================================================================
# PARTE A – Comparación a nivel de vías Hallmark
# =============================================================================

message("=== PARTE A: Comparación de vías Hallmark (I-SPY2 vs GSE241876) ===")

# 2) Cargar resultados de vías de I-SPY2 (ssGSEA + limma)
#    Archivo generado en el script 5:
#    H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv (con columna 'pathway')
ispy2_path <- readr::read_tsv(
  file.path(path_tabs, "H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv"),
  show_col_types = FALSE
)

# Aseguramos columna 'pathway'
if (!("pathway" %in% names(ispy2_path))) {
  stop("El archivo H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv no tiene columna 'pathway'.")
}

ispy2_clean <- ispy2_path %>%
  dplyr::select(pathway, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_ISPY2 = logFC,
    FDR_ISPY2   = adj.P.Val
  )

# 3) Cargar resultados de vías de GSE241876 (fgsea Hallmark)
#    Archivo generado en el script 5:
#    H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv
gse_path <- readr::read_tsv(
  file.path(path_tabs, "H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv"),
  show_col_types = FALSE
)

# Aseguramos columna 'pathway'
if (!("pathway" %in% names(gse_path))) {
  stop("El archivo H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv no tiene columna 'pathway'.")
}

gse_clean <- gse_path %>%
  dplyr::select(pathway, NES, padj) %>%
  dplyr::rename(
    NES_GSE = NES,
    FDR_GSE = padj
  )

# 4) Unimos por vía Hallmark
pathways_merged <- dplyr::inner_join(
  ispy2_clean,
  gse_clean,
  by = "pathway"
)

message("Vías comunes Hallmark: ", nrow(pathways_merged))

# Definimos variables lógicas de significación
pathways_merged <- pathways_merged %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_GSE   = FDR_GSE   < 0.05,
    sig_both  = sig_ISPY2 & sig_GSE,
    categoria = dplyr::case_when(
      sig_both ~ "Significativa en ambas",
      sig_ISPY2 & !sig_GSE ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_GSE ~ "Solo GSE241876",
      TRUE ~ "No significativa"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativa en ambas",
        "Solo I-SPY2",
        "Solo GSE241876",
        "No significativa"
      )
    )
  )

# Estadísticos de correlación
cor_paths_all <- suppressWarnings(
  cor(pathways_merged$logFC_ISPY2,
      pathways_merged$NES_GSE,
      method = "spearman",
      use = "complete.obs")
)

summary_paths <- tibble::tibble(
  n_pathways_comunes = nrow(pathways_merged),
  n_sig_ISPY2        = sum(pathways_merged$sig_ISPY2, na.rm = TRUE),
  n_sig_GSE          = sum(pathways_merged$sig_GSE,   na.rm = TRUE),
  n_sig_ambas        = sum(pathways_merged$sig_both,  na.rm = TRUE),
  cor_spearman_logFC_ISPY2_vs_NES_GSE = cor_paths_all
)

write_tsv_safe(pathways_merged, "H6_pathways_ISPY2_vs_GSE241876_merged.tsv")
write_tsv_safe(summary_paths,   "H6_pathways_ISPY2_vs_GSE241876_resumen.tsv")

message("Resumen vías (I-SPY2 vs GSE241876):")
print(summary_paths)

# 5) Figura: diagrama de dispersión logFC_ISPY2 vs NES_GSE
p_paths <- plot_scatter_comparison(
  df     = pathways_merged,
  x_var  = "logFC_ISPY2",
  y_var  = "NES_GSE",
  col_var = "categoria",
  x_lab  = "I-SPY2 – log2(FC) (R vs NR) en ssGSEA (Hallmark)",
  y_lab  = "GSE241876 – NES (fgsea, Hallmark)",
  title  = "Comparación de vías Hallmark entre I-SPY2 y GSE241876"
)

save_plot(p_paths, "H6_ISPY2_vs_GSE241876_Hallmark_scatter.png")

message("PARTE A COMPLETADA.")


# =============================================================================
# PARTE B – Comparación a nivel de gen (limma)
# =============================================================================

message("=== PARTE B: Comparación a nivel de gen (I-SPY2 vs GSE241876) ===")

#### 6) Cargar resultados limma a nivel de gen ####

# a) I-SPY2: archivo de limma R vs NR
#    En tu pipeline previo seguramente se llamó algo como:
#      H3_ISPY2_limma_R_vs_NR.tsv   o   H4_ISPY2_limma_R_vs_NR.tsv
#
#    Aquí intentamos primero con H3 y, si no existe, con H4.
file_ispy_genes_H3 <- file.path(path_tabs, "H3_ISPY2_limma_R_vs_NR.tsv")
file_ispy_genes_H4 <- file.path(path_tabs, "H4_ISPY2_limma_R_vs_NR.tsv")

if (file.exists(file_ispy_genes_H3)) {
  file_ispy_genes <- file_ispy_genes_H3
} else if (file.exists(file_ispy_genes_H4)) {
  file_ispy_genes <- file_ispy_genes_H4
} else {
  stop("No se encontró ni H3_ISPY2_limma_R_vs_NR.tsv ni H4_ISPY2_limma_R_vs_NR.tsv en results/tablas.")
}

message("Usando resultados de I-SPY2 en: ", basename(file_ispy_genes))

ispy2_genes <- readr::read_tsv(
  file_ispy_genes,
  show_col_types = FALSE
)

# Aseguramos columna 'Gene' (símbolo de gen)
if (!("Gene" %in% names(ispy2_genes))) {
  # Si no existe, intentamos sacarla del rowname:
  ispy2_genes <- ispy2_genes %>%
    tibble::rownames_to_column("Gene")
}

ispy2_genes_clean <- ispy2_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_ISPY2 = logFC,
    FDR_ISPY2   = adj.P.Val
  )

# b) GSE241876: archivo limma generado en el script 5
#    H5_GSE241876_limma_R_vs_NR.tsv
gse_genes <- readr::read_tsv(
  file.path(path_tabs, "H5_GSE241876_limma_R_vs_NR.tsv"),
  show_col_types = FALSE
)

if (!("Gene" %in% names(gse_genes))) {
  gse_genes <- gse_genes %>%
    tibble::rownames_to_column("Gene")
}

gse_genes_clean <- gse_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_GSE = logFC,
    FDR_GSE   = adj.P.Val
  )

#### 7) Unimos por símbolo de gen ####

genes_merged <- dplyr::inner_join(
  ispy2_genes_clean,
  gse_genes_clean,
  by = "Gene"
)

message("Genes comunes entre I-SPY2 y GSE241876: ", nrow(genes_merged))

# Definimos significación y categorías
genes_merged <- genes_merged %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_GSE   = FDR_GSE   < 0.05,
    sig_both  = sig_ISPY2 & sig_GSE,
    categoria = dplyr::case_when(
      sig_both ~ "Significativo en ambas",
      sig_ISPY2 & !sig_GSE ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_GSE ~ "Solo GSE241876",
      TRUE ~ "No significativo"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativo en ambas",
        "Solo I-SPY2",
        "Solo GSE241876",
        "No significativo"
      )
    )
  )

#### 8) Estadísticos de correlación ####

# Correlación global de logFC
cor_genes_all <- suppressWarnings(
  cor(genes_merged$logFC_ISPY2,
      genes_merged$logFC_GSE,
      method = "spearman",
      use   = "complete.obs")
)

# Correlación solo en genes significativos en ambas (si existieran)
if (any(genes_merged$sig_both, na.rm = TRUE)) {
  df_sig_both <- genes_merged %>%
    dplyr::filter(sig_both) %>%
    dplyr::filter(
      is.finite(logFC_ISPY2),
      is.finite(logFC_GSE)
    )
  
  if (nrow(df_sig_both) > 1) {
    cor_genes_sig_both <- suppressWarnings(
      cor(df_sig_both$logFC_ISPY2,
          df_sig_both$logFC_GSE,
          method = "spearman",
          use   = "complete.obs")
    )
  } else {
    cor_genes_sig_both <- NA_real_
  }
} else {
  cor_genes_sig_both <- NA_real_
}

# Resumen numérico
summary_genes <- tibble::tibble(
  n_genes_comunes          = nrow(genes_merged),
  n_sig_ISPY2              = sum(genes_merged$sig_ISPY2, na.rm = TRUE),
  n_sig_GSE                = sum(genes_merged$sig_GSE,   na.rm = TRUE),
  n_sig_ambas              = sum(genes_merged$sig_both,  na.rm = TRUE),
  cor_spearman_logFC_global = cor_genes_all,
  cor_spearman_logFC_sig_both = cor_genes_sig_both
)

write_tsv_safe(genes_merged, "H6_genes_ISPY2_vs_GSE241876_merged.tsv")
write_tsv_safe(summary_genes, "H6_genes_ISPY2_vs_GSE241876_resumen.tsv")

message("Resumen genes (I-SPY2 vs GSE241876):")
print(summary_genes)

#### 9) Figura: diagrama de dispersión de logFC por gen ####

# Nos quedamos solo con genes con logFC finito en ambas cohortes
genes_plot <- genes_merged %>%
  dplyr::filter(
    is.finite(logFC_ISPY2),
    is.finite(logFC_GSE)
  )

p_genes <- plot_scatter_comparison(
  df     = genes_plot,
  x_var  = "logFC_ISPY2",
  y_var  = "logFC_GSE",
  col_var = "categoria",
  x_lab  = "I-SPY2 – log2(FC) (R vs NR)",
  y_lab  = "GSE241876 – log2(FC) (R vs NR)",
  title  = "Comparación de efectos (logFC) entre I-SPY2 y GSE241876 a nivel de gen"
)

save_plot(p_genes, "H6_ISPY2_vs_GSE241876_genes_logFC_scatter.png")

message("=== Script 06 (comparación de cohortes) COMPLETADO ===")
