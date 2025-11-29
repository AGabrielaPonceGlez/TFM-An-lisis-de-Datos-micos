# =============================================================================
# 06_comparacion_cohortes.R – Integración multi-cohorte
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# Objetivo general:
#   Comparar los patrones de respuesta entre cohortes a distintos niveles:
#
#   PARTE A – Vías Hallmark (nivel funcional)
#     A1) I-SPY2 (ssGSEA Hallmark + limma) vs GSE241876 (fgsea Hallmark)
#     A2) I-SPY2 (ssGSEA) vs HTG_BrighTNess (ssGSEA)
#     A3) GSE241876 (fgsea) vs HTG_BrighTNess (ssGSEA)
#     -> Tablas de vías comunes, correlaciones de efectos y scatterplots.
#
#   PARTE B – Comparación a nivel de gen (limma)
#     B1) I-SPY2 vs GSE241876 (ya planteado previamente, lo formalizamos aquí)
#     B2) (Opcional) I-SPY2 vs HTG_BrighTNess (si se desea explorar panel común)
#
#   PARTE C – Vista integrativa de IFN-γ entre cohortes
#     C1) Resumen del rendimiento de IFN-γ en I-SPY2 (AUC ROC).
#     C2) Resumen del comportamiento de IFN-γ en HTG_BrighTNess (logFC de vía).
#     C3) Resumen del nivel basal de IFN-γ en HTG_SCANB (distribución de scores).
#     -> Tabla sintética comparando la “firma IFN-γ” entre cohortes.
#
# Todas las tablas se guardan en results/tablas y las figuras en results/figuras.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggplot2)
})

# -------------------------
# Rutas y helpers generales
# -------------------------
path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")

invisible(lapply(c(path_tabs, path_figs),
                 dir.create, recursive = TRUE, showWarnings = FALSE))

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
      x     = x_lab,
      y     = y_lab,
      color = "Categoría"
    ) +
    theme_bw(base_size = 11)
}

# Función auxiliar para añadir flags de significación
add_sig_flags <- function(df, col_p, prefix) {
  df %>%
    dplyr::mutate(
      !!paste0("sig_", prefix) :=
        .data[[col_p]] < 0.05
    )
}

# =============================================================================
# PARTE A – COMPARACIÓN A NIVEL DE VÍAS HALLMARK
# =============================================================================

message("=== PARTE A: Comparación de vías Hallmark entre cohortes ===")

# -----------------------------
# A1) I-SPY2 vs GSE241876
# -----------------------------

# I-SPY2: H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv
file_ispy_path <- file.path(path_tabs, "H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv")
stopifnot(file.exists(file_ispy_path))

ispy2_path <- readr::read_tsv(file_ispy_path, show_col_types = FALSE)

if (!("pathway" %in% names(ispy2_path))) {
  stop("El archivo H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv no tiene columna 'pathway'.")
}

ispy2_clean <- ispy2_path %>%
  dplyr::select(pathway, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_ISPY2 = logFC,
    FDR_ISPY2   = adj.P.Val
  )

# GSE241876: H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv
file_gse_path <- file.path(path_tabs, "H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv")
stopifnot(file.exists(file_gse_path))

gse_path <- readr::read_tsv(file_gse_path, show_col_types = FALSE)

if (!("pathway" %in% names(gse_path))) {
  stop("El archivo H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv no tiene columna 'pathway'.")
}

gse_clean <- gse_path %>%
  dplyr::select(pathway, NES, padj) %>%
  dplyr::rename(
    NES_GSE = NES,
    FDR_GSE = padj
  )

# Unimos por pathway
path_ISPY2_vs_GSE <- dplyr::inner_join(
  ispy2_clean,
  gse_clean,
  by = "pathway"
) %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_GSE   = FDR_GSE   < 0.05,
    sig_ambas = sig_ISPY2 & sig_GSE,
    categoria = dplyr::case_when(
      sig_ambas             ~ "Significativa en ambas",
      sig_ISPY2 & !sig_GSE  ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_GSE  ~ "Solo GSE241876",
      TRUE                  ~ "No significativa"
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

cor_paths_all <- suppressWarnings(
  cor(path_ISPY2_vs_GSE$logFC_ISPY2,
      path_ISPY2_vs_GSE$NES_GSE,
      method = "spearman",
      use    = "complete.obs")
)

summary_paths_ISPY2_GSE <- tibble::tibble(
  comparacion = "ISPY2_vs_GSE241876",
  n_pathways_comunes       = nrow(path_ISPY2_vs_GSE),
  n_sig_ISPY2              = sum(path_ISPY2_vs_GSE$sig_ISPY2, na.rm = TRUE),
  n_sig_GSE                = sum(path_ISPY2_vs_GSE$sig_GSE,   na.rm = TRUE),
  n_sig_ambas              = sum(path_ISPY2_vs_GSE$sig_ambas, na.rm = TRUE),
  cor_spearman_logFC_vs_NES = cor_paths_all
)

write_tsv_safe(path_ISPY2_vs_GSE,
               "H6_pathways_ISPY2_vs_GSE241876_merged.tsv")
write_tsv_safe(summary_paths_ISPY2_GSE,
               "H6_pathways_ISPY2_vs_GSE241876_resumen.tsv")

p_paths_ISPY2_GSE <- plot_scatter_comparison(
  df      = path_ISPY2_vs_GSE,
  x_var   = "logFC_ISPY2",
  y_var   = "NES_GSE",
  col_var = "categoria",
  x_lab   = "I-SPY2 – log2(FC) (R vs NR) en ssGSEA (Hallmark)",
  y_lab   = "GSE241876 – NES (fgsea, Hallmark)",
  title   = "Vías Hallmark: I-SPY2 vs GSE241876"
)

save_plot(p_paths_ISPY2_GSE,
          "H6_ISPY2_vs_GSE241876_Hallmark_scatter.png")


# -----------------------------
# A2) I-SPY2 vs HTG_BrighTNess
# -----------------------------

file_brtn_path <- file.path(path_tabs,
                            "H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways.tsv")
stopifnot(file.exists(file_brtn_path))

brtn_path <- readr::read_tsv(file_brtn_path, show_col_types = FALSE)

if (!("pathway" %in% names(brtn_path))) {
  stop("El archivo H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways.tsv no tiene columna 'pathway'.")
}

brtn_clean <- brtn_path %>%
  dplyr::select(pathway, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_BRTN = logFC,
    FDR_BRTN   = adj.P.Val
  )

path_ISPY2_vs_BRTN <- dplyr::inner_join(
  ispy2_clean,
  brtn_clean,
  by = "pathway"
) %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_BRTN  = FDR_BRTN  < 0.05,
    sig_ambas = sig_ISPY2 & sig_BRTN,
    categoria = dplyr::case_when(
      sig_ambas            ~ "Significativa en ambas",
      sig_ISPY2 & !sig_BRTN ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_BRTN ~ "Solo BrighTNess",
      TRUE                  ~ "No significativa"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativa en ambas",
        "Solo I-SPY2",
        "Solo BrighTNess",
        "No significativa"
      )
    )
  )

cor_paths_ISPY2_BRTN <- suppressWarnings(
  cor(path_ISPY2_vs_BRTN$logFC_ISPY2,
      path_ISPY2_vs_BRTN$logFC_BRTN,
      method = "spearman",
      use    = "complete.obs")
)

summary_paths_ISPY2_BRTN <- tibble::tibble(
  comparacion = "ISPY2_vs_HTG_BrighTNess",
  n_pathways_comunes       = nrow(path_ISPY2_vs_BRTN),
  n_sig_ISPY2              = sum(path_ISPY2_vs_BRTN$sig_ISPY2, na.rm = TRUE),
  n_sig_BRTN               = sum(path_ISPY2_vs_BRTN$sig_BRTN,  na.rm = TRUE),
  n_sig_ambas              = sum(path_ISPY2_vs_BRTN$sig_ambas, na.rm = TRUE),
  cor_spearman_logFC_vs_logFC = cor_paths_ISPY2_BRTN
)

write_tsv_safe(path_ISPY2_vs_BRTN,
               "H6_pathways_ISPY2_vs_HTG_BrighTNess_merged.tsv")
write_tsv_safe(summary_paths_ISPY2_BRTN,
               "H6_pathways_ISPY2_vs_HTG_BrighTNess_resumen.tsv")

p_paths_ISPY2_BRTN <- plot_scatter_comparison(
  df      = path_ISPY2_vs_BRTN,
  x_var   = "logFC_ISPY2",
  y_var   = "logFC_BRTN",
  col_var = "categoria",
  x_lab   = "I-SPY2 – log2(FC) (R vs NR) en vías Hallmark",
  y_lab   = "HTG_BrighTNess – log2(FC) (R vs NR) en vías Hallmark",
  title   = "Vías Hallmark: I-SPY2 vs HTG_BrighTNess"
)

save_plot(p_paths_ISPY2_BRTN,
          "H6_ISPY2_vs_HTG_BrighTNess_Hallmark_scatter.png")


# -----------------------------
# A3) GSE241876 vs HTG_BrighTNess
# -----------------------------

path_GSE_vs_BRTN <- dplyr::inner_join(
  gse_clean,
  brtn_clean,
  by = "pathway"
) %>%
  dplyr::mutate(
    sig_GSE  = FDR_GSE  < 0.05,
    sig_BRTN = FDR_BRTN < 0.05,
    sig_ambas = sig_GSE & sig_BRTN,
    categoria = dplyr::case_when(
      sig_ambas            ~ "Significativa en ambas",
      sig_GSE & !sig_BRTN  ~ "Solo GSE241876",
      !sig_GSE & sig_BRTN  ~ "Solo BrighTNess",
      TRUE                 ~ "No significativa"
    ),
    categoria = factor(
      categoria,
      levels = c(
        "Significativa en ambas",
        "Solo GSE241876",
        "Solo BrighTNess",
        "No significativa"
      )
    )
  )

cor_paths_GSE_BRTN <- suppressWarnings(
  cor(path_GSE_vs_BRTN$NES_GSE,
      path_GSE_vs_BRTN$logFC_BRTN,
      method = "spearman",
      use    = "complete.obs")
)

summary_paths_GSE_BRTN <- tibble::tibble(
  comparacion = "GSE241876_vs_HTG_BrighTNess",
  n_pathways_comunes       = nrow(path_GSE_vs_BRTN),
  n_sig_GSE                = sum(path_GSE_vs_BRTN$sig_GSE,   na.rm = TRUE),
  n_sig_BRTN               = sum(path_GSE_vs_BRTN$sig_BRTN,  na.rm = TRUE),
  n_sig_ambas              = sum(path_GSE_vs_BRTN$sig_ambas, na.rm = TRUE),
  cor_spearman_NES_vs_logFC = cor_paths_GSE_BRTN
)

write_tsv_safe(path_GSE_vs_BRTN,
               "H6_pathways_GSE241876_vs_HTG_BrighTNess_merged.tsv")
write_tsv_safe(summary_paths_GSE_BRTN,
               "H6_pathways_GSE241876_vs_HTG_BrighTNess_resumen.tsv")

p_paths_GSE_BRTN <- plot_scatter_comparison(
  df      = path_GSE_vs_BRTN,
  x_var   = "NES_GSE",
  y_var   = "logFC_BRTN",
  col_var = "categoria",
  x_lab   = "GSE241876 – NES (fgsea, Hallmark)",
  y_lab   = "HTG_BrighTNess – log2(FC) (R vs NR) en vías Hallmark",
  title   = "Vías Hallmark: GSE241876 vs HTG_BrighTNess"
)

save_plot(p_paths_GSE_BRTN,
          "H6_GSE241876_vs_HTG_BrighTNess_Hallmark_scatter.png")


# =============================================================================
# PARTE B – COMPARACIÓN A NIVEL DE GEN (limma)
# =============================================================================

message("=== PARTE B: Comparación a nivel de gen (I-SPY2 vs GSE241876) ===")

# -----------------------------
# B1) I-SPY2 vs GSE241876 (genes)
# -----------------------------

# I-SPY2: limma R vs NR
file_ispy_genes <- file.path(path_tabs, "H4_ISPY2_limma_R_vs_NR.tsv")
if (!file.exists(file_ispy_genes)) {
  # fallback por si el nombre es H3_...
  file_ispy_genes_alt <- file.path(path_tabs, "H3_ISPY2_limma_R_vs_NR.tsv")
  if (file.exists(file_ispy_genes_alt)) {
    file_ispy_genes <- file_ispy_genes_alt
  } else {
    stop("No se encontró ni H4_ISPY2_limma_R_vs_NR.tsv ni H3_ISPY2_limma_R_vs_NR.tsv.")
  }
}

ispy2_genes <- readr::read_tsv(file_ispy_genes, show_col_types = FALSE)

if (!("Gene" %in% names(ispy2_genes))) {
  ispy2_genes <- ispy2_genes %>%
    tibble::rownames_to_column("Gene")
}

ispy2_genes_clean <- ispy2_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_ISPY2 = logFC,
    FDR_ISPY2   = adj.P.Val
  )

# GSE241876: limma-voom R vs NR
file_gse_genes <- file.path(path_tabs, "H5_GSE241876_limma_R_vs_NR.tsv")
stopifnot(file.exists(file_gse_genes))

gse_genes <- readr::read_tsv(file_gse_genes, show_col_types = FALSE)
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

genes_ISPY2_vs_GSE <- dplyr::inner_join(
  ispy2_genes_clean,
  gse_genes_clean,
  by = "Gene"
) %>%
  dplyr::mutate(
    sig_ISPY2 = FDR_ISPY2 < 0.05,
    sig_GSE   = FDR_GSE   < 0.05,
    sig_ambas = sig_ISPY2 & sig_GSE,
    categoria = dplyr::case_when(
      sig_ambas             ~ "Significativo en ambas",
      sig_ISPY2 & !sig_GSE  ~ "Solo I-SPY2",
      !sig_ISPY2 & sig_GSE  ~ "Solo GSE241876",
      TRUE                  ~ "No significativo"
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

cor_genes_all <- suppressWarnings(
  cor(genes_ISPY2_vs_GSE$logFC_ISPY2,
      genes_ISPY2_vs_GSE$logFC_GSE,
      method = "spearman",
      use    = "complete.obs")
)

if (any(genes_ISPY2_vs_GSE$sig_ambas, na.rm = TRUE)) {
  df_sig_both <- genes_ISPY2_vs_GSE %>%
    dplyr::filter(sig_ambas,
                  is.finite(logFC_ISPY2),
                  is.finite(logFC_GSE))
  if (nrow(df_sig_both) > 1) {
    cor_genes_sig_both <- suppressWarnings(
      cor(df_sig_both$logFC_ISPY2,
          df_sig_both$logFC_GSE,
          method = "spearman",
          use    = "complete.obs")
    )
  } else {
    cor_genes_sig_both <- NA_real_
  }
} else {
  cor_genes_sig_both <- NA_real_
}

summary_genes_ISPY2_GSE <- tibble::tibble(
  comparacion               = "ISPY2_vs_GSE241876",
  n_genes_comunes           = nrow(genes_ISPY2_vs_GSE),
  n_sig_ISPY2               = sum(genes_ISPY2_vs_GSE$sig_ISPY2, na.rm = TRUE),
  n_sig_GSE                 = sum(genes_ISPY2_vs_GSE$sig_GSE,   na.rm = TRUE),
  n_sig_ambas               = sum(genes_ISPY2_vs_GSE$sig_ambas, na.rm = TRUE),
  cor_spearman_logFC_global = cor_genes_all,
  cor_spearman_logFC_sig_both = cor_genes_sig_both
)

write_tsv_safe(genes_ISPY2_vs_GSE,
               "H6_genes_ISPY2_vs_GSE241876_merged.tsv")
write_tsv_safe(summary_genes_ISPY2_GSE,
               "H6_genes_ISPY2_vs_GSE241876_resumen.tsv")

genes_plot <- genes_ISPY2_vs_GSE %>%
  dplyr::filter(is.finite(logFC_ISPY2),
                is.finite(logFC_GSE))

p_genes_ISPY2_GSE <- plot_scatter_comparison(
  df      = genes_plot,
  x_var   = "logFC_ISPY2",
  y_var   = "logFC_GSE",
  col_var = "categoria",
  x_lab   = "I-SPY2 – log2(FC) (R vs NR)",
  y_lab   = "GSE241876 – log2(FC) (R vs NR)",
  title   = "Efectos a nivel de gen: I-SPY2 vs GSE241876"
)

save_plot(p_genes_ISPY2_GSE,
          "H6_ISPY2_vs_GSE241876_genes_logFC_scatter.png")


# =============================================================================
# PARTE C – VISTA INTEGRATIVA DE LA FIRMA IFN-γ ENTRE COHORTES
# =============================================================================

message("=== PARTE C: Vista integrativa de IFN-γ (Hallmark) entre cohortes ===")

# I-SPY2: AUC IFN-γ
file_ifng_ispy_auc <- file.path(path_tabs, "H7_ISPY2_IFNG_AUC.tsv")
if (file.exists(file_ifng_ispy_auc)) {
  ifng_ispy_auc <- readr::read_tsv(file_ifng_ispy_auc, show_col_types = FALSE)
} else {
  warning("No se encontró H7_ISPY2_IFNG_AUC.tsv; se dejará NA para AUC IFN-γ en I-SPY2.")
  ifng_ispy_auc <- tibble::tibble(
    marker = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    AUC    = NA_real_,
    cohort = "GSE173839_ISPY2",
    n_R    = NA_real_,
    n_NR   = NA_real_
  )
}

# HTG_BrighTNess: efecto IFN-γ a nivel de vía
path_ifng <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
ifng_brtn_row <- brtn_clean %>%
  dplyr::filter(pathway == path_ifng)

if (nrow(ifng_brtn_row) == 0) {
  warning("En HTG_BrighTNess no se encontró la vía IFN-γ en la tabla de vías.")
  ifng_brtn_summary <- tibble::tibble(
    cohort     = "HTG_BrighTNess",
    pathway    = path_ifng,
    logFC      = NA_real_,
    FDR        = NA_real_
  )
} else {
  ifng_brtn_summary <- tibble::tibble(
    cohort  = "HTG_BrighTNess",
    pathway = path_ifng,
    logFC   = ifng_brtn_row$logFC_BRTN,
    FDR     = ifng_brtn_row$FDR_BRTN
  )
}

# HTG_SCANB: resumen descriptivo de scores IFN-γ
file_ifng_scanb <- file.path(path_tabs, "H5_HTG_SCANB_IFNG_score_resumen.tsv")
if (file.exists(file_ifng_scanb)) {
  ifng_scanb_summary <- readr::read_tsv(file_ifng_scanb, show_col_types = FALSE)
} else {
  warning("No se encontró H5_HTG_SCANB_IFNG_score_resumen.tsv; se dejarán NA para SCAN-B.")
  ifng_scanb_summary <- tibble::tibble(
    cohort     = "HTG_SCANB",
    pathway    = path_ifng,
    n_samples  = NA_real_,
    mean_score = NA_real_,
    sd_score   = NA_real_,
    q25        = NA_real_,
    q50        = NA_real_,
    q75        = NA_real_
  )
}

# Construimos una tabla integrativa compacta para la firma IFN-γ
ifng_integrativo <- tibble::tibble(
  cohorte = c("GSE173839_ISPY2", "HTG_BrighTNess", "HTG_SCANB"),
  tipo_medida = c("AUC ROC (score vía IFN-γ)",
                  "log2(FC) R vs NR (vía IFN-γ)",
                  "Score basal ssGSEA IFN-γ"),
  valor_principal = c(
    ifng_ispy_auc$AUC[1],
    ifng_brtn_summary$logFC[1],
    ifng_scanb_summary$mean_score[1]
  ),
  info_adicional = c(
    paste0("n_R = ", ifng_ispy_auc$n_R[1],
           ", n_NR = ", ifng_ispy_auc$n_NR[1]),
    paste0("FDR = ",
           if (!is.na(ifng_brtn_summary$FDR[1]))
             signif(ifng_brtn_summary$FDR[1], 3) else "NA"),
    paste0("n = ", ifng_scanb_summary$n_samples[1],
           "; q25/q50/q75 = ",
           if (!is.na(ifng_scanb_summary$q25[1]))
             paste0(
               signif(ifng_scanb_summary$q25[1], 3), "/",
               signif(ifng_scanb_summary$q50[1], 3), "/",
               signif(ifng_scanb_summary$q75[1], 3)
             ) else "NA")
  )
)

write_tsv_safe(ifng_integrativo,
               "H6_IFNG_integrativo_cohortes.tsv")
print(ifng_integrativo)

message("=== Script 06 (comparación de cohortes, vías + genes + IFN-γ) COMPLETADO ===")
