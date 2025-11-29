#!/usr/bin/env Rscript

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

message("=== Script 06.1 – Comparaciones a nivel de gen con cohortes HTG ===")

# ------------------------------------------------------
# 1) Carga de resultados limma a nivel de gen por cohorte
# ------------------------------------------------------

# I-SPY2: limma R vs NR (H4 o H3 como alternativa)
file_ispy_genes <- file.path(path_tabs, "H4_ISPY2_limma_R_vs_NR.tsv")
if (!file.exists(file_ispy_genes)) {
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

# HTG_MEDI4736: limma R vs NR
file_medi_genes <- file.path(path_tabs, "H4_HTG_MEDI4736_limma_R_vs_NR.tsv")
stopifnot(file.exists(file_medi_genes))
medi_genes <- readr::read_tsv(file_medi_genes, show_col_types = FALSE)

if (!("Gene" %in% names(medi_genes))) {
  medi_genes <- medi_genes %>%
    tibble::rownames_to_column("Gene")
}

medi_genes_clean <- medi_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_MEDI = logFC,
    FDR_MEDI   = adj.P.Val
  )

# HTG_BrighTNess: limma R vs NR
file_brtn_genes <- file.path(path_tabs, "H4_HTG_BrighTNess_limma_R_vs_NR.tsv")
stopifnot(file.exists(file_brtn_genes))
brtn_genes <- readr::read_tsv(file_brtn_genes, show_col_types = FALSE)

if (!("Gene" %in% names(brtn_genes))) {
  brtn_genes <- brtn_genes %>%
    tibble::rownames_to_column("Gene")
}

brtn_genes_clean <- brtn_genes %>%
  dplyr::select(Gene, logFC, adj.P.Val) %>%
  dplyr::rename(
    logFC_BRTN = logFC,
    FDR_BRTN   = adj.P.Val
  )

# ------------------------------------------------------
# 2) Función auxiliar para comparación par a par de genes
# ------------------------------------------------------

comparar_genes <- function(df1, df2,
                           cohort1_label, cohort2_label,
                           logFC1, FDR1, logFC2, FDR2,
                           prefix_out) {
  
  # Unimos por gen y creamos categorías de significación
  merged <- dplyr::inner_join(df1, df2, by = "Gene") %>%
    dplyr::mutate(
      sig_1   = .data[[FDR1]] < 0.05,
      sig_2   = .data[[FDR2]] < 0.05,
      sig_ambas = sig_1 & sig_2,
      categoria = dplyr::case_when(
        sig_ambas            ~ "Significativo en ambas",
        sig_1 & !sig_2       ~ paste0("Solo ", cohort1_label),
        !sig_1 & sig_2       ~ paste0("Solo ", cohort2_label),
        TRUE                 ~ "No significativo"
      ),
      categoria = factor(
        categoria,
        levels = c(
          "Significativo en ambas",
          paste0("Solo ", cohort1_label),
          paste0("Solo ", cohort2_label),
          "No significativo"
        )
      )
    )
  
  # Correlación global de efectos
  cor_global <- suppressWarnings(
    cor(merged[[logFC1]],
        merged[[logFC2]],
        method = "spearman",
        use    = "complete.obs")
  )
  
  # Correlación solo en genes significativos en ambas
  if (any(merged$sig_ambas, na.rm = TRUE)) {
    df_sig_both <- merged %>%
      dplyr::filter(sig_ambas,
                    is.finite(.data[[logFC1]]),
                    is.finite(.data[[logFC2]]))
    if (nrow(df_sig_both) > 1) {
      cor_sig_both <- suppressWarnings(
        cor(df_sig_both[[logFC1]],
            df_sig_both[[logFC2]],
            method = "spearman",
            use    = "complete.obs")
      )
    } else {
      cor_sig_both <- NA_real_
    }
  } else {
    cor_sig_both <- NA_real_
  }
  
  resumen <- tibble::tibble(
    comparacion                 = paste0(cohort1_label, "_vs_", cohort2_label),
    n_genes_comunes             = nrow(merged),
    n_sig_1                     = sum(merged$sig_1,   na.rm = TRUE),
    n_sig_2                     = sum(merged$sig_2,   na.rm = TRUE),
    n_sig_ambas                 = sum(merged$sig_ambas, na.rm = TRUE),
    cor_spearman_logFC_global   = cor_global,
    cor_spearman_logFC_sig_both = cor_sig_both
  )
  
  # Guardamos tablas
  write_tsv_safe(
    merged,
    paste0("H6_genes_", prefix_out, "_merged.tsv")
  )
  write_tsv_safe(
    resumen,
    paste0("H6_genes_", prefix_out, "_resumen.tsv")
  )
  
  # Scatter
  merged_plot <- merged %>%
    dplyr::filter(
      is.finite(.data[[logFC1]]),
      is.finite(.data[[logFC2]])
    )
  
  p <- plot_scatter_comparison(
    df      = merged_plot,
    x_var   = logFC1,
    y_var   = logFC2,
    col_var = "categoria",
    x_lab   = paste0(cohort1_label, " – log2(FC) (R vs NR)"),
    y_lab   = paste0(cohort2_label, " – log2(FC) (R vs NR)"),
    title   = paste0("Efectos a nivel de gen: ",
                     cohort1_label, " vs ", cohort2_label)
  )
  
  save_plot(
    p,
    paste0("H6_", prefix_out, "_genes_logFC_scatter.png")
  )
  
  invisible(list(merged = merged, resumen = resumen, plot = p))
}

# ------------------------------------------------------
# 3) Comparaciones específicas que nos interesan
#    - I-SPY2 vs HTG_MEDI4736
#    - I-SPY2 vs HTG_BrighTNess
#    - GSE241876 vs HTG_MEDI4736
#    - GSE241876 vs HTG_BrighTNess
#    - HTG_MEDI4736 vs HTG_BrighTNess
# ------------------------------------------------------

message(">>> Comparación a nivel de gen: I-SPY2 vs HTG_MEDI4736")
res_ISPY2_MEDI <- comparar_genes(
  df1           = ispy2_genes_clean,
  df2           = medi_genes_clean,
  cohort1_label = "I-SPY2",
  cohort2_label = "HTG_MEDI4736",
  logFC1        = "logFC_ISPY2",
  FDR1          = "FDR_ISPY2",
  logFC2        = "logFC_MEDI",
  FDR2          = "FDR_MEDI",
  prefix_out    = "ISPY2_vs_HTG_MEDI4736"
)

message(">>> Comparación a nivel de gen: I-SPY2 vs HTG_BrighTNess")
res_ISPY2_BRTN <- comparar_genes(
  df1           = ispy2_genes_clean,
  df2           = brtn_genes_clean,
  cohort1_label = "I-SPY2",
  cohort2_label = "HTG_BrighTNess",
  logFC1        = "logFC_ISPY2",
  FDR1          = "FDR_ISPY2",
  logFC2        = "logFC_BRTN",
  FDR2          = "FDR_BRTN",
  prefix_out    = "ISPY2_vs_HTG_BrighTNess"
)

message(">>> Comparación a nivel de gen: GSE241876 vs HTG_MEDI4736")
res_GSE_MEDI <- comparar_genes(
  df1           = gse_genes_clean,
  df2           = medi_genes_clean,
  cohort1_label = "GSE241876",
  cohort2_label = "HTG_MEDI4736",
  logFC1        = "logFC_GSE",
  FDR1          = "FDR_GSE",
  logFC2        = "logFC_MEDI",
  FDR2          = "FDR_MEDI",
  prefix_out    = "GSE241876_vs_HTG_MEDI4736"
)

message(">>> Comparación a nivel de gen: GSE241876 vs HTG_BrighTNess")
res_GSE_BRTN <- comparar_genes(
  df1           = gse_genes_clean,
  df2           = brtn_genes_clean,
  cohort1_label = "GSE241876",
  cohort2_label = "HTG_BrighTNess",
  logFC1        = "logFC_GSE",
  FDR1          = "FDR_GSE",
  logFC2        = "logFC_BRTN",
  FDR2          = "FDR_BRTN",
  prefix_out    = "GSE241876_vs_HTG_BrighTNess"
)

message(">>> Comparación a nivel de gen: HTG_MEDI4736 vs HTG_BrighTNess")
res_MEDI_BRTN <- comparar_genes(
  df1           = medi_genes_clean,
  df2           = brtn_genes_clean,
  cohort1_label = "HTG_MEDI4736",
  cohort2_label = "HTG_BrighTNess",
  logFC1        = "logFC_MEDI",
  FDR1          = "FDR_MEDI",
  logFC2        = "logFC_BRTN",
  FDR2          = "FDR_BRTN",
  prefix_out    = "HTG_MEDI4736_vs_HTG_BrighTNess"
)

message("=== Script 06.1 (comparación genes con HTG_MEDI4736 / HTG_BrighTNess) COMPLETADO ===")
