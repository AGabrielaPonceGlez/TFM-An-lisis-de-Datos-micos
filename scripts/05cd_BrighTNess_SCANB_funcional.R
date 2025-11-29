# =============================================================================
# 05c_HTG_BrighTNess_SCANB_funcional.R
# Hito 5c y 5d – Análisis funcional en cohortes HTG:
#   - 5c) HTG_BrighTNess: ssGSEA Hallmark + limma a nivel de vía (R vs NR)
#   - 5d) HTG_SCANB: ssGSEA Hallmark (descriptivo, sin R/NR)
#
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
#
# Requisitos previos:
#   - Haber ejecutado 02_curacion_datos.R (meta_master.tsv generado).
#   - Tener disponible data/G9_HTG/valid_datasets_HTG.RData.
#   - Haber ejecutado (opcional pero recomendable) 05_analisis_funcional.R
#     para mantener coherencia en el uso de Hallmark.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(limma)
  library(GSVA)
  library(msigdbr)
})

# -----------------------------------------------------------------------------
# 0) Rutas y funciones auxiliares
# -----------------------------------------------------------------------------

path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")
invisible(lapply(c(path_tabs, path_figs), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

to_upper_trim <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)
}

write_tsv_safe <- function(df, filename) {
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))
}

save_plot <- function(p, filename, w = 8, h = 6, dpi = 300) {
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))
}

# A partir de un objeto HTG (matrix o data.frame) devolvemos una matriz
# genes x muestras con rownames = símbolo de gen (o ID disponible)
get_htg_matrix <- function(obj) {
  if (is.matrix(obj)) {
    # Asumimos que las filas son genes y las rownames ya contienen IDs
    if (is.null(rownames(obj))) {
      stop("La matriz HTG no tiene rownames; no podemos identificar los genes.")
    }
    return(obj)
  }
  if (is.data.frame(obj)) {
    if (ncol(obj) < 2) {
      stop("El data.frame HTG tiene menos de 2 columnas.")
    }
    genes <- obj[[1]]
    mat   <- as.matrix(obj[, -1, drop = FALSE])
    rownames(mat) <- as.character(genes)
    suppressWarnings(mode(mat) <- "numeric")
    return(mat)
  }
  stop("Objeto HTG no soportado (ni matrix ni data.frame).")
}

# Colapsa genes duplicados usando la media de filas (como en script 05)
collapse_by_gene <- function(X, genes) {
  ord   <- order(genes)
  genes <- genes[ord]
  X     <- X[ord, , drop = FALSE]
  X2    <- rowsum(X, group = genes) / as.vector(table(genes))
  list(expr = X2, genes = rownames(X2))
}

# Gráfico tipo barplot para vías (logFC en eje Y, top_n vías)
plot_pathway_bar <- function(tt, title, top_n = 15) {
  tt2 <- tt %>%
    dplyr::mutate(
      pathway = rownames(.),
      sig     = adj.P.Val < 0.05
    ) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  ggplot(tt2, aes(x = pathway, y = logFC, fill = sig)) +
    geom_col() +
    coord_flip() +
    labs(
      title = title,
      x     = "Vía Hallmark",
      y     = "log2(FC) (R vs NR)"
    ) +
    theme_bw(base_size = 11)
}

# -----------------------------------------------------------------------------
# 1) Cargar meta_master y entorno HTG
# -----------------------------------------------------------------------------

meta_master <- readr::read_tsv(
  file.path(path_tabs, "meta_master.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(across(
    c(cohort, technology, sample_id, arm, response, timepoint),
    as.character
  ))

# Cargamos el RData con las matrices HTG (BrighTNess, MEDI, SCAN-B, etc.)
htg_env <- new.env()
path_htg_rdata <- file.path("data", "G9_HTG", "valid_datasets_HTG.RData")
if (file.exists(path_htg_rdata)) {
  load(path_htg_rdata, envir = htg_env)
  message("HTG cargado desde: ", path_htg_rdata)
} else {
  stop("No se encontró valid_datasets_HTG.RData; ejecutar 01/02 o revisar ruta.")
}

# -----------------------------------------------------------------------------
# 2) Colecciones Hallmark (MSigDB) – mismas que en 05_analisis_funcional.R
# -----------------------------------------------------------------------------

msig_h <- msigdbr::msigdbr(
  species    = "Homo sapiens",
  collection = "H"
)

hallmark_list <- msig_h %>%
  dplyr::select(gs_name, gene_symbol) %>%
  dplyr::group_by(gs_name) %>%
  dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  { setNames(.$genes, .$gs_name) }

message("Vías Hallmark cargadas: ", length(hallmark_list))

path_ifng <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
if (!(path_ifng %in% names(hallmark_list))) {
  warning("No se encontró la vía HALLMARK_INTERFERON_GAMMA_RESPONSE en Hallmark.")
}

# =============================================================================
# 5c) HTG_BrighTNess – ssGSEA Hallmark + limma a nivel de vía (R vs NR)
# =============================================================================

message("=== 5c) HTG_BrighTNess – ssGSEA Hallmark + limma vías (R vs NR) ===")

if (!exists("brtn.htg.Rseq", envir = htg_env)) {
  warning("No se encontró 'brtn.htg.Rseq' en htg_env; se omite 5c (BrighTNess).")
} else {
  # -----------------------------
  # 5c.1) Matriz de expresión HTG BrighTNess
  # -----------------------------
  brtn_htg   <- get("brtn.htg.Rseq", envir = htg_env)
  X_brtn_raw <- get_htg_matrix(brtn_htg)   # genes x muestras
  
  # Homogeneizamos IDs de muestra: mayúsculas, sin 'X' inicial (como en script 04)
  colnames(X_brtn_raw) <- to_upper_trim(colnames(X_brtn_raw))
  colnames(X_brtn_raw) <- sub("^X", "", colnames(X_brtn_raw))
  
  genes_brtn <- rownames(X_brtn_raw)
  
  # Colapsamos posibles genes duplicados
  collapsed_brtn <- collapse_by_gene(X_brtn_raw, genes_brtn)
  X_brtn         <- collapsed_brtn$expr
  genes_brtn     <- collapsed_brtn$genes
  
  message("HTG_BrighTNess – genes tras colapso: ", nrow(X_brtn))
  message("HTG_BrighTNess – muestras totales: ", ncol(X_brtn))
  
  # -----------------------------
  # 5c.2) Fenodata: R vs NR en meta_master
  # -----------------------------
  pheno_brtn <- meta_master %>%
    dplyr::filter(
      cohort   == "HTG_BrighTNess",
      response %in% c("R", "NR")
    ) %>%
    dplyr::mutate(
      sample_id = to_upper_trim(sample_id),
      response  = factor(response, levels = c("NR", "R"))
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  # Alineamos por intersección de IDs
  common_ids <- intersect(colnames(X_brtn), pheno_brtn$sample_id)
  message("HTG_BrighTNess – muestras con fenotipo R/NR: ", length(common_ids))
  
  if (length(common_ids) < 4 ||
      length(unique(pheno_brtn$response[pheno_brtn$sample_id %in% common_ids])) < 2) {
    warning("HTG_BrighTNess: muy pocas muestras R/NR tras emparejar; ",
            "no se realizará limma a nivel de vía.")
  } else {
    # Restricción a muestras emparejadas
    X_brtn  <- X_brtn[, common_ids, drop = FALSE]
    pheno_brtn <- pheno_brtn[match(common_ids, pheno_brtn$sample_id),
                             , drop = FALSE]
    
    stopifnot(identical(colnames(X_brtn), pheno_brtn$sample_id))
    
    # -----------------------------
    # 5c.3) ssGSEA Hallmark (GSVA API nueva)
    # -----------------------------
    message("Calculando ssGSEA (Hallmark) para HTG_BrighTNess...")
    
    # Aseguramos que rownames(X_brtn) son símbolos de gen (HTG OBP lo suele ser)
    # y que no hay NAs
    keep_genes <- !is.na(rownames(X_brtn)) & rownames(X_brtn) != ""
    X_brtn_ssg <- X_brtn[keep_genes, , drop = FALSE]
    
    ssgsea_param_brtn <- GSVA::ssgseaParam(
      exprData = X_brtn_ssg,
      geneSets = hallmark_list
    )
    
    scores_brtn <- GSVA::gsva(
      ssgsea_param_brtn,
      verbose = FALSE
    )
    
    message("HTG_BrighTNess – ssGSEA dim (vías x muestras): ",
            paste(dim(scores_brtn), collapse = " x "))
    
    # Guardamos tabla de scores vías x muestras
    scores_brtn_df <- as.data.frame(scores_brtn) %>%
      tibble::rownames_to_column("pathway")
    
    write_tsv_safe(scores_brtn_df,
                   "H5_HTG_BrighTNess_ssGSEA_scores.tsv")
    
    # -----------------------------
    # 5c.4) limma a nivel de vía (~ response)
    # -----------------------------
    design_brtn <- stats::model.matrix(~ response, data = pheno_brtn)
    
    fit_brtn <- limma::lmFit(scores_brtn, design_brtn)
    fit_brtn <- limma::eBayes(fit_brtn, trend = TRUE, robust = TRUE)
    
    coef_name_brtn <- "responseR"
    if (!(coef_name_brtn %in% colnames(coef(fit_brtn)))) {
      stop("HTG_BrighTNess: no se encontró el coeficiente 'responseR' en el modelo de vías.")
    }
    
    tt_path_brtn <- limma::topTable(
      fit_brtn,
      coef   = coef_name_brtn,
      number = Inf,
      sort.by = "P"
    )
    
    tt_path_brtn_out <- tt_path_brtn %>%
      tibble::rownames_to_column("pathway")
    
    write_tsv_safe(tt_path_brtn_out,
                   "H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways.tsv")
    write_tsv_safe(
      tt_path_brtn_out %>% dplyr::slice_head(n = 30),
      "H5_HTG_BrighTNess_ssGSEA_R_vs_NR_pathways_top30.tsv"
    )
    
    # Resumen de número de vías significativas
    summary_brtn_path <- tibble::tibble(
      cohort                 = "HTG_BrighTNess",
      n_pathways             = nrow(tt_path_brtn),
      n_sig_fdr_05           = sum(tt_path_brtn$adj.P.Val < 0.05, na.rm = TRUE),
      n_sig_absLFC_0.2_FDR05 = sum(
        tt_path_brtn$adj.P.Val < 0.05 &
          abs(tt_path_brtn$logFC) >= 0.2,
        na.rm = TRUE
      )
    )
    
    write_tsv_safe(summary_brtn_path,
                   "H5_HTG_BrighTNess_ssGSEA_resumen.tsv")
    print(summary_brtn_path)
    
    # -----------------------------
    # 5c.5) Barplot top15 vías por logFC (R vs NR)
    # -----------------------------
    p_bar_brtn <- plot_pathway_bar(
      tt_path_brtn,
      "HTG_BrighTNess – ssGSEA Hallmark: vías asociadas a R vs NR",
      top_n = 15
    )
    
    save_plot(p_bar_brtn,
              "H5_HTG_BrighTNess_ssGSEA_pathways_barplot_top15.png")
    
    message("5c) HTG_BrighTNess – ssGSEA + limma vías COMPLETADO.")
  }
}

# =============================================================================
# 5d) HTG_SCANB – ssGSEA Hallmark descriptivo (sin R/NR)
# =============================================================================

message("=== 5d) HTG_SCANB – ssGSEA Hallmark (descriptivo) ===")

if (!exists("scanb.htg.Rseq", envir = htg_env)) {
  warning("No se encontró 'scanb.htg.Rseq' en htg_env; se omite 5d (SCAN-B).")
} else {
  # -----------------------------
  # 5d.1) Matriz de expresión HTG SCAN-B
  # -----------------------------
  scanb_htg   <- get("scanb.htg.Rseq", envir = htg_env)
  X_scanb_raw <- get_htg_matrix(scanb_htg)  # genes x muestras
  
  colnames(X_scanb_raw) <- to_upper_trim(colnames(X_scanb_raw))
  colnames(X_scanb_raw) <- sub("^X", "", colnames(X_scanb_raw))
  
  genes_scanb <- rownames(X_scanb_raw)
  
  collapsed_scanb <- collapse_by_gene(X_scanb_raw, genes_scanb)
  X_scanb         <- collapsed_scanb$expr
  genes_scanb     <- collapsed_scanb$genes
  
  message("HTG_SCANB – genes tras colapso: ", nrow(X_scanb))
  message("HTG_SCANB – muestras totales: ", ncol(X_scanb))
  
  # -----------------------------
  # 5d.2) ssGSEA Hallmark (sin fenotipo, solo descriptivo)
  # -----------------------------
  keep_genes <- !is.na(rownames(X_scanb)) & rownames(X_scanb) != ""
  X_scanb_ssg <- X_scanb[keep_genes, , drop = FALSE]
  
  message("Calculando ssGSEA (Hallmark) para HTG_SCANB...")
  
  ssgsea_param_scanb <- GSVA::ssgseaParam(
    exprData = X_scanb_ssg,
    geneSets = hallmark_list
  )
  
  scores_scanb <- GSVA::gsva(
    ssgsea_param_scanb,
    verbose = FALSE
  )
  
  message("HTG_SCANB – ssGSEA dim (vías x muestras): ",
          paste(dim(scores_scanb), collapse = " x "))
  
  scores_scanb_df <- as.data.frame(scores_scanb) %>%
    tibble::rownames_to_column("pathway")
  
  write_tsv_safe(scores_scanb_df,
                 "H5_HTG_SCANB_ssGSEA_scores.tsv")
  
  # -----------------------------
  # 5d.3) Resumen descriptivo (ej. distribución IFN-γ)
  # -----------------------------
  if (path_ifng %in% rownames(scores_scanb)) {
    ifng_scanb <- as.numeric(scores_scanb[path_ifng, ])
    
    summary_ifng_scanb <- tibble::tibble(
      cohort   = "HTG_SCANB",
      pathway  = path_ifng,
      n_samples = length(ifng_scanb),
      mean_score = mean(ifng_scanb, na.rm = TRUE),
      sd_score   = sd(ifng_scanb, na.rm = TRUE),
      q25        = quantile(ifng_scanb, 0.25, na.rm = TRUE),
      q50        = quantile(ifng_scanb, 0.50, na.rm = TRUE),
      q75        = quantile(ifng_scanb, 0.75, na.rm = TRUE)
    )
    
    write_tsv_safe(summary_ifng_scanb,
                   "H5_HTG_SCANB_IFNG_score_resumen.tsv")
    print(summary_ifng_scanb)
    
    # Histograma simple de los scores IFN-γ en SCAN-B
    df_ifng_scanb <- tibble::tibble(
      sample_id = colnames(scores_scanb),
      IFNG_score = ifng_scanb
    )
    
    p_ifng_scanb <- ggplot(df_ifng_scanb, aes(x = IFNG_score)) +
      geom_histogram(bins = 30, fill = "grey70", color = "black") +
      labs(
        title = "HTG_SCANB – Distribución de scores IFN-γ (Hallmark, ssGSEA)",
        x     = "Score ssGSEA IFN-γ",
        y     = "Frecuencia"
      ) +
      theme_bw(base_size = 11)
    
    save_plot(p_ifng_scanb,
              "H5_HTG_SCANB_IFNG_ssGSEA_hist.png")
    
  } else {
    warning("En HTG_SCANB no se encontró la vía HALLMARK_INTERFERON_GAMMA_RESPONSE ",
            "en los scores ssGSEA; no se genera resumen IFN-γ.")
  }
  
  message("5d) HTG_SCANB – ssGSEA Hallmark descriptivo COMPLETADO.")
}

message("=== Script 05c/05d (HTG BrighTNess + HTG SCAN-B funcional) COMPLETADO ===")
