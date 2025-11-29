# =============================================================================
# 04_analisis_diferencial.R – R vs NR por cohorte
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# En este script desarrollamos el análisis diferencial para:
#   4.a) Cohorte I-SPY2 (GSE173839, microarrays Agilent)
#   4.b) Cohorte HTG MEDI4736 – cohorte principal (HTG EdgeSeq OBP, datos NORMALIZADOS)
#   4.c) Cohorte HTG BrighTNess – validación adicional (HTG EdgeSeq OBP, datos NORMALIZADOS)
#
# En las tres cohortes contrastamos la expresión génica entre respondedores (R)
# y no respondedores (NR):
#   - I-SPY2: limma clásico sobre intensidades Agilent, ajustando por brazo
#             de tratamiento (arm: control vs durvalumab/olaparib).
#   - HTG MEDI4736: limma directo sobre matriz HTG normalizada (~ response).
#   - HTG BrighTNess: limma directo sobre matriz HTG normalizada (~ response).
#
# Nota importante:
#   - La cohorte GSE241876 (RNA-seq) NO se analiza con limma aquí porque solo
#     dispone de 2 NR. La usaremos más adelante de forma exploratoria con
#     ssGSEA/GSVA en otro script.
# =============================================================================


# Cargamos las librerías:
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(limma)
  library(edgeR)  # ya no lo usamos en HTG, pero lo dejamos cargado por si hiciera falta en otros scripts
})

# Definimos las rutas de trabajo:
path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
path_figs    <- file.path(path_results, "figuras")
invisible(lapply(c(path_tabs, path_figs), dir.create,
                 recursive = TRUE, showWarnings = FALSE))


#### Funciones de apoyo ####

# Función para limpiar y homogeneizar identificadores:
to_upper_trim <- function(x){
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)
}

# Función para guardar tablas en formato TSV:
write_tsv_safe <- function(df, filename){
  readr::write_tsv(df, file.path(path_tabs, filename))
  message("Tabla guardada: ", file.path(path_tabs, filename))
}

# Función para guardar gráficos en PNG:
save_plot <- function(p, filename, w = 8, h = 6, dpi = 300){
  ggsave(file.path(path_figs, filename), p, width = w, height = h, dpi = dpi)
  message("Figura guardada: ", file.path(path_figs, filename))
}

# Preparamos lista con matriz de expresión y vectores de genes/muestras
# (primera columna = genes, resto = muestras):
prepare_expr <- function(tbl, gene_col = 1){
  stopifnot(is.data.frame(tbl))
  genes <- tbl[[gene_col]]
  mat   <- as.matrix(tbl[, -gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes)
}

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
    # Si no hay coincidencias directas, probamos a eliminar sufijos típicos:
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

# Función para generar un gráfico tipo “volcano plot”:
volcano_plot <- function(tt, title, lfc_thr = log2(1.5), p_thr = 0.05){
  tt2 <- tt %>%
    dplyr::mutate(
      sig = ifelse(
        adj.P.Val < p_thr & abs(logFC) >= lfc_thr,
        "Significativo", "No significativo"
      )
    )
  
  ggplot(tt2, aes(x = logFC, y = -log10(adj.P.Val))) +
    geom_point(aes(shape = sig), alpha = 0.85) +
    geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = 2) +
    geom_hline(yintercept = -log10(p_thr), linetype = 2) +
    labs(title = title, x = "log2(FC)", y = "-log10(FDR)") +
    theme_bw(base_size = 12)
}

# Función para generar un gráfico MA:
ma_plot <- function(tt, title){
  ggplot(tt, aes(x = AveExpr, y = logFC)) +
    geom_point(alpha = .6) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = title,
         x = "Media log2 expresión",
         y = "log2(FC) R vs NR") +
    theme_bw(base_size = 12)
}

# Función para sacar matriz de expresión HTG (genes x muestras) de htg_env:
get_htg_matrix <- function(obj){
  if (is.matrix(obj)) {
    # Asumimos que las filas son genes (rownames ya contienen IDs)
    return(obj)
  }
  if (is.data.frame(obj)) {
    # Asumimos primera columna = identificador de gen, resto = columnas de expresión
    genes <- obj[[1]]
    mat   <- as.matrix(obj[, -1, drop = FALSE])
    rownames(mat) <- as.character(genes)
    return(mat)
  }
  stop("Objeto HTG no soportado (ni matrix ni data.frame).")
}


#### Cargamos metadatos armonizados ####

meta_master <- readr::read_tsv(
  file.path(path_tabs, "meta_master.tsv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(across(
    c(cohort, technology, sample_id, arm, response, timepoint),
    as.character
  ))

# Cargamos HTG (para MEDI y BrighTNess) desde el RData:
htg_env <- new.env()
path_htg_rdata <- file.path("data", "G9_HTG", "valid_datasets_HTG.RData")
if (file.exists(path_htg_rdata)) {
  load(path_htg_rdata, envir = htg_env)
  message("HTG cargado desde: ", path_htg_rdata)
} else {
  warning("No se encontró valid_datasets_HTG.RData; no se podrán analizar MEDI/BrighTNess.")
}


# ---------------------------------------------------------------------------
# 4.a) Cohorte I-SPY2 (GSE173839) – microarrays Agilent
# ---------------------------------------------------------------------------

# Cargamos la matriz de expresión de I-SPY2 (microarrays Agilent):
ispy_tbl <- data.table::fread(
  file.path(
    "data", "GSE173839",
    "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"
  ),
  sep = "\t", header = TRUE
) %>% tibble::as_tibble()

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

# 3. Construimos fenodata y ordenamos para que coincida con X.
#    Usamos únicamente las muestras clasificadas como R/NR.
pheno <- meta_master %>%
  dplyr::filter(
    cohort    == "GSE173839_ISPY2",
    sample_id %in% samples,
    response  %in% c("R","NR")   # solo respondedores/no respondedores
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR","R")),   # coef 'responseR' = R > NR
    arm      = to_upper_trim(arm),
    arm      = forcats::fct_lump_min(arm, min = 5),      # agrupamos brazos muy raros
    arm      = stats::relevel(factor(arm), ref = "CONTROL")
  )

# Reordenamos filas para que sample_id == columnas de X:
pheno <- pheno[match(samples, pheno$sample_id), , drop = FALSE]
stopifnot(identical(pheno$sample_id, samples))

# 4. Ajustamos limma. Incluimos la variable "response" (efecto principal)
#    y "arm" (covariable de tratamiento).
design <- stats::model.matrix(~ response + arm, data = pheno)
fit    <- limma::lmFit(X, design)
fit    <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

# Extraemos resultados para el contraste R vs NR:
coef_name <- "responseR"
if (!(coef_name %in% colnames(coef(fit)))) {
  stop("No se encontró el coeficiente 'responseR' en el modelo para I-SPY2")
}

tt <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "P") %>%
  tibble::rownames_to_column("rowid") %>%
  dplyr::mutate(Gene = genes[as.integer(rowid)]) %>%
  dplyr::select(Gene, logFC, AveExpr, t, P.Value, adj.P.Val, B)

# 5. Guardamos tablas y figuras para I-SPY2:
write_tsv_safe(tt, "H4_ISPY2_limma_R_vs_NR.tsv")
write_tsv_safe(tt %>% dplyr::slice_head(n = 50),
               "H4_ISPY2_limma_R_vs_NR_top50.tsv")

save_plot(
  volcano_plot(tt, "I-SPY2 – R vs NR (limma, ~ response + arm)"),
  "H4_ISPY2_volcano.png"
)
save_plot(
  ma_plot(tt, "I-SPY2 – MA plot (limma)"),
  "H4_ISPY2_MA.png"
)

# 6. Resumen del número de genes diferenciales significativos en I-SPY2:
summary_ispy <- tibble::tibble(
  cohort                 = "GSE173839_ISPY2",
  n_genes                = nrow(tt),
  n_sig_fdr_05           = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
  n_sig_absLFC_1.5_FDR05 = sum(tt$adj.P.Val < 0.05 &
                                 abs(tt$logFC) >= log2(1.5), na.rm = TRUE)
)
write_tsv_safe(summary_ispy, "H4_ISPY2_resumen.tsv")
print(summary_ispy)

message("Análisis diferencial para I-SPY2 completado.")


# ---------------------------------------------------------------------------
# 4.b) Cohorte HTG MEDI4736 – limma directo (datos HTG normalizados)
# ---------------------------------------------------------------------------

if (exists("medi.htg.Rseq", envir = htg_env)) {
  
  medi_htg   <- get("medi.htg.Rseq", envir = htg_env)
  X_medi_htg <- get_htg_matrix(medi_htg)   # genes x muestras
  
  # Homogeneizamos nombres de columnas (sample IDs de expresión)
  colnames(X_medi_htg) <- to_upper_trim(colnames(X_medi_htg))
  # Por si hubiera prefijos tipo "X12345" (como en BrighTNess), los eliminamos:
  colnames(X_medi_htg) <- sub("^X", "", colnames(X_medi_htg))
  
  # Fenodata a partir de meta_master: solo R/NR en HTG_MEDI4736
  pheno_medi <- meta_master %>%
    dplyr::filter(
      cohort   == "HTG_MEDI4736",
      response %in% c("R","NR")
    ) %>%
    dplyr::mutate(
      sample_id = to_upper_trim(sample_id),
      response  = factor(response, levels = c("NR","R"))  # NR referencia
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  # Alineamos por IDs comunes
  common_ids_medi <- intersect(colnames(X_medi_htg), pheno_medi$sample_id)
  
  if (length(common_ids_medi) == 0) {
    warning("HTG MEDI4736: no hay IDs comunes entre expresión y meta_master.")
  } else {
    
    X_medi_htg  <- X_medi_htg[, common_ids_medi, drop = FALSE]
    pheno_medi  <- pheno_medi[match(common_ids_medi, pheno_medi$sample_id), , drop = FALSE]
    
    message("HTG MEDI4736 – muestras R/NR para limma: ", ncol(X_medi_htg))
    
    if (ncol(X_medi_htg) >= 4 &&
        length(unique(pheno_medi$response)) == 2) {
      
      # Diseño ~ response (datos normalizados, limma directo)
      design_medi <- stats::model.matrix(~ response, data = pheno_medi)
      
      fit_medi <- limma::lmFit(X_medi_htg, design_medi)
      fit_medi <- limma::eBayes(fit_medi, trend = TRUE, robust = TRUE)
      
      coef_name_medi <- "responseR"
      if (!(coef_name_medi %in% colnames(coef(fit_medi)))) {
        stop("No se encontró el coeficiente 'responseR' en el modelo para HTG MEDI4736")
      }
      
      tt_medi <- limma::topTable(
        fit_medi, coef = coef_name_medi,
        number = Inf, sort.by = "P"
      ) %>%
        tibble::rownames_to_column("Gene")
      
      # Guardamos tablas y figuras para HTG MEDI4736:
      write_tsv_safe(tt_medi, "H4_HTG_MEDI4736_limma_R_vs_NR.tsv")
      write_tsv_safe(
        tt_medi %>% dplyr::slice_head(n = 50),
        "H4_HTG_MEDI4736_limma_R_vs_NR_top50.tsv"
      )
      
      save_plot(
        volcano_plot(tt_medi, "HTG MEDI4736 – R vs NR (limma)"),
        "H4_HTG_MEDI4736_volcano.png"
      )
      save_plot(
        ma_plot(tt_medi, "HTG MEDI4736 – MA plot (limma)"),
        "H4_HTG_MEDI4736_MA.png"
      )
      
      summary_medi <- tibble::tibble(
        cohort                 = "HTG_MEDI4736",
        n_genes                = nrow(tt_medi),
        n_sig_fdr_05           = sum(tt_medi$adj.P.Val < 0.05, na.rm = TRUE),
        n_sig_absLFC_1.5_FDR05 = sum(
          tt_medi$adj.P.Val < 0.05 & abs(tt_medi$logFC) >= log2(1.5),
          na.rm = TRUE
        )
      )
      write_tsv_safe(summary_medi, "H4_HTG_MEDI4736_resumen.tsv")
      print(summary_medi)
      
    } else {
      warning("HTG MEDI4736 no tiene suficientes muestras R/NR para limma (tras emparejar IDs).")
    }
  }
  
} else {
  message("Objeto 'medi.htg.Rseq' no encontrado en htg_env; se omite análisis HTG MEDI4736.")
}


# ---------------------------------------------------------------------------
# 4.c) Cohorte HTG BrighTNess – limma directo (datos HTG normalizados)
# ---------------------------------------------------------------------------

if (exists("brtn.htg.Rseq", envir = htg_env)) {
  
  brtn_htg   <- get("brtn.htg.Rseq", envir = htg_env)
  X_brtn_htg <- get_htg_matrix(brtn_htg)   # genes x muestras
  
  # Homogeneizamos nombres de columnas (sample IDs de expresión)
  colnames(X_brtn_htg) <- to_upper_trim(colnames(X_brtn_htg))
  # En BrighTNess ya vimos que venían como "X102001", etc:
  colnames(X_brtn_htg) <- sub("^X", "", colnames(X_brtn_htg))
  
  # Fenodata a partir de meta_master: solo R/NR en HTG_BrighTNess
  pheno_brtn <- meta_master %>%
    dplyr::filter(
      cohort   == "HTG_BrighTNess",
      response %in% c("R","NR")
    ) %>%
    dplyr::mutate(
      sample_id = to_upper_trim(sample_id),
      response  = factor(response, levels = c("NR","R"))
    ) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE)
  
  # Alineamos por intersección de IDs entre expresión y fenodata
  common_ids <- intersect(colnames(X_brtn_htg), pheno_brtn$sample_id)
  
  if (length(common_ids) == 0) {
    warning("HTG BrighTNess: no hay IDs comunes entre la matriz de expresión y meta_master.")
  } else {
    # Reordenamos ambos objetos usando los IDs comunes
    X_brtn_htg <- X_brtn_htg[, common_ids, drop = FALSE]
    pheno_brtn <- pheno_brtn[match(common_ids, pheno_brtn$sample_id), , drop = FALSE]
    
    message("HTG BrighTNess – muestras R/NR para limma (IDs comunes): ", ncol(X_brtn_htg))
    
    if (ncol(X_brtn_htg) >= 4 &&
        length(unique(pheno_brtn$response)) == 2) {
      
      # Diseño ~ response (datos HTG normalizados, limma directo)
      design_brtn <- stats::model.matrix(~ response, data = pheno_brtn)
      
      fit_b <- limma::lmFit(X_brtn_htg, design_brtn)
      fit_b <- limma::eBayes(fit_b, trend = TRUE, robust = TRUE)
      
      coef_name_b <- "responseR"
      if (!(coef_name_b %in% colnames(coef(fit_b)))) {
        stop("No se encontró el coeficiente 'responseR' en el modelo para HTG BrighTNess")
      }
      
      tt_brtn <- limma::topTable(
        fit_b, coef = coef_name_b,
        number = Inf, sort.by = "P"
      ) %>%
        tibble::rownames_to_column("Gene")
      
      write_tsv_safe(tt_brtn, "H4_HTG_BrighTNess_limma_R_vs_NR.tsv")
      write_tsv_safe(
        tt_brtn %>% dplyr::slice_head(n = 50),
        "H4_HTG_BrighTNess_limma_R_vs_NR_top50.tsv"
      )
      
      save_plot(
        volcano_plot(tt_brtn, "HTG BrighTNess – R vs NR (limma)"),
        "H4_HTG_BrighTNess_volcano.png"
      )
      save_plot(
        ma_plot(tt_brtn, "HTG BrighTNess – MA plot (limma)"),
        "H4_HTG_BrighTNess_MA.png"
      )
      
      summary_brtn <- tibble::tibble(
        cohort                 = "HTG_BrighTNess",
        n_genes                = nrow(tt_brtn),
        n_sig_fdr_05           = sum(tt_brtn$adj.P.Val < 0.05, na.rm = TRUE),
        n_sig_absLFC_1.5_FDR05 = sum(
          tt_brtn$adj.P.Val < 0.05 & abs(tt_brtn$logFC) >= log2(1.5),
          na.rm = TRUE
        )
      )
      write_tsv_safe(summary_brtn, "H4_HTG_BrighTNess_resumen.tsv")
      print(summary_brtn)
      
    } else {
      warning("HTG BrighTNess no tiene suficientes muestras R/NR para limma (tras emparejar IDs).")
    }
  }
  
} else {
  message("Objeto 'brtn.htg.Rseq' no encontrado en htg_env; se omite análisis HTG BrighTNess.")
}


# ---------------------------------------------------------------------------
# Nota sobre GSE241876
# ---------------------------------------------------------------------------

# GSE241876 (RNA-seq) se usará únicamente como validación exploratoria en el
# script de ssGSEA/GSVA (debido al bajo número de NR, n = 2). Por tanto, aquí
# no ejecutamos limma para esa cohorte.

message("Análisis diferencial COMPLETADO para I-SPY2, HTG MEDI4736 y HTG BrighTNess (limma / limma directo HTG).")





