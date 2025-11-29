# =============================================================================
# 05_analisis_funcional.R – Análisis funcional por cohorte
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# =============================================================================
#
# En este script realizamos:
#
#   PARTE A – I-SPY2 (microarrays Agilent, GSE173839_ISPY2)
#     - Cálculo de actividades de vías Hallmark mediante ssGSEA (GSVA, API nueva).
#     - Modelo limma a nivel de vía (R vs NR) ajustando por brazo de tratamiento.
#
#   PARTE B – MEDI (RNA-seq crudo, GSE241876)
#     - Análisis diferencial génico (limma-voom, R vs NR).
#     - GSEA clásica (fgsea) usando las estadísticas de gen (t).
#
# Flujo general:
#   1) Cargar meta_master y expresión de I-SPY2.
#   2) Cargar colecciones Hallmark de MSigDB (msigdbr).
#   3) PARTE A: ssGSEA + limma sobre scores (I-SPY2).
#   4) PARTE B: limma-voom en GSE241876 + GSEA con fgsea.
# =============================================================================


#### 0) Carga de librerías ####

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(ggplot2)
  library(limma)
  library(edgeR)
  library(GSVA)
  library(msigdbr)
  library(fgsea)
})


#### 1) Rutas y funciones auxiliares ####

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

# Prepara lista con matriz de expresión y vector de genes
# (primera columna = gen, resto = muestras)
prepare_expr <- function(tbl, gene_col = 1) {
  stopifnot(is.data.frame(tbl))
  if (ncol(tbl) < 2) stop("La tabla de expresión no tiene columnas suficientes.")
  genes <- tbl[[gene_col]]
  mat   <- as.matrix(tbl[, -gene_col, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  colnames(mat) <- colnames(tbl)[-gene_col]
  list(expr = mat, genes = genes)
}

# Empareja expresión con meta_master para una cohorte concreta
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
    # Si no hay coincidencias directas, probamos a eliminar sufijos típicos:
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

# Colapsa genes duplicados usando la media de filas
collapse_by_gene <- function(X, genes) {
  ord   <- order(genes)
  genes <- genes[ord]
  X     <- X[ord, , drop = FALSE]
  X2    <- rowsum(X, group = genes) / as.vector(table(genes))
  list(expr = X2, genes = rownames(X2))
}

# Gráfico de barras para vías (usa columnas logFC y adj.P.Val)
plot_pathway_bar <- function(tt, title, top_n = 15) {
  tt2 <- tt %>%
    dplyr::mutate(
      pathway = rownames(.),
      sig = adj.P.Val < 0.05
    ) %>%
    dplyr::arrange(adj.P.Val) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  ggplot(tt2, aes(x = pathway, y = logFC, fill = sig)) +
    geom_col() +
    coord_flip() +
    labs(title = title,
         x = "Vía Hallmark",
         y = "log2(FC) (R vs NR)") +
    theme_bw(base_size = 11)
}


#### 2) Cargar meta_master e I-SPY2 ####

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


#### 3) Obtener colecciones Hallmark de MSigDB ####

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


# =============================================================================
# PARTE A – I-SPY2: ssGSEA (Hallmark) + limma a nivel de vía
# =============================================================================

#### A1) Preparación de expresión I-SPY2 ####

prep_ispy <- prepare_expr(ispy_tbl, gene_col = 1)
X_all     <- prep_ispy$expr
genes_all <- prep_ispy$genes

# Filtro mínimo por varianza (quitamos genes constantes)
v_all <- apply(X_all, 1, stats::var)
keep  <- is.finite(v_all) & v_all > 0
X     <- X_all[keep, , drop = FALSE]
genes <- genes_all[keep]

# Emparejamos con meta_master para la cohorte I-SPY2
matched_ispy <- match_expr_meta(
  list(expr = X, genes = genes),
  meta_master,
  "GSE173839_ISPY2"
)

X     <- matched_ispy$expr
genes <- matched_ispy$genes
pheno_ispy <- matched_ispy$pheno

# Nos quedamos solo con R/NR y preparamos 'response' y 'arm'
pheno_ispy <- pheno_ispy %>%
  dplyr::filter(response %in% c("R", "NR")) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR", "R")),   # coef 'responseR' = R > NR
    arm      = to_upper_trim(arm),
    arm      = forcats::fct_lump_min(arm, min = 5),       # agrupamos brazos minoritarios
    arm      = stats::relevel(factor(arm), ref = "CONTROL")
  )

samples_keep_ispy <- pheno_ispy$sample_id
X <- X[, samples_keep_ispy, drop = FALSE]

pheno_ispy <- pheno_ispy[match(colnames(X), pheno_ispy$sample_id), , drop = FALSE]
stopifnot(identical(pheno_ispy$sample_id, colnames(X)))

# Colapsamos posibles genes duplicados por símbolo
collapsed_ispy <- collapse_by_gene(X, genes)
X_ssg          <- collapsed_ispy$expr
rownames(X_ssg) <- collapsed_ispy$genes

message("I-SPY2: genes tras colapso: ", nrow(X_ssg))
message("I-SPY2: muestras (R/NR) tras emparejar: ", ncol(X_ssg))


#### A2) ssGSEA con GSVA (API nueva) ####

# API nueva de GSVA: construimos un objeto de parámetros con ssgseaParam()
ssgsea_param_ispy <- GSVA::ssgseaParam(
  exprData = X_ssg,
  geneSets = hallmark_list
  # En tu versión, los argumentos abs.ranking / mx.diff no están,
  # así que usamos los valores por defecto del método ssGSEA.
)

scores_ispy <- GSVA::gsva(
  ssgsea_param_ispy,
  verbose = FALSE
)

message("Dim scores ssGSEA I-SPY2 (vías x muestras): ",
        paste(dim(scores_ispy), collapse = " x "))

scores_ispy_df <- as.data.frame(scores_ispy) %>%
  tibble::rownames_to_column("pathway")

write_tsv_safe(scores_ispy_df, "H5_ISPY2_ssGSEA_scores.tsv")


#### A3) Modelo limma a nivel de vía (R vs NR, ajustando por 'arm') ####

design_ispy <- stats::model.matrix(~ response + arm, data = pheno_ispy)

fit_ispy <- limma::lmFit(scores_ispy, design_ispy)
fit_ispy <- limma::eBayes(fit_ispy, trend = TRUE, robust = TRUE)

coef_name <- "responseR"
if (!(coef_name %in% colnames(coef(fit_ispy)))) {
  stop("No se encontró el coeficiente 'responseR' en el modelo para I-SPY2 (vías).")
}

tt_path_ispy <- limma::topTable(
  fit_ispy,
  coef   = coef_name,
  number = Inf,
  sort.by = "P"
)

# Guardamos tabla de vías (añadimos columna explícita con el nombre de la vía)
tt_path_ispy_out <- tt_path_ispy %>%
  tibble::rownames_to_column("pathway")

write_tsv_safe(tt_path_ispy_out, "H5_ISPY2_ssGSEA_R_vs_NR_pathways.tsv")
write_tsv_safe(tt_path_ispy_out %>% dplyr::slice_head(n = 30),
               "H5_ISPY2_ssGSEA_R_vs_NR_pathways_top30.tsv")

summary_ispy_path <- tibble::tibble(
  cohort                 = "GSE173839_ISPY2",
  n_pathways             = nrow(tt_path_ispy),
  n_sig_fdr_05           = sum(tt_path_ispy$adj.P.Val < 0.05),
  n_sig_absLFC_0.2_FDR05 = sum(
    tt_path_ispy$adj.P.Val < 0.05 &
      abs(tt_path_ispy$logFC) >= 0.2
  )
)

write_tsv_safe(summary_ispy_path, "H5_ISPY2_ssGSEA_resumen.tsv")
print(summary_ispy_path)

# Figura: barplot de las top vías por logFC
p_bar_ispy <- plot_pathway_bar(
  tt_path_ispy,
  "I-SPY2 – ssGSEA Hallmark: vías asociadas a R vs NR (ajustado por brazo)",
  top_n = 15
)

save_plot(p_bar_ispy, "H5_ISPY2_ssGSEA_pathways_barplot_top15.png")

message("PARTE A (I-SPY2 ssGSEA + limma vías) COMPLETADA.")


# =============================================================================
# PARTE B – GSE241876 (MEDI): limma-voom + GSEA clásica (fgsea)
# =============================================================================

message("=== Iniciando análisis funcional para GSE241876 (MEDI) ===")

#### B1) Cargar matriz de expresión cruda (raw counts) ####

gse_tbl <- data.table::fread(
  file.path("data", "GSE241876", "GSE241876_raw_count_matrix.csv.gz")
) %>% tibble::as_tibble()

message("GSE241876 raw cargado. Dim: ", paste(dim(gse_tbl), collapse = " x "))

# Nos quedamos con símbolo de gen + columnas de muestras (R020, R022, etc.)
gse_expr_tbl <- gse_tbl %>%
  dplyr::select(
    GeneSymbol,
    dplyr::starts_with("R")
  )

message("GSE241876 – expresión (GeneSymbol + muestras). Dim: ",
        paste(dim(gse_expr_tbl), collapse = " x "))


#### B2) Construir matriz genes x muestras y filtrar ####

genes_g_all <- gse_expr_tbl$GeneSymbol
X_gse_all   <- as.matrix(gse_expr_tbl[, -1, drop = FALSE])
mode(X_gse_all) <- "numeric"
rownames(X_gse_all) <- genes_g_all

message("GSE241876 – matriz inicial: ",
        paste(dim(X_gse_all), collapse = " x "))

# Filtro por varianza > 0
v_gse      <- apply(X_gse_all, 1, stats::var)
keep_var   <- is.finite(v_gse) & v_gse > 0
X_gse      <- X_gse_all[keep_var, , drop = FALSE]
genes_g    <- genes_g_all[keep_var]

message("GSE241876 – genes tras varianza > 0: ", nrow(X_gse))

# Eliminamos genes sin símbolo (NA) antes de colapsar
keep_not_na <- !is.na(genes_g)
X_gse       <- X_gse[keep_not_na, , drop = FALSE]
genes_g     <- genes_g[keep_not_na]

message("GSE241876 – genes tras eliminar NA: ", nrow(X_gse))

# Colapsamos genes duplicados por símbolo
collapsed_gse <- collapse_by_gene(X_gse, genes_g)
X_gse_coll    <- collapsed_gse$expr
genes_g_coll  <- collapsed_gse$genes
rownames(X_gse_coll) <- genes_g_coll

message("GSE241876 – genes tras colapso por símbolo: ", nrow(X_gse_coll))
message("GSE241876 – muestras (todas): ", ncol(X_gse_coll))


#### B3) Emparejar con meta_master y filtrar por R/NR ####

samples_expr_gse <- colnames(X_gse_coll)

pheno_gse <- meta_master %>%
  dplyr::filter(
    cohort    == "GSE241876",
    sample_id %in% samples_expr_gse,
    response  %in% c("R", "NR")
  ) %>%
  dplyr::distinct(sample_id, .keep_all = TRUE) %>%
  dplyr::mutate(
    response = factor(response, levels = c("NR", "R"))
  )

samples_keep_gse <- pheno_gse$sample_id
X_gse_coll       <- X_gse_coll[, samples_keep_gse, drop = FALSE]

pheno_gse <- pheno_gse[match(colnames(X_gse_coll), pheno_gse$sample_id), , drop = FALSE]
stopifnot(identical(pheno_gse$sample_id, colnames(X_gse_coll)))

message("GSE241876 – muestras con R/NR: ", ncol(X_gse_coll))
print(table(pheno_gse$response))


#### B4) limma-voom (análisis diferencial génico) ####

dge_gse <- edgeR::DGEList(
  counts = X_gse_coll,
  genes  = data.frame(Gene = rownames(X_gse_coll))
)

keep_edge <- edgeR::filterByExpr(dge_gse, group = pheno_gse$response)
dge_gse   <- dge_gse[keep_edge, , keep.lib.sizes = FALSE]
dge_gse   <- edgeR::calcNormFactors(dge_gse, method = "TMM")

message("GSE241876 – genes tras filterByExpr: ", nrow(dge_gse))

design_gse <- stats::model.matrix(~ response, data = pheno_gse)

v_gse_voom <- limma::voom(dge_gse, design_gse, plot = FALSE)
fit_gse    <- limma::lmFit(v_gse_voom, design_gse)
fit_gse    <- limma::eBayes(fit_gse, trend = TRUE, robust = TRUE)

coef_name_gse <- "responseR"
if (!(coef_name_gse %in% colnames(coef(fit_gse)))) {
  stop("No se encontró el coeficiente 'responseR' en el modelo para GSE241876.")
}

tt_gse <- limma::topTable(
  fit_gse,
  coef   = coef_name_gse,
  number = Inf,
  sort.by = "P"
)

# Aseguramos que exista columna 'Gene' (la debería crear topTable a partir de dge_gse$genes)
if (!("Gene" %in% colnames(tt_gse))) {
  tt_gse <- tt_gse %>% tibble::rownames_to_column("Gene")
}

write_tsv_safe(tt_gse, "H5_GSE241876_limma_R_vs_NR.tsv")
write_tsv_safe(tt_gse %>% dplyr::slice_head(n = 50),
               "H5_GSE241876_limma_R_vs_NR_top50.tsv")

summary_gse <- tibble::tibble(
  cohort                 = "GSE241876",
  n_genes                = nrow(tt_gse),
  n_sig_fdr_05           = sum(tt_gse$adj.P.Val < 0.05, na.rm = TRUE),
  n_sig_absLFC_1.5_FDR05 = sum(
    tt_gse$adj.P.Val < 0.05 &
      abs(tt_gse$logFC) >= log2(1.5),
    na.rm = TRUE
  )
)

write_tsv_safe(summary_gse, "H5_GSE241876_limma_resumen.tsv")
print(summary_gse)


#### B5) GSEA clásica (fgsea) con estadísticas de gen ####

# Vector de estadísticas (t-statistic) nombrado por símbolo de gen
stats_gse <- tt_gse$t
names(stats_gse) <- tt_gse$Gene

# Quitamos genes sin nombre y ordenamos
stats_gse <- stats_gse[!is.na(names(stats_gse))]
stats_gse <- sort(stats_gse, decreasing = TRUE)

fgsea_res <- fgsea::fgsea(
  pathways = hallmark_list,
  stats    = stats_gse,
  minSize  = 10,
  maxSize  = 500,
  nperm    = 10000
)

# Ordenamos por FDR (padj)
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

fgsea_res_tbl <- fgsea_res %>% as_tibble()
write_tsv_safe(fgsea_res_tbl, "H5_GSE241876_fgsea_Hallmark_R_vs_NR.tsv")

summary_gsea <- tibble::tibble(
  cohort                  = "GSE241876",
  n_pathways              = nrow(fgsea_res_tbl),
  n_sig_fdr_05            = sum(fgsea_res_tbl$padj < 0.05, na.rm = TRUE),
  n_sig_absNES_1.5_FDR05  = sum(
    fgsea_res_tbl$padj < 0.05 &
      abs(fgsea_res_tbl$NES) >= 1.5,
    na.rm = TRUE
  )
)

write_tsv_safe(summary_gsea, "H5_GSE241876_fgsea_resumen.tsv")
print(summary_gsea)


#### B6) Barplot de vías top según fgsea (NES) ####

# Reutilizamos plot_pathway_bar "engañándolo" con NES y padj
tt_gsea_path <- fgsea_res_tbl %>%
  dplyr::select(pathway, NES, padj) %>%
  dplyr::rename(logFC = NES, adj.P.Val = padj) %>%
  tibble::column_to_rownames("pathway")

p_bar_gse <- plot_pathway_bar(
  tt_gsea_path,
  "GSE241876 – GSEA Hallmark: vías asociadas a R vs NR",
  top_n = 15
)

save_plot(p_bar_gse, "H5_GSE241876_fgsea_Hallmark_barplot_top15.png")

message("=== Script 05 (I-SPY2 ssGSEA + GSE241876 limma-voom + GSEA) COMPLETADO ===")



