# =======================================================================
# 02_curacion_datos.R (Carga, curación y diccionario de metadatos)
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.

# En este segundo script vamos a llevar a cabo toda la curación e integración de
# los metadatos de las distintas cohortes, la generación de un diccionario de
# metadatos, y resúmenes descriptivos:
# =======================================================================



# Cargamos las librerías necesarias:
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(readxl)
  library(SummarizedExperiment)
  library(POMA)})


# Definimos las rutas principales:
path_data    <- "data"
path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
dir.create(path_tabs, showWarnings = FALSE, recursive = TRUE)


# A continuación, hacemos uso de los helpers. El primero guarda las tablas en TSV;
# mientras que el segundo estandariza los identificadores:
write_tsv_safe <- function(df, path) {
  tryCatch({ readr::write_tsv(df, path); message("✔️ Guardado: ", path) },
           error = function(e) warning("No se pudo guardar ", path, " -> ", e$message))}
to_upper_trim <- function(x) { x <- as.character(x); x <- trimws(x); toupper(x) }


# =======================================================================
# 1) COHORTES HTG (MEDI4736, BrighTNess, SCAN-B) --> *datos ya procesados*
# =======================================================================

# En este primer bloque vamos a cargar los objetos procesados del paquete HTG que
# incluyen tres cohortes diferentes: BrighTNess, MEDI4736 y SCAN-B.
# Cada una de ellas contiene información clínica y de respuesta que armonizamos
# bajo un mismo formato para integrarlas:
path_htg_rdata <- file.path(path_data, "G9_HTG", "valid_datasets_HTG.RData")
htg_loaded <- FALSE; htg_env <- new.env()

# Verificamos la existencia del archivo y lo cargamos en un entorno temporal:
if (file.exists(path_htg_rdata)) {
  load(path_htg_rdata, envir = htg_env)
  htg_loaded <- TRUE
  message("Cargado HTG RData: ", path_htg_rdata)
  message("Objetos HTG disponibles: ", paste(ls(htg_env), collapse = ", "))
} else {warning("No se encontró: ", path_htg_rdata)}

# Integramos las tres cohortes:
htg_meta_h <- NULL
if (htg_loaded) {
  combo <- list()
  
  
# BrighTNess. En esta cohorte extraemos los metadatos clínicos:

  if (exists("brtn.htg.cli.dat", envir = htg_env)) {
    df <- get("brtn.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG BrighTNess: ", nrow(df), " muestras")
    combo[["BrighTNess"]] <- df %>%
      transmute(
        cohort     = "HTG_BrighTNess",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(Sample),
        arm        = to_upper_trim(planned_arm),
        response   = case_when(
          pathologic_complete_response %in% c("1","PCR","R","RESPONDER","YES","SI","SÍ") ~ "R",
          pathologic_complete_response %in% c("0","NO_PCR","NONRESPONDER","NR","NO","NO RESPONDER","-1") ~ "NR",
          TRUE ~ NA_character_
        ),
        timepoint = NA_character_)}
  

# MEDI4736 (HTG). Realizamos el mismo procedimiento:

  if (exists("medi.htg.cli.dat", envir = htg_env)) {
    df <- get("medi.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG MEDI4736: ", nrow(df), " muestras")
    combo[["MEDI4736_HTG"]] <- df %>%
      transmute(
        cohort     = "HTG_MEDI4736",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(sample_name),
        arm        = to_upper_trim(treatment),
        response   = case_when(
          pathologic.complete.response %in% c("1","PCR","R","RESPONDER","YES","SI","SÍ") ~ "R",
          pathologic.complete.response %in% c("0","NO_PCR","NONRESPONDER","NR","NO","NO RESPONDER","-1") ~ "NR",
          TRUE ~ NA_character_
        ),
        timepoint = NA_character_)}
  
  
# SCAN-B. inalmente, cargamos la cohorte SCAN-B:
  
  if (exists("scanb.htg.cli.dat", envir = htg_env)) {
    df <- get("scanb.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG SCAN-B: ", nrow(df), " muestras")
    combo[["SCANB"]] <- df %>%
      transmute(
        cohort     = "HTG_SCANB",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(SAMPLE),
        arm        = NA_character_,
        response   = NA_character_,
        timepoint  = NA_character_)}
  
# Una vez procesadas, combinamos todas las cohortes en una tabla única de metadatos HTG:
  if (length(combo)) {
    htg_meta_h <- dplyr::bind_rows(combo) %>%
      filter(!is.na(sample_id) & sample_id != "")
    
    message("Metadatos combinados HTG: ", nrow(htg_meta_h), " muestras totales.")
    print(htg_meta_h %>% group_by(cohort, response) %>% summarise(n = n(), .groups = "drop"))} else {
    message("No se identificaron cohortes HTG con metadatos válidos.")}}

# =======================================================================
# 2) GSE173839 (I-SPY2, Agilent) – matriz + metadatos
# =======================================================================

# Aquí cargamos la información correspondiente al estudio I-SPY2:

path_ispy_expr <- file.path(path_data, "GSE173839",
                            "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz")
path_ispy_meta <- file.path(path_data, "GSE173839",
                            "GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz")

ispy_expr <- ispy_meta <- NULL


# Cargamos la matriz de expresión y mostramos un fragmento para comprobar su estructura:
if (file.exists(path_ispy_expr)) {
  ispy_expr <- data.table::fread(path_ispy_expr, sep = "\t", header = TRUE) |> as_tibble()
  message("I-SPY2 matriz (Agilent) cargada. Dim: ", paste(dim(ispy_expr), collapse = " x "))
  print(ispy_expr %>% select(1:min(6, ncol(.))) %>% head(3))
} else warning("No se encontró matriz I-SPY2: ", path_ispy_expr)

# # Cargamos los metadatos clínicos:
if (file.exists(path_ispy_meta)) {
  ispy_meta <- data.table::fread(path_ispy_meta) |> as_tibble()
  message("I-SPY2 metadatos cargados. Dim: ", paste(dim(ispy_meta), collapse = " x "))
  print(ispy_meta %>% select(1:min(8, ncol(.))) %>% head(5))
} else warning("No se encontró metadatos I-SPY2: ", path_ispy_meta)

# Ahora vamos a identificar dinámicamente las columnas equivalentes a “sample_id”,
# “arm”, “response” y “timepoint” (estandarización):
ispy_meta_h <- NULL
if (!is.null(ispy_meta)) {
  cand_id   <- intersect(names(ispy_meta), c("sample_id","Sample_ID","ID","GEX_ID","ArrayID","array_id","ResearchID"))
  cand_arm  <- intersect(names(ispy_meta), c("treatment","arm","Arm","Treatment","Group","ARM"))
  cand_resp <- intersect(names(ispy_meta), c("pCR","PCR","response","Response","R_NR","resp","pCR.status"))
  cand_tp   <- intersect(names(ispy_meta), c("timepoint","Timepoint","biopsy_time","Visit"))

  # Transformamos y armonizamos las variables clínicas clave:
  ispy_meta_h <- ispy_meta %>%
    mutate(
      sample_id = if (length(cand_id)) .[[cand_id[1]]] else NA,
      arm       = if (length(cand_arm)) .[[cand_arm[1]]] else NA,
      response  = if (length(cand_resp)) .[[cand_resp[1]]] else NA,
      timepoint = if (length(cand_tp)) .[[cand_tp[1]]] else NA
    ) %>%
    transmute(
      cohort = "GSE173839_ISPY2",
      technology = "Agilent_microarray",
      sample_id = to_upper_trim(sample_id),
      arm       = ifelse(is.na(arm), NA, to_upper_trim(arm)),
      response  = ifelse(is.na(response), NA, toupper(as.character(response))),
      timepoint = ifelse(is.na(timepoint), NA, as.character(timepoint))
    ) %>%
    mutate(response = case_when(
      response %in% c("1","PCR","R","RESPONDER","YES","SI","SÍ") ~ "R",
      response %in% c("0","NO_PCR","NONRESPONDER","NR","NO","NO RESPONDER","-1") ~ "NR",
      TRUE ~ response
    ))
}

# =======================================================================
# 3) GSE241876 (MEDI4736, RNA-seq) – crudo + DESeq normalizado
# =======================================================================
# Finalmente, procesamos la cohorte GSE241876:

path_medi_raw   <- file.path(path_data, "GSE241876", "GSE241876_raw_count_matrix.csv.gz")
path_medi_norm  <- file.path(path_data, "GSE241876", "GSE241876_DeseqNormalizedCount.csv.gz")
path_medi_readme<- file.path(path_data, "GSE241876", "GSE241876_readme.txt")

medi_raw  <- medi_norm <- NULL

# Cargamos las matrices de expresión (raw y normalizada):
if (file.exists(path_medi_raw)) {
  medi_raw <- data.table::fread(path_medi_raw) |> as_tibble()
  message("MEDI4736 raw cargado. Dim: ", paste(dim(medi_raw), collapse = " x "))
  print(medi_raw %>% select(1:min(6, ncol(.))) %>% head(3))
} else warning("No se encontró MEDI4736 raw: ", path_medi_raw)

if (file.exists(path_medi_norm)) {
  medi_norm <- data.table::fread(path_medi_norm) |> as_tibble()
  message("MEDI4736 DESeq2 cargado. Dim: ", paste(dim(medi_norm), collapse = " x "))
  print(medi_norm %>% select(1:min(6, ncol(.))) %>% head(3))
} else warning("No se encontró MEDI4736 norm: ", path_medi_norm)

if (file.exists(path_medi_readme)) message("Readme MEDI4736: ", path_medi_readme)

# Eliminamos columnas de anotaciones génicas para quedarnos solo con las muestras:
drop_annot_cols <- function(tbl) {
  if (is.null(tbl)) return(tbl)
  bad <- c("ENSEMBLEID","ENTREZID","GENESYMBOL","SYMBOL","GENE","DESCRIPTION")
  keep <- !(toupper(names(tbl)) %in% bad)
  tbl[, keep]
}
medi_raw  <- drop_annot_cols(medi_raw)
medi_norm <- drop_annot_cols(medi_norm)

# Si no contamos con metadatos clínicos, construimos un esqueleto básico:
extract_sample_ids <- function(expr_tbl) {
  if (is.null(expr_tbl)) return(NULL)
  if (ncol(expr_tbl) <= 1) return(NULL)  # asumimos col1 = gene
  colnames(expr_tbl)[-1]
}
medi_ids <- extract_sample_ids(medi_norm); if (is.null(medi_ids)) medi_ids <- extract_sample_ids(medi_raw)
medi_meta_h <- if (!is.null(medi_ids)) tibble(
  cohort = "GSE241876_MEDI4736",
  technology = "RNAseq",
  sample_id = to_upper_trim(medi_ids),
  arm = NA_character_, response = NA_character_, timepoint = NA_character_
) else NULL

# =======================================================================
# 4) TABLA MAESTRA + DICCIONARIO + RESUMEN
# =======================================================================
# Por último, unificamos todos los metadatos en una tabla maestra ("meta_master"),
# creamos un diccionario de las variables y guardamos los resultados en formato TSV
# dentro de la carpeta "results/tablas":
meta_list <- Filter(Negate(is.null), list(ispy_meta_h, medi_meta_h, htg_meta_h))

if (length(meta_list) == 0) {
  meta_master <- tibble::tibble()
} else {
  meta_master <- dplyr::bind_rows(meta_list) |>
    dplyr::distinct(cohort, technology, sample_id, .keep_all = TRUE) |>
    dplyr::arrange(cohort, technology, sample_id) |>
    tibble::as_tibble()
}

message("Dim meta_master: ", paste(dim(meta_master), collapse = " x "))
print(utils::head(meta_master, 10))
write_tsv_safe(meta_master, file.path(path_tabs, "meta_master.tsv"))

dict_meta <- tibble::tibble(
  campo = c("cohort","technology","sample_id","arm","response","timepoint"),
  descripcion = c(
    "Cohorte/estudio de origen",
    "Tecnología (RNAseq / Agilent_microarray / HTG_EdgeSeq_OBP)",
    "Identificador de muestra (armonizado en mayúsculas y sin espacios)",
    "Brazo de tratamiento (si aplica)",
    "Respuesta: R/NR (pCR u otros criterios, según cohorte)",
    "Momento de biopsia / visita (si aplica)"
  )
)
write_tsv_safe(dict_meta, file.path(path_tabs, "diccionario_metadatos.tsv"))

# Resumen por cohorte/tecnología/response:
if (nrow(meta_master) > 0) {
  resumen1 <- meta_master |>
    dplyr::count(cohort, technology, response, name = "n_muestras") |>
    dplyr::arrange(cohort, technology, dplyr::desc(n_muestras))
  write_tsv_safe(resumen1, file.path(path_tabs, "resumen_cohorte_tecnologia_response.tsv"))
  print(resumen1)
}

message("Carga y curación inicial COMPLETADA.")


