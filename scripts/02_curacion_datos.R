# =======================================================================
# 02_curacion_datos.R (Carga, curación y diccionario de metadatos)
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# =======================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(readr)
  library(readxl)
  library(SummarizedExperiment)
  library(POMA)
})

# -----------------------------
# Rutas principales y helpers
# -----------------------------
path_data    <- "data"
path_results <- "results"
path_tabs    <- file.path(path_results, "tablas")
dir.create(path_tabs, showWarnings = FALSE, recursive = TRUE)

write_tsv_safe <- function(df, path) {
  tryCatch({
    readr::write_tsv(df, path)
    message("Guardado: ", path)
  }, error = function(e) {
    warning("No se pudo guardar ", path, " -> ", e$message)
  })
}

to_upper_trim <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  toupper(x)
}

# =======================================================================
# 1) COHORTES HTG (BrighTNess, MEDI4736, SCAN-B)
# =======================================================================

path_htg_rdata <- file.path(path_data, "G9_HTG", "valid_datasets_HTG.RData")
htg_env   <- new.env()
htg_meta_h <- NULL

if (file.exists(path_htg_rdata)) {
  load(path_htg_rdata, envir = htg_env)
  message("HTG cargado desde: ", path_htg_rdata)
  
  combo <- list()
  
  # ---------------------------
  # HTG BrighTNess
  # ---------------------------
  if (exists("brtn.htg.cli.dat", envir = htg_env)) {
    df <- get("brtn.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG BrighTNess: ", nrow(df), " muestras")
    
    combo[["HTG_BrighTNess"]] <- df %>%
      mutate(
        resp_raw = toupper(trimws(pathologic_complete_response))
      ) %>%
      transmute(
        cohort     = "HTG_BrighTNess",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(Sample),
        arm        = to_upper_trim(planned_arm),
        response   = case_when(
          resp_raw == "PCR" ~ "R",
          resp_raw %in% c("RD", "NON-PCR", "NON PCR", "NO PCR") ~ "NR",
          TRUE ~ NA_character_
        ),
        timepoint  = NA_character_
      )
  }
  
  # ---------------------------
  # HTG MEDI4736
  # ---------------------------
  if (exists("medi.htg.cli.dat", envir = htg_env)) {
    df <- get("medi.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG MEDI4736: ", nrow(df), " muestras")
    
    combo[["HTG_MEDI4736"]] <- df %>%
      mutate(
        resp_raw = toupper(trimws(pathologic.complete.response))
      ) %>%
      transmute(
        cohort     = "HTG_MEDI4736",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(sample_name),
        arm        = to_upper_trim(treatment),
        response   = case_when(
          resp_raw %in% c("YES","Y","1","PCR","R","RESPONDER") ~ "R",
          resp_raw %in% c("NO","N","0","NR","NONRESPONDER")   ~ "NR",
          TRUE ~ NA_character_
        ),
        timepoint  = NA_character_
      )
  }
  
  # ---------------------------
  # HTG SCAN-B (sin R/NR; para supervivencia)
  # ---------------------------
  if (exists("scanb.htg.cli.dat", envir = htg_env)) {
    df <- get("scanb.htg.cli.dat", envir = htg_env) %>% as_tibble()
    message("Cohorte HTG SCAN-B: ", nrow(df), " muestras")
    
    combo[["HTG_SCANB"]] <- df %>%
      transmute(
        cohort     = "HTG_SCANB",
        technology = "HTG_EdgeSeq_OBP",
        sample_id  = to_upper_trim(SAMPLE),
        arm        = NA_character_,
        response   = NA_character_,
        timepoint  = NA_character_
      )
  }
  
  combo_df <- Filter(\(x) is.data.frame(x) || tibble::is_tibble(x), combo)
  
  if (length(combo_df) > 0) {
    htg_meta_h <- dplyr::bind_rows(combo_df) %>%
      tibble::as_tibble() %>%
      dplyr::filter(!is.na(sample_id) & sample_id != "")
    
    message("Metadatos combinados HTG: ", nrow(htg_meta_h), " muestras totales.")
    htg_meta_h %>%
      dplyr::count(cohort, response) %>%
      print()
  } else {
    message("No se identificaron cohortes HTG con metadatos válidos.")
  }
  
} else {
  warning("No se encontró archivo HTG: ", path_htg_rdata)
}

# =======================================================================
# 2) GSE173839 (I-SPY2, Agilent) – expresión + metadatos clínicos
# =======================================================================

path_ispy_expr <- file.path(
  path_data, "GSE173839",
  "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"
)
path_ispy_meta <- file.path(
  path_data, "GSE173839",
  "GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz"
)

ispy_expr <- ispy_meta <- NULL

# Matriz de expresión
if (file.exists(path_ispy_expr)) {
  ispy_expr <- data.table::fread(path_ispy_expr, sep = "\t", header = TRUE) %>%
    tibble::as_tibble()
  message("I-SPY2 matriz (Agilent) cargada. Dim: ", paste(dim(ispy_expr), collapse = " x "))
  print(ispy_expr %>% dplyr::select(1:min(6, ncol(.))) %>% head(3))
} else {
  warning("No se encontró matriz I-SPY2: ", path_ispy_expr)
}

# Metadatos clínicos I-SPY2
if (file.exists(path_ispy_meta)) {
  ispy_meta <- data.table::fread(path_ispy_meta) %>% tibble::as_tibble()
  message("I-SPY2 metadatos cargados. Dim: ", paste(dim(ispy_meta), collapse = " x "))
  print(ispy_meta %>% dplyr::select(1:min(8, ncol(.))) %>% head(5))
} else {
  warning("No se encontró metadatos I-SPY2: ", path_ispy_meta)
}

# Estandarización de columnas clave para meta_master
ispy_meta_h <- NULL
if (!is.null(ispy_meta)) {
  cand_id   <- intersect(names(ispy_meta),
                         c("sample_id","Sample_ID","ID","GEX_ID","ArrayID","array_id","ResearchID"))
  cand_arm  <- intersect(names(ispy_meta),
                         c("treatment","arm","Arm","Treatment","Group","ARM"))
  cand_resp <- intersect(names(ispy_meta),
                         c("pCR","PCR","response","Response","R_NR","resp","pCR.status"))
  cand_tp   <- intersect(names(ispy_meta),
                         c("timepoint","Timepoint","biopsy_time","Visit"))
  
  ispy_meta_h <- ispy_meta %>%
    dplyr::mutate(
      sample_id = if (length(cand_id)) .[[cand_id[1]]] else NA,
      arm       = if (length(cand_arm)) .[[cand_arm[1]]] else NA,
      response  = if (length(cand_resp)) .[[cand_resp[1]]] else NA,
      timepoint = if (length(cand_tp)) .[[cand_tp[1]]] else NA
    ) %>%
    dplyr::transmute(
      cohort     = "GSE173839_ISPY2",
      technology = "Agilent_microarray",
      sample_id  = to_upper_trim(sample_id),
      arm        = dplyr::if_else(is.na(arm), NA_character_, to_upper_trim(arm)),
      response   = dplyr::if_else(is.na(response), NA_character_,
                                  toupper(as.character(response))),
      timepoint  = dplyr::if_else(is.na(timepoint), NA_character_,
                                  as.character(timepoint))
    ) %>%
    dplyr::mutate(
      response = dplyr::case_when(
        response %in% c("1","PCR","R","RESPONDER","YES","SI","SÍ") ~ "R",
        response %in% c("0","NO_PCR","NONRESPONDER","NR","NO","NO RESPONDER","-1") ~ "NR",
        TRUE ~ response
      )
    )
}

# =======================================================================
# 3) GSE241876 (RNA-seq) – expresión + clínico manual (15 pacientes evaluables)
# =======================================================================

path_medi_raw    <- file.path(path_data, "GSE241876", "GSE241876_raw_count_matrix.csv.gz")
path_medi_norm   <- file.path(path_data, "GSE241876", "GSE241876_DeseqNormalizedCount.csv.gz")
path_medi_readme <- file.path(path_data, "GSE241876", "GSE241876_readme.txt")
path_medi_clin   <- file.path(path_data, "GSE241876", "GSE241876_clinical_manual.csv")

medi_raw  <- medi_norm <- NULL
medi_clin <- NULL

# --- Cargar matrices de expresión ---
if (file.exists(path_medi_raw)) {
  medi_raw <- data.table::fread(path_medi_raw) |> as_tibble()
  message("GSE241876 raw cargado. Dim: ", paste(dim(medi_raw), collapse = " x "))
}
if (file.exists(path_medi_norm)) {
  medi_norm <- data.table::fread(path_medi_norm) |> as_tibble()
  message("GSE241876 DESeq2 cargado. Dim: ", paste(dim(medi_norm), collapse = " x "))
}
if (file.exists(path_medi_readme)) {
  message("Readme GSE241876: ", path_medi_readme)
}

# --- Cargar CSV clínico (YA LIMPIO: solo 15 pacientes evaluables R/NR) ---
if (file.exists(path_medi_clin)) {
  medi_clin <- readr::read_delim(
    file = path_medi_clin,
    delim = ";",
    show_col_types = FALSE,
    trim_ws = TRUE
  ) |>
    dplyr::mutate(
      patient_id  = to_upper_trim(patient_id),
      response    = toupper(trimws(response)),   # R / NR
      recist      = toupper(trimws(recist)),
      PDL1_status = toupper(trimws(PDL1_status)),
      alive       = toupper(trimws(alive))
    )
  
  message("Metadatos clínicos manuales (GSE241876) cargados. n = ", nrow(medi_clin))
  print(medi_clin)
} else {
  warning("No se encontró el fichero clínico manual: ", path_medi_clin)
}

# --- Limpieza de columnas de anotación ---
drop_annot_cols <- function(tbl) {
  if (is.null(tbl)) return(tbl)
  bad  <- c("ENSEMBLEID","ENTREZID","GENESYMBOL","SYMBOL","GENE","DESCRIPTION","COLUMN1")
  keep <- !(toupper(names(tbl)) %in% bad)
  tbl[, keep]
}

medi_raw  <- drop_annot_cols(medi_raw)
medi_norm <- drop_annot_cols(medi_norm)

# --- Extraer IDs de muestra ---
extract_sample_ids <- function(expr_tbl) {
  if (is.null(expr_tbl)) return(NULL)
  if (ncol(expr_tbl) <= 1) return(NULL)
  colnames(expr_tbl)[-1]
}

medi_ids <- extract_sample_ids(medi_norm)
if (is.null(medi_ids)) medi_ids <- extract_sample_ids(medi_raw)

# --- Crear metadatos armonizados SOLO PARA PACIENTES EVALUABLES (15) ---
medi_meta_h <- NULL

if (!is.null(medi_ids) && !is.null(medi_clin)) {
  
  medi_meta_h <- tibble(
    cohort     = "GSE241876",
    technology = "RNAseq",
    sample_id  = to_upper_trim(medi_ids)
  ) |>
    mutate(
      patient_id = gsub("_(PRE|POST)$", "", sample_id),
      timepoint  = case_when(
        grepl("_PRE$",  sample_id) ~ "PRE",
        grepl("_POST$", sample_id) ~ "POST",
        TRUE ~ NA_character_
      ),
      arm        = "CNP_PEMBROLIZUMAB+CARBOPLATINO+NAB_PACLITAXEL"
    ) |>
    # Unimos con el clínico (15 pacientes evaluables R/NR)
    left_join(medi_clin, by = "patient_id") |>
    # Nos quedamos SOLO con pacientes del CSV y sin NA en response
    filter(
      patient_id %in% medi_clin$patient_id,
      !is.na(response)
    ) |>
    # Priorizamos PRE > POST > NA para quedarnos con 1 muestra por paciente
    arrange(
      patient_id,
      match(timepoint, c("PRE","POST"))
    ) |>
    distinct(patient_id, .keep_all = TRUE) |>
    # Normalizamos response explícitamente a R/NR por seguridad
    mutate(
      response = case_when(
        response %in% c("R","NR") ~ response,
        TRUE ~ NA_character_
      )
    )
  
  message("Metadatos armonizados GSE241876 (solo pacientes evaluables): ",
          nrow(medi_meta_h), " filas.")
  medi_meta_h %>%
    dplyr::count(response, name = "n_pacientes") %>%
    print()
  
} else {
  message("No se pudieron construir metadatos armonizados para GSE241876 (faltan IDs o clínico).")
}

# =======================================================================
# 4) TABLA MAESTRA + DICCIONARIO + RESUMEN
# =======================================================================

meta_list <- Filter(
  Negate(is.null),
  list(ispy_meta_h, medi_meta_h, htg_meta_h)
)

if (length(meta_list) == 0) {
  meta_master <- tibble::tibble()
} else {
  meta_master <- dplyr::bind_rows(meta_list) %>%
    dplyr::distinct(cohort, technology, sample_id, .keep_all = TRUE) %>%
    dplyr::arrange(cohort, technology, sample_id) %>%
    tibble::as_tibble()
}

message("Dim meta_master: ", paste(dim(meta_master), collapse = " x "))
print(utils::head(meta_master, 10))

write_tsv_safe(meta_master, file.path(path_tabs, "meta_master.tsv"))

# Diccionario de metadatos
dict_meta <- tibble::tibble(
  campo = c("cohort","technology","sample_id","arm","response","timepoint"),
  descripcion = c(
    "Cohorte/estudio de origen",
    "Tecnología (RNAseq / Agilent_microarray / HTG_EdgeSeq_OBP)",
    "Identificador de muestra (armonizado en mayúsculas y sin espacios)",
    "Brazo de tratamiento (si aplica)",
    "Respuesta: R/NR/NE (pCR, RECIST u otros criterios, según cohorte)",
    "Momento de biopsia / visita (si aplica; p.ej. PRE/POST tratamiento)"
  )
)

write_tsv_safe(dict_meta, file.path(path_tabs, "diccionario_metadatos.tsv"))

# Resumen por cohorte/tecnología/response
if (nrow(meta_master) > 0) {
  resumen1 <- meta_master %>%
    dplyr::count(cohort, technology, response, name = "n_muestras") %>%
    dplyr::arrange(cohort, technology, dplyr::desc(n_muestras))
  
  write_tsv_safe(resumen1, file.path(path_tabs, "resumen_cohorte_tecnologia_response.tsv"))
  print(resumen1)
}

message("Carga y curación inicial COMPLETADA.")


