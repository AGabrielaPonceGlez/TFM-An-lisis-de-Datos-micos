# ============================================================
# 01_setup.R – Configuración del entorno reproducible del TFM
# Proyecto: Análisis transcriptómico para la identificación de patrones de
# respuesta a inmunoterapia en cáncer de mama triple negativo.
# ============================================================

# En este primer Script vamos a preparar el entorno de trabajo para el análisis.
# De esta forma se garantiza la reproducibilidad de proyecto. Para ello,
# configuramos el entorno con "renv", instalamos los paquetes necesarios de CRAN
# y Bioconductor, y creamos la estructura de carpetas de todo el análisis:


# Inicialización de "renv":
if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")
if (!dir.exists("renv")) {
  renv::init(bare = TRUE)
} else {renv::activate()}

# A continuación, configuramos BiocManager y los repositorios de Bioconductor:
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
options(repos = BiocManager::repositories())

# En Windows --> binarios (acelerar la instalación):
if (tolower(Sys.info()[["sysname"]]) == "windows") options(pkgType = "binary")


# =========================================
# Instalación y carga de paquetes de CRAN
# =========================================

paquetes_cran <- c(
  "tidyverse",
  "data.table",
  "readr",
  "readxl",
  "ggplot2",
  "ROCR",
  "msigdbr", "pheatmap")

for (p in paquetes_cran) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  suppressPackageStartupMessages(library(p, character.only = TRUE))}

# ===============================================
# Instalación y carga de paquetes de Bioconductor
# ===============================================

bioc_base <- c("S4Vectors", "IRanges", "GenomeInfoDb")

for (pkg in bioc_base) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = TRUE)}
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))}

# Paquetes para objetos de expresión y análisis exploratorio:
bioc_obj <- c("SummarizedExperiment", "POMA")

for (pkg in bioc_obj) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = TRUE)}
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))}

# Paquete para análisis funcional (ssGSEA / GSVA):
bioc_func <- c("GSVA")

for (pkg in bioc_func) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = TRUE)}
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))}

# ====================================
# Estructura de carpetas del proyecto
# ====================================

# Ahora vamos a definir las rutas del análisis creando carpetas para organizar
# todos los archivos:

path_data    <- file.path("data")
path_scripts <- file.path("scripts")
path_results <- file.path("results")
path_figs    <- file.path("results", "figuras")
path_tabs    <- file.path("results", "tablas")
path_docs    <- file.path("docs")

dirs <- c(path_data, path_scripts, path_results, path_figs, path_tabs, path_docs)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# ============================================
# Comprobación de archivos de datos esperados
# ============================================

# En este punto verificamos que los ficheros que vamos a usar a lo largo del TFM
# estén:

esperados <- c(
  file.path(path_data, "G9_HTG", "valid_datasets_HTG.RData"),
  file.path(path_data, "G9_HTG", "G9_combined_results_all_genes_24JUL2023.xlsx"),
  file.path(path_data, "GSE173839",
            "GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz"),
  file.path(path_data, "GSE173839",
            "GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz"),
  file.path(path_data, "GSE241876", "GSE241876_raw_count_matrix.csv.gz"),
  file.path(path_data, "GSE241876", "GSE241876_DeseqNormalizedCount.csv.gz"),
  file.path(path_data, "GSE241876", "GSE241876_readme.txt"),
  file.path(path_data, "GSE241876", "GSE241876_clinical_manual.csv"))

faltan <- esperados[!file.exists(esperados)]
if (length(faltan) > 0) {
  warning("Los siguientes archivos esperados no se encontraron en 'data':\n",
          paste(" -", faltan, collapse = "\n"))
} else {
  message("Todos los archivos esperados están disponibles en 'data'.")}

# ==================================================
# Snapshot de renv para garantizar reproducibilidad
# ==================================================

set.seed(123)
try(renv::snapshot(prompt = FALSE, force = TRUE), silent = TRUE)

message("Setup COMPLETADO. 'renv.lock' actualizado.")


