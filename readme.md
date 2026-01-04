# Análisis de Datos Ómicos

**Análisis transcriptómico para la identificación de patrones de respuesta a inmunoterapia en cáncer de mama triple negativo**

Este repositorio contiene el flujo completo de análisis desarrollado en el Trabajo de Fin de Máster (TFM), centrado en el análisis transcriptómico de pacientes con cáncer de mama triple negativo (TNBC) tratados con inmunoterapia.

El objetivo del proyecto es explorar diferencias transcriptómicas y funcionales entre pacientes respondedores (R) y no respondedores (NR), a través de datos procedentes de distintas cohortes y tecnologías.

---

## Requisitos

- R (versión ≥ 4.x)
- RStudio (recomendado)
- Sistema operativo: Windows, macOS o Linux

Este proyecto utiliza el paquete **renv** para gestionar un entorno reproducible.

---

## Reproducibilidad del entorno (renv)

### Abrir el proyecto

Se recomienda abrir el archivo de proyecto de RStudio ubicado en la carpeta principal:

```
/TFM.Rproj
```

---

### Activar el entorno reproducible

Desde la raíz del proyecto, ejecutar:

```r
source("renv/activate.R")
```

Si el archivo `renv.lock` está disponible, restaurar las dependencias con:

```r
renv::restore()
```

---

## Datos de entrada

Los datos utilizados proceden de estudios públicos y deben estar disponibles en las rutas indicadas antes de ejecutar los scripts.

---

### I-SPY2 – GSE173839 (Agilent)

**Ubicación:**

```
data/GSE173839/
```

**Archivos esperados:**

- `GSE173839_ISPY2_AgilentGeneExp_durvaPlusCtr_FFPE_meanCol_geneLevel_n105.txt.gz`
- `GSE173839_ISPY2_DurvalumabOlaparibArm_biomarkers.csv.gz`

---

### HTG EdgeSeq

**Ubicación:**

```
data/G9_HTG/
```

**Archivos esperados:**

- `valid_datasets_HTG.RData`
- `G9_combined_results_all_genes_24JUL2023.xlsx`

---

### GSE241876 (RNA-seq)

**Ubicación:**

```
data/GSE241876/
```

**Archivos esperados:**

- `GSE241876_raw_count_matrix.csv.gz`
- `GSE241876_DeseqNormalizedCount.csv.gz`
- `GSE241876_readme.txt`
- `GSE241876_clinical_manual.csv`

---

## Orden de ejecución de los scripts

1. `01_setup.R` – Configuración del entorno
2. `02_curacion_datos.R` – Curación e integración de metadatos
3. `03_exploracion_inicial.R` – Análisis exploratorio y control de calidad
4. `04_analisis_diferencial.R` – Análisis de expresión diferencial
5. `05_analisis_funcional.R` – Análisis funcional
6. `06_validacion_sensibilidad.R` – Validación y análisis de sensibilidad
7. `07_comparacion_cohortes.R` – Integración y comparación multi-cohorte

---

## Resultados

Los resultados generados se almacenan en:

- `results/tablas/`
- `results/figuras/`

Los archivos siguen la convención:

```
H{n}_DESCRIPCION.ext
```

---

## Nota

Este repositorio forma parte del Trabajo de Fin de Máster con fines académicos y reproducibles. Los datos originales pertenecen a sus respectivos estudios y consorcios.

