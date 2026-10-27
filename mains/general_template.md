# general_template.qmd — notas

Pipeline de proteómica LFQ independiente del dataset. Solo hay que editar la sección **Configuration** (tres chunks).

## Qué hay en el template

- **Configuración en tres bloques al inicio**
  - Proyecto y datos: `dataset_name`, `organism_key` (human, mouse, zebrafish, rat), `software` (`diann` o `maxquant`), `prot_file`, `pep_file`, `min_peptides`.
  - Diseño: `design_file`.
  - Decisiones: `analysis_mode`, `min_valid_fraction`, `norm_method`, `imputation_method`, `imputation_params`.
- **Diseño en xlsx**
  - Hoja `samples`: `label` y `condition` obligatorias; `replicate` y `batch` opcionales.
  - Hoja `comparisons`: `group1` y `group2`. Si falta, se generan todas las comparaciones por pares.
  - `label` puede ser un sufijo único de la columna, así que no hace falta copiar el prefijo de ruta.
  - Valida duplicados, mínimo de 2 réplicas por condición y que ninguna condición contenga `_vs_`.
  - Si el xlsx no existe, escribe una plantilla y detiene el render.
- **Modo `exploratory` / `production`**
  - Exploratorio: calcula todos los métodos de normalización aplicables y los 5 de imputación, con sus plots comparativos.
  - Producción: calcula solo `norm_method` e `imputation_method`.
  - Los plots comparativos de imputación solo se ejecutan en exploratorio (`#| eval: !expr explore`).
  - Los plots de normalización escalan su altura con el número de assays.
- **Carga de datos**
  - Zero a NA solo en columnas de intensidad.
  - Mapa de columnas por software (DIA-NN / MaxQuant); MaxQuant elimina contaminantes, reversos y site-only.
  - Filtro de péptidos opcional, sin `pivot_longer`.
- **Estilo**: `flatly`, `toc` lateral, `code-fold: show`, `embed-resources`, referencias cruzadas `@sec-`.

## Código de apoyo: `code/prone_helpers.R`

- `as_prone_se()`: adapta el SE a PRONE (assays `raw` y `log2` sin rownames, `Protein.IDs`, `Column`, metadata de condición, batch y refs).
- `prone_select_methods()`: quita `limBE` si hay menos de 2 batches e `IRS` si faltan referencias TMT.
- `prone_valid_assays()`: sustituye el bucle de depuración de `prcomp`.
- `run_prone_normalization()`: envuelve `as_prone_se()` y `PRONE::normalize_se()`.
- En `qc_proteomics.R` solo cambió `normalize_prone()`, que ahora delega en estos helpers.

## Decisiones a tener en cuenta

- **Config en chunks, no `params:`**: permite vectores y listas (`imputation_params`). Se puede migrar a `quarto render -P` si se prefiere.
- **`output_path` termina en `/`**: `06_GO.R` aún usa `paste0(output_path, "tables/...")`. Quitar cuando se refactoricen esos scripts.
- **Orden de secciones**: "Experiment Design" va antes de "Importing Data", porque las columnas de muestra salen del xlsx.
- **Defaults**: `norm_method = "Median"` e `imputation_method = "knn"` son defaults de plantilla; el MAD del análisis anterior era específico de aquel dataset.
- **Otras correcciones menores**
  - `plot_sample_numbers <- plot_sample_numbers(...)` sobrescribía la función; la variable se renombró.
  - Typo `proreins` corregido.
  - Eliminados `rm(list = ls())` y el reload del RDS de normalización.
- **`global_variables.R`** sigue fijando 4 colores por condición. El `.qmd` los redefine tras leer el diseño; conviene moverlo allí al tocar ese archivo.

## Pendiente (no tocado)

- DE (`de_analysis.R`): `run_de_analysis()` crea `data_results` y `significant_*_data` en el entorno global con `assign()`; debe devolverlos.
- `07_Strings.R` carga `dep_analysis.RData`, que ya no se genera.
- Contrastes: `de_plots.R` parte con `strsplit("_vs_")`.
- Secciones de plots, GO, STRING, GSEA, rbioapi, Venn y summary: copiadas tal cual; `pAdjustMethod = "none"` en GO a revisar.

## Verificación realizada

Con datos sintéticos en R 4.5.3: el código del `.qmd` parsea, y funcionan los helpers de PRONE (excepto `normalize_se()` real), la lectura y validación del diseño, el emparejamiento de columnas, zero a NA y el filtro de péptidos. No se ha renderizado el documento completo (PRONE no estaba instalado y no hay datos reales).
