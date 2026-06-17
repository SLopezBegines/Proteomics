
# Save Plots

# FUNCIÓN AUXILIAR: GUARDAR GRÁFICOS EN MÚLTIPLES FORMATOS


save_plot <- function(plotname,
                      plot,
                      output_dir = output_path,
                      image_number = NULL,
                      units = "in", #or "cm"
                      width = NULL,
                      height = NULL,
                      dpi = 300,
                      formats = c("tiff", "pdf"),
                      print_plot = TRUE,
                      cleanup = FALSE) {
  #' Guarda gráficos ggplot2 en múltiples formatos
  #'
  #' @param plotname Nombre base del archivo (sin extensión)
  #' @param plot Objeto ggplot2
  #' @param output_dir Directorio de salida
  #' @param image_number Número de imagen para prefijo (opcional)
  #' @param width Ancho en cm
  #' @param height Alto en cm
  #' @param dpi Resolución para formatos raster
  #' @param formats Vector de formatos ("tiff", "pdf", "png")
  #' @param print_plot Imprimir plot en consola
  #' @param cleanup Eliminar objeto plot después de guardar
  #' @return Ruta del primer archivo guardado (invisible)
  
  # Verificar objeto ggplot
  if (!inherits(plot, "ggplot")) {
    warning(paste("El objeto", plotname, "no es un ggplot válido. Omitiendo."))
    return(invisible(NULL))
  }
  
  # Construir nombre base del archivo
  # Usar variable global si no se especifica
  if (is.null(image_number)) {
    image_number <- get("image_number", envir = .GlobalEnv)
  }
  
  base_filename <- sprintf("%03d_%s", image_number, plotname)
  
  
  # Guardar en cada formato
  saved_paths <- c()
  
  for (format in formats) {
    extension <- paste0(".", format)
    filepath <- file.path(output_dir, paste0(base_filename, extension))
    
    tryCatch(
      {
        # Construir argumentos condicionalmente
        save_args <- list(
          filename = filepath,
          plot = plot,
          units = units
        )
        # Agregar dimensiones solo si se especifican
        if (!is.null(width)) save_args$width <- width
        if (!is.null(height)) save_args$height <- height
        
        # Solo agregar dpi para formatos raster
        if (format %in% c("tiff", "png")) {
          save_args$dpi <- dpi
        }
        
        do.call(ggsave, save_args)
        saved_paths <- c(saved_paths, filepath)
      },
      error = function(e) {
        warning(paste("No se pudo guardar:", filepath, "- Error:", e$message))
      }
    )
  }
  
  # Incrementar contador global
  assign("image_number", image_number + 1, envir = .GlobalEnv)
  # Imprimir plot si se solicita
  if (print_plot) {
    print(plot)
  }
  
  # Limpiar objeto si se solicita
  if (cleanup && exists("plot", inherits = FALSE)) {
    rm(plot, envir = parent.frame())
  }
  
  return(invisible(saved_paths[1]))
}


# Prepare DE gene sets ####
# Reusable helper to build per-comparison gene-set lists (UP/DOWN/UPDOWN)
# from a `data_results` table. Used as input preparation for GO, STRING,
# gseGO, gseKEGG, etc.

prepare_DE_gene_sets <- function(data_results,
                                  comparisons,
                                  direction = c("UP", "DOWN"),
                                  independent_UPDOWN = TRUE,
                                  id_cols = "ID",
                                  flatten = TRUE,
                                  diffexpressed_type = c("qval", "pval", "padj_holm_pval")) {
  #' @param data_results Data frame produced by data_analysis() (04_data_analysis.R),
  #'   with `<comparison>_diffexpressed_<diffexpressed_type>` columns
  #'   ("UP"/"DOWN"/"NO") for each comparison in `comparisons`.
  #' @param comparisons Character vector of comparison names (e.g. "CTRL_vs_WT").
  #' @param direction Which categories to keep: "UP", "DOWN", or both.
  #' @param independent_UPDOWN If TRUE, UP and DOWN gene sets are returned
  #'   separately (one element per comparison and direction); if FALSE, a
  #'   single combined set per comparison is returned.
  #' @param id_cols Column(s) of `data_results` to keep for each gene set.
  #' @param flatten If TRUE, each one-column gene set is unlisted into a
  #'   plain named vector (input format expected by enrichGO). If FALSE,
  #'   each gene set is kept as a data frame (e.g. for STRINGdb).
  #' @param diffexpressed_type Which significance criterion to use for the
  #'   "diffexpressed" column: "qval" (Storey's q-value, default), "pval"
  #'   (raw p-value), or "padj_holm_pval" (adjusted p-value via `adjusted_method`).
  #' @return A flat named list of gene sets (vectors if `flatten = TRUE`,
  #'   data frames otherwise), named "names_<comparison>_UP"/"_DN"/"_UPDOWN".

  direction <- match.arg(direction, choices = c("UP", "DOWN"), several.ok = TRUE)
  diffexpressed_type <- match.arg(diffexpressed_type)

  gene_sets <- lapply(comparisons, function(comp) {
    diffexpressed_var <- paste0(comp, "_diffexpressed_", diffexpressed_type)

    if (independent_UPDOWN) {
      sets <- lapply(direction, function(dir) {
        sub_df <- data_results %>% dplyr::filter(!!sym(diffexpressed_var) == dir)
        as.data.frame(sub_df[, id_cols, drop = FALSE])
      })
      names(sets) <- paste0("names_", comp, "_", ifelse(direction == "UP", "UP", "DN"))
    } else {
      sub_df <- data_results %>% dplyr::filter(!!sym(diffexpressed_var) %in% direction)
      sets <- list(as.data.frame(sub_df[, id_cols, drop = FALSE]))
      names(sets) <- paste0("names_", comp, "_UPDOWN")
    }
    sets
  })

  # Flatten comparisons -> flat list of gene sets
  gene_sets <- unlist(gene_sets, recursive = FALSE)

  if (flatten) {
    # Flatten one-column data frames -> plain named vectors
    gene_sets <- unlist(gene_sets, recursive = FALSE)
  }

  return(gene_sets)
}

