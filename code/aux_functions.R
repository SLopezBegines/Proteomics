
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

