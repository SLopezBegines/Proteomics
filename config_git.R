#' Configura Git en un proyecto de R ya existente
#'
#' Esta función automatiza la configuración de Git en un RProject existente.
#' Inicializa el repositorio Git (si no está), configura el usuario,
#' crea o actualiza el archivo .gitignore, vincula el repositorio con GitHub
#' y realiza un primer commit y push.
#'
#' @param nombre_usuario Nombre del usuario (como aparecerá en los commits)
#' @param email_usuario Email del usuario (asociado a GitHub)
#' @param url_repositorio URL del repositorio remoto en GitHub (ej. https://github.com/usuario/repositorio.git)
#' 
#' @return Mensajes informativos del progreso del proceso.
#' @export
#'
#' @examples
#' configurar_git_proyecto(
#'   nombre_usuario = "Santiago López Begines",
#'   email_usuario = "slopez@example.com",
#'   url_repositorio = "https://github.com/SLopezBegines/AD_ML_A4_study.git"
#' )
configurar_git_proyecto <- function(nombre_usuario, email_usuario, url_repositorio) {
  # Detectar archivo .Rproj en el directorio actual
  proyecto_rproj <- list.files(pattern = "\\.Rproj$")
  if (length(proyecto_rproj) == 0) stop("❌ No se encontró un archivo .Rproj en el directorio actual.")
  
  proyecto_path <- normalizePath(getwd())
  message(paste0("📂 Proyecto detectado en: ", proyecto_path))
  
  # Configuración global de Git
  system(paste('git config --global user.name "', nombre_usuario, '"', sep = ""))
  system(paste('git config --global user.email "', email_usuario, '"', sep = ""))
  
  # Inicializar Git si aún no existe
  if (!dir.exists(file.path(proyecto_path, ".git"))) {
    system("git init")
    message("✅ Repositorio Git inicializado.")
  } else {
    message("ℹ️ El repositorio Git ya estaba inicializado.")
  }
  
  # Crear o actualizar .gitignore
  gitignore_path <- file.path(proyecto_path, ".gitignore")
  lineas_a_ignorar <- c(
    ".Rhistory", ".RData", ".Rproj.user/", ".Ruserdata", ".DS_Store",
    "*.Rproj", "input_data/", "mains/output/","mains/output/figures/", "renv/", "renv.lock", ".Rprofile", "rawdata/"
  )
  
  if (!file.exists(gitignore_path)) {
    writeLines(lineas_a_ignorar, gitignore_path)
    message("📝 Archivo .gitignore creado.")
  } else {
    existente <- readLines(gitignore_path)
    nuevas <- setdiff(lineas_a_ignorar, existente)
    if (length(nuevas) > 0) {
      write(nuevas, gitignore_path, append = TRUE)
      message("📝 Archivo .gitignore actualizado.")
    } else {
      message("ℹ️ .gitignore ya contiene todas las reglas necesarias.")
    }
  }
  
  # Añadir remoto si no existe
  remotos <- system("git remote", intern = TRUE)
  if (!"origin" %in% remotos) {
    system(paste("git remote add origin", url_repositorio))
    message("🔗 Remoto origin añadido.")
  } else {
    message("ℹ️ Remoto origin ya está configurado.")
  }
  
  # Commit inicial si hay archivos no registrados
  system("git add .")
  estado <- system("git status --porcelain", intern = TRUE)
  if (length(estado) > 0) {
    system('git commit -m "Commit inicial del proyecto"')
    message("✅ Commit inicial realizado.")
  } else {
    message("ℹ️ No hay cambios para commitear.")
  }
  
  # Establecer rama y hacer push
  rama_actual <- system("git branch --show-current", intern = TRUE)
  if (length(rama_actual) == 0) {
    system("git branch -M main")
    rama_actual <- "main"
  }
  
  # Push al remoto
  system(paste("git push -u origin", rama_actual))
  message("📤 Proyecto sincronizado con GitHub.")
}

configurar_git_proyecto(
  nombre_usuario = "Santiago López Begines",
  email_usuario = "santiago.lopez.begines@gmail.com",
  url_repositorio = "https://github.com/SLopezBegines/Proteomics.git"
)

