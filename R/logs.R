library(log4r)

read_environment_variable <- function(name, default_value) {
  v <- Sys.getenv(name)
  if (is.na(v) || v == "") {
    v <- default_value
  }
  return(v)
}

log_file_name <- read_environment_variable(
  "CAP_SEURAT_CONVERTER_LOG_FILE_NAME",
  paste(getwd(), "logfile.txt",
    sep = "/src/"
  )
)
log_threshold <- read_environment_variable(
  "CAP_SEURAT_CONVERTER_LOG_LEVEL",
  "DEBUG"
)

log_layout <- default_log_layout()
log_console_appender <- console_appender(layout = log_layout)
log_file_appender <- file_appender(log_file_name,
  append = TRUE,
  layout = log_layout
)
log <- log4r::logger(
  threshold = log_threshold,
  appenders = list(log_console_appender, log_file_appender)
)

log_info <- function(msg) {
  log4r::info(log, msg)
}

log_error <- function(msg) {
  log4r::error(log, msg)
}

log_debug <- function(msg) {
  log4r::debug(log, msg)
}

log_warn <- function(msg) {
  log4r::warn(log, msg)
}
