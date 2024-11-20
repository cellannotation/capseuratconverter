log_layout <- log4r::default_log_layout()
log_console_appender <- log4r::console_appender(layout = log_layout)

appenders <- list(log_console_appender)

if (exists("log_file_name")) {
  log_file_appender <- log4r::file_appender(
    log_file_name,
    append = TRUE,
    layout = log_layout
  )
  appenders <- list(log_console_appender, log_file_appender)
}

log <- log4r::logger(
  threshold = "DEBUG",
  appenders = appenders
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

#' Convert H5AD file to Seurat object
#'
#' @description This function creates a Seurat object from AnnData (h5ad) file.
#'
#' @details This function reads the h5ad file using 'rhdf5' package
#' and create named list with all fields existing in the file. After that,
#' Seurat5 object is created from the list.
#' The package was created for Cell Annotation Platform (CAP) project and expects that AnnData file was created with
#' python AnnData package of version 0.8+
#'
#' Resulting Seurat object will have the only assay named 'RNA'.
#'
#' @param h5ad_path The file path to the H5AD file.
#'
#' @return A Seurat5 object.
#'
#' @examples
#' seurat_obj <- h5ad_to_seurat("/path/to/h5ad_file.h5ad")
#' str(seurat_obj)
#'
#' @export
h5ad_to_seurat <- function(h5ad_path) {
  options(Seurat.object.assay.version = "v5")

  adata <- read_h5ad(h5ad_path)

  if (is.null(adata$raw)) {
    # No raw layer
    log_debug("No raw layer found in h5ad file! Use assay@layers$counts=adata.X")
    main_assay <- SeuratObject::CreateAssay5Object(counts = adata$X)
    log_debug("Add var as assay@meta.data")
    main_assay <- SeuratObject::AddMetaData(main_assay, adata$var) # use raw as it is wider
  } else {
    # Raw layer exists
    log_debug("Raw layer found in h5ad file! Use assay@layers$counts=adata.raw.X, assay@layers$data=adata.X")
    main_assay <- SeuratObject::CreateAssay5Object(counts = adata$raw$X, data = adata$X)
    log_debug("Add var as assay@meta.data")
    main_assay <- SeuratObject::AddMetaData(main_assay, adata$var) # use raw as it is wider
  }
  log_debug("Create Seurat ojbect from Assay5")
  seurat_obj <- SeuratObject::CreateSeuratObject(main_assay)
  log_debug("Add obs section as @meta.data")
  seurat_obj <- SeuratObject::AddMetaData(seurat_obj, adata$obs)
  log_debug("Add uns section as seurat@misc$uns")
  SeuratObject::Misc(seurat_obj, "uns") <- adata$uns

  log_debug("Add obsm section as seurat@reductions")
  for (emb in names(adata$obsm)) {
    matrix <- adata$obsm[[emb]]
    colnames(matrix) <- paste0(emb, seq_len(ncol(adata$obsm[[emb]])))
    rownames(matrix) <- rownames(adata$obs)
    seurat_obj[[emb]] <- SeuratObject::CreateDimReducObject(embeddings = matrix, key = paste0(emb, "_"), assay = "RNA")
    log_debug(paste0("Added embeddings: ", emb))
  }
  log_debug("Finish cap_h5ad_to_seurat!")
  return(seurat_obj)
}


#' Convert H5AD file to RDS format
#'
#' This function converts a H5AD file to RDS format.
#'
#' @param h5ad_path The path to the H5AD file.
#'
#' @details This function creates Seurat5 object via h5ad_to_seurat and saves it to RDS file named
#' the same way as input file but with .rds extension.
#'
#' @return String, the path to the generated RDS file.
#'
#' @examples
#' h5ad2rds("/path/to/h5ad_file.h5ad")
#'
#' @export
h5ad2rds <- function(h5ad_path) {
  srt <- h5ad_to_seurat(h5ad_path)
  rds_path <- gsub(".h5ad$", ".rds", h5ad_path)
  saveRDS(srt, rds_path)
  log_info(paste0("Convertion done! File saved to: ", rds_path))
  return(rds_path)
}


assert <- function(a, b, message) {
  if (a != b) {
    message <- paste("AssertionError: ", message)
    message <- paste0(message, "\n", a, " != ", b)
    stop(message)
  }
}


read_attr <- function(file, path, attr) {
  ats <- rhdf5::h5readAttributes(file, path)
  return(ats[attr])
}


read_encoding_version <- function(file, path) {
  return(read_attr(file, path, "encoding-version"))
}


read_encoding_type <- function(file, path) {
  return(read_attr(file, path, "encoding-type"))
}


read_sparse_matrix <- function(file, path, format, transpose = TRUE) {
  log_debug(paste0("Start read_sparse_matrix: ", path, " transpose=", transpose, "..."))
  assert(read_encoding_version(file, path), "0.1.0", "Sparse matrix must have encoding version 0.1.0")

  shape <- unlist(read_attr(file, path, "shape"))
  data <- rhdf5::h5read(file, paste0(path, "/data")) # non zero values
  indices <- rhdf5::h5read(file, paste0(path, "/indices")) # column indices for csr or row indices for csc
  indptr <- rhdf5::h5read(file, paste0(path, "/indptr")) # index pointers for csr or column pointers for csc

  # Ensure indices and indptr are integers
  indices <- as.integer(indices)
  indptr <- as.integer(indptr)

  if (!is(data, "numeric")) {
    log_debug("data is not numeric, convert it now...")
    data <- as.numeric(data)
  }

  if (format == "csr") {
    matrix <- Matrix::sparseMatrix(j = indices + 1, p = indptr, x = data, dims = shape)
  } else if (format == "csc") {
    matrix <- Matrix::sparseMatrix(i = indices + 1, p = indptr, x = data, dims = shape)
  } else {
    stop(paste0("Unknown format", format))
  }
  # NOTE: The X matrix in AnnData is [n_cells x n_genes]
  # But in Seurat it must be [n_genes x n_cells]
  # So, we have to convert it in most cases
  if (transpose) {
    log_debug("transpose matrix")
    matrix <- Matrix::t(matrix)
  }
  log_debug("Finish read_sparse_matrix!")
  return(matrix)
}


read_dense_matrix <- function(file, path, transpose) {
  log_debug(paste0("Start read_dense_matrix: ", path, " transpose=", transpose, "..."))
  assert(read_encoding_version(file, path), "0.2.0", "Dense array must have encoding version 0.2.0")
  x <- rhdf5::h5read(file, path)
  # It returs transposed matrix for some reason
  # Looks like it is a feature of R, so transpose if it is not actually needed
  if (transpose == FALSE) {
    log_debug("transpose matrix")
    x <- t(x)
  }
  log_debug("Finish read_dense_matrix!")
  return(x)
}


get_X <- function(file, path) {
  log_info(paste0("Start get_X at path ", path, " ..."))
  encoding_type <- read_encoding_type(file, path)

  if (encoding_type == "csr_matrix") {
    matrix <- read_sparse_matrix(file, path, "csr")
  } else if (encoding_type == "csc_matrix") {
    matrix <- read_sparse_matrix(file, path, "csc")
  } else if (encoding_type == "array") {
    matrix <- read_dense_matrix(file, path, TRUE)
  } else {
    stop(paste0("Unknown encoding type", encoding_type))
  }
  log_info("Finish get_X!")
  return(matrix)
}


read_df_col_array <- function(file, path) {
  log_debug(paste0("Start read_df_col_array: ", path, "..."))
  assert(read_encoding_version(file, path), "0.2.0", "The encoding version of <array> must be 0.2.0")
  col <- rhdf5::h5read(file, path)
  return(col)
}

read_df_col_str_array <- function(file, path) {
  log_debug(paste0("Start read_df_col_str_array: ", path, "..."))
  assert(read_encoding_version(file, path), "0.2.0", "The encoding version of <string-array> must be 0.2.0")
  col <- rhdf5::h5read(file, path)
  # strig arrays are stored as
  # variable length string data type
  # with a utf-8 encoding
  # I bet the next line is not needed but I will keep it for now
  col <- as.character(col)
  return(col)
}

read_df_col_cat <- function(file, path) {
  log_debug(paste0("Start read_df_col_cat: ", path, "..."))

  assert(read_encoding_version(file, path), "0.2.0", "The encoding version of <categorical> must be 0.2.0")

  code_key <- paste0(path, "/codes")
  cat_key <- paste0(path, "/categories")

  codes <- rhdf5::h5read(file, code_key)
  categories <- rhdf5::h5read(file, cat_key)

  # codes is a vector from 0 to n_categories - 1
  # need add -1 to handle empty values
  # TODO: check if recognized as NaN correctly
  categories <- c(NULL, categories)
  levels <- c(-1:(length(categories) - 2))

  col <- factor(codes, levels = levels, labels = categories)
  return(col)
}


read_df_col <- function(file, path) {
  # Assumes that column exists in the file
  # make check before the call!
  col_type <- read_encoding_type(file, path)

  if (col_type == "array") {
    col <- read_df_col_array(file, path)
  } else if (col_type == "categorical") {
    col <- read_df_col_cat(file, path)
  } else if (col_type == "string-array") {
    col <- read_df_col_str_array(file, path)
  } else {
    stop(paste0("Unknown column type", path))
  }
  log_debug(paste0("Finish read_df_col: ", path, "..."))
  return(col)
}


read_df <- function(file, path) {
  log_debug(paste0("Start read_df: ", path, "..."))
  # Code assumes that df exists, so check it before call the funciton!
  assert(read_encoding_type(file, path), "dataframe", "The encoding type of AnnData file must be dataframe")
  assert(read_encoding_version(file, path), "0.2.0", "The encoding version of AnnData file must be 0.2.0")

  index_col <- read_attr(file, path, "_index")
  col_order <- read_attr(file, path, "column-order")
  col_order <- unlist(col_order)

  index <- read_df_col(file, paste0(path, "/", index_col))
  df <- data.frame(row.names = index)

  for (col in col_order) {
    col_path <- paste0(path, "/", col)
    col_values <- read_df_col(file, col_path)
    df <- cbind(df, col = col_values)
  }

  colnames(df) <- col_order
  log_debug("Finish read_df!")
  return(df)
}


get_obs <- function(file) {
  log_info("Start get_obs ...")
  path <- "/obs"
  obs <- read_df(file, path)
  log_info("Finish get_obs!")
  return(obs)
}


get_var <- function(file, layer = NaN) {
  log_info("Start get_var ...")
  path <- "/var"
  if (layer == "raw") {
    path <- paste0("/raw", path)
  }
  var <- read_df(file, path)
  log_info("Finish get_var!")
  return(var)
}


adapt_naming <- function(name) {
  # Update the name of the embedding to be compatible with Seurat
  # Remove 'X_' prefix + replace all '_' to "."
  # X_pca -> pca; X_fake_obsm -> fake.obsm
  if (substr(name, 1, 2) == "X_") {
    name <- substr(name, 3, nchar(name))
  }
  name <- gsub("_", ".", name)
  return(name)
}


get_obsm <- function(file) {
  log_info("Start get_obsm ...")
  path <- "/obsm"
  assert(read_encoding_type(file, path), "dict", "The encoding type of AnnData obsm section must be <dict>")
  assert(read_encoding_version(file, path), "0.1.0", "The encoding version of AnnData <dict> must be 0.1.0")

  obsm_group <- rhdf5::H5Gopen(file, path)
  on.exit(rhdf5::H5Gclose(obsm_group))

  obsm_list <- rhdf5::h5ls(obsm_group, recursive = FALSE)
  obsm <- list()

  for (emb_name in obsm_list[["name"]]) {
    emb_path <- paste0(path, "/", emb_name)
    type <- read_encoding_type(file, emb_path)
    if (type == "array") {
      matrix <- read_dense_matrix(file, emb_path, FALSE)
    } else if (type == "dataframe") {
      # TODO: do we need to implement it?
      log_warn("DataFrame in obsm is not supported yet, skip it...")
    } else {
      stop(paste0("Unknown encoding type", emb_path))
    }
    adapted_name <- adapt_naming(emb_name)
    obsm[[adapted_name]] <- matrix
  }
  log_info("Finish get_obsm!")
  return(obsm)
}


get_raw_X_var <- function(file) {
  log_info("Start get_raw_X_var ...")
  path <- "/raw"

  raw_group <- NULL

  # TODO: perhaps it is better to explicitly get the list of entiities in the file and check if the raw exists
  tryCatch(
    {
      raw_group <- rhdf5::H5Gopen(file, path)
      on.exit(rhdf5::H5Gclose(raw_group))
    },
    error = function(e) {
      log_info("No raw data found")
    }
  )

  # If 'raw_group' is NULL, return early
  if (is.null(raw_group)) {
    return(raw_group)
  }

  raw.X <- get_X(file, paste0(path, "/X"))
  raw.var <- get_var(file, layer = "raw")
  res <- list(X = raw.X, var = raw.var)
  log_info("Finish get_raw_X_var!")
  return(res)
}


read_mapping <- function(file, path, transpose) {
  # Function assumes that the mapping exists in the file
  # Check it before the call!
  log_debug(paste0("Start read_mapping: ", path, "..."))

  assert(read_encoding_type(file, path), "dict", paste0("The encoding type of AnnData mapping section ", path, " must be <dict>"))
  assert(read_encoding_version(file, path), "0.1.0", paste0("The encoding version of AnnData <dict> ", path, " must be 0.1.0"))

  group <- rhdf5::H5Gopen(file, path)
  on.exit(rhdf5::H5Gclose(group))

  # recursively read all element in the group

  elements <- rhdf5::h5ls(group, recursive = FALSE)[["name"]]
  result <- list()

  for (element in elements) {
    element_path <- paste0(path, "/", element)
    element_type <- read_encoding_type(file, element_path)
    log_debug(paste0("Read element: ", element, " with type: ", element_type))
    if (element_type == "dict") {
      values <- read_mapping(file, element_path)
    } else if (element_type == "dataframe") {
      values <- read_df(file, element_path)
    } else if (element_type == "array") {
      values <- read_dense_matrix(file, element_path, transpose = transpose)
    } else if (element_type == "csr_matrix") {
      values <- read_sparse_matrix(file, element_path, "csr", transpose = transpose)
    } else if (element_type == "csc_matrix") {
      values <- read_sparse_matrix(file, element_path, "csc", transpose = transpose)
    } else if (element_type == "numeric-scalar" || element_type == "string") {
      values <- rhdf5::h5read(file, element_path)
    } else if (element_type == "string-array") {
      values <- read_df_col_str_array(file, element_path)
    } else {
      stop(paste0("Unsupported encoding type ", element_type))
    }
    new_name <- adapt_naming(element)
    result[[new_name]] <- values
  }
  log_debug("Finish read_mapping!")
  return(result)
}

get_layers <- function(file) {
  log_info("Start get_layers ...")
  path <- "/layers"
  layers_group <- NULL

  # TODO: perhaps it is better to explicitly get the list of entiities in the file and check if the layers exists
  tryCatch(
    {
      layers_group <- rhdf5::H5Gopen(file, path)
      on.exit(rhdf5::H5Gclose(layers_group))
    },
    error = function(e) {
      log_info("No layers data found")
    }
  )

  if (is.null(layers_group)) {
    return(NaN)
  }

  assert(read_encoding_type(file, path), "dict", "The encoding type of AnnData layers section must be <dict>")
  assert(read_encoding_version(file, path), "0.1.0", "The encoding version of AnnData <dict> must be 0.1.0")

  result <- read_mapping(file, path, transpose = TRUE)
  log_info("Finish get_layers!")
  return(result)
}


get_uns <- function(file) {
  log_info("Start get_uns ...")
  path <- "/uns"

  assert(read_encoding_type(file, path), "dict", "The encoding type of AnnData uns section must be <dict>")
  assert(read_encoding_version(file, path), "0.1.0", "The encoding version of AnnData <dict> must be 0.1.0")

  result <- read_mapping(file, path, transpose = FALSE)
  log_info("Finish get_uns!")
  return(result)
}


read_h5ad <- function(path) {
  log_info(paste0("Start read_h5ad: path=", path, " ..."))
  f <- rhdf5::H5Fopen(path)
  on.exit(rhdf5::H5Fclose(f))

  assert(read_encoding_version(f, "/"), "0.1.0", "The encoding version of AnnData file must be 0.1.0")
  assert(read_encoding_type(f, "/"), "anndata", "The encoding type of AnnData file must be anndata")

  x <- get_X(f, "/X")
  obs <- get_obs(f)
  var <- get_var(f)
  rownames(x) <- rownames(var)
  colnames(x) <- rownames(obs)
  log_debug("Rownames and colnames for X data updated")
  obsm <- get_obsm(f)
  raw <- get_raw_X_var(f)
  if (!is.null(raw)) {
    log_debug("Update rownames and colnames for raw data")
    rownames(raw$X) <- rownames(var)
    colnames(raw$X) <- rownames(obs)
  }
  layers <- get_layers(f)
  uns <- get_uns(f)

  anndata <- list(X = x, obs = obs, var = var, obsm = obsm, raw = raw, layers = layers, uns = uns)
  log_info("Finish read_h5ad!")
  return(anndata)
}
