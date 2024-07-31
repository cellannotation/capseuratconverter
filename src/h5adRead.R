# The code is prepared for AnnData 0.8+ version
# https://anndata.readthedocs.io/en/latest/fileformat-prose.html

library(rhdf5)
library(Matrix)


assert <- function(a, b, message){
    if (a != b) {
        message <- paste("AssertionError: ", message)
        message <- paste0(message, "\n", a, " != ", b)
        stop(message)
    }
}


read_attr <- function(file, path, attr){
    ats <- h5readAttributes(file, path)
    return (ats[attr])
}


read_encoding_version <- function(file, path){
    return (read_attr(file, path, 'encoding-version'))
}


read_encoding_type <- function(file, path){
    return (read_attr(file, path, 'encoding-type'))
}


read_sparse_matrix <- function(file, path, format, transpose=TRUE){
    assert(read_encoding_version(file, path), "0.1.0", "Sparse matrix must have encoding version 0.1.0")
    
    shape <- unlist(read_attr(file, path, 'shape'))
    data <- h5read(file, paste0(path, "/data"))  # non zero values
    indices <- h5read(file, paste0(path, "/indices"))  # column indices for csr or row indices for csc
    indptr <- h5read(file, paste0(path, "/indptr"))  # index pointers for csr or column pointers for csc

    if (format == "csr") {
        matrix <- sparseMatrix(j = indices + 1, p = indptr, x = data, dims = shape)
    } else if (format == "csc") {
        matrix <- sparseMatrix(i = indices + 1, p = indptr , x = data, dims = shape)
    } else {
        stop(paste0("Unknown format", format))
    }
    # NOTE: The X matrix in AnnData is [n_cells x n_genes]
    # But in Seurat it must be [n_genes x n_cells]
    # So, we have to convert it in most cases
    if (transpose) { matrix <- t(matrix) }
    return (matrix)
}


read_dense_matrix <- function(file, path, transpose){
    assert(read_encoding_version(file, path), "0.2.0", "Dense array must have encoding version 0.2.0")
    x <- h5read(file, path)
    # It returs transposed matrix for some reason
    # Looks like it is a feature of R, so transpose if it is not actually needed
    if (transpose == FALSE) { x <- t(x)}
    return (x)
}


#' Must get an X object (either a h5 Dataset or h5 Group)
#' @param file is a path to h5ad file
#' @param path is a path to the X object in h5ad file structure
#' @return a matrix object
get_X <- function(file, path){    
    encoding_type <- read_encoding_type(file, path)
    
    if (encoding_type == 'csr_matrix') {
        matrix <- read_sparse_matrix(file, path, "csr")
    } else if (encoding_type == 'csc_matrix') {
        matrix <- read_sparse_matrix(file, path, "csc")
    } else if (encoding_type == 'array') {
        matrix <- read_dense_matrix(file, path, TRUE)
    } else {
        stop(paste0("Unknown encoding type", encoding_type))
    }
    return (matrix)
}


read_df_col_array <- function(file, path){
    assert(read_encoding_version(file, path), '0.2.0', "The encoding version of <array> must be 0.2.0")
    col <- h5read(file, path)
    return (col)
}

read_df_col_str_array <- function(file, path){
    assert(read_encoding_version(file, path), '0.2.0', "The encoding version of <string-array> must be 0.2.0")
    col <- h5read(file, path)
    # strig arrays are stored as
    # variable length string data type
    # with a utf-8 encoding 
    # I bet the next line is not needed but I will keep it for now
    col <- as.character(col)  
    return (col)
}

read_df_col_cat <- function(file, path){
    assert(read_encoding_version(file, path), '0.2.0', "The encoding version of <categorical> must be 0.2.0")

    code_key <- paste0(path, "/codes")
    cat_key <- paste0(path, "/categories")

    codes <- h5read(file, code_key)  
    categories <- h5read(file, cat_key)

    # codes is a vector from 0 to n_categories - 1
    # need add -1 to handle empty values
    # TODO: check if recognized as NaN correctly
    categories <- c(NaN, categories)
    levels <- c(-1:(length(categories)-2))

    col <- factor(codes, levels = levels, labels = categories)
    return (col)
}


read_df_col <- function(file, path){
    # Assumes that column exists in the file
    # make check before the call!
    col_type <- read_encoding_type(file, path)

    if (col_type == 'array'){
        print("Get numberic array")
        col <- read_df_col_array(file, path)
    } else if (col_type == 'categorical'){
        print("Get categorical")
        col <- read_df_col_cat(file, path)
    } else if (col_type == "string-array"){
        print("Get string array")
        col <- read_df_col_str_array(file, path)
    } else {
        stop(paste0("Unknown column type", path))
    }
    return (col)
}


read_df <- function(file, path){
    # Code assumes that df exists, so check it before call the funciton!
    assert(read_encoding_type(file, path), 'dataframe', "The encoding type of AnnData file must be dataframe")
    assert(read_encoding_version(file, path), '0.2.0', "The encoding version of AnnData file must be 0.2.0")

    index_col <- read_attr(file, path, '_index')
    print(index_col)
    col_order <- read_attr(file, path, 'column-order')
    col_order <- unlist(col_order)

    index <- read_df_col(file, paste0(path, "/", index_col))
    df <- data.frame(row.names = index)
        
    for (col in col_order){
        print(paste0("try to handle column: ", col))
        col_path <- paste0(path, "/", col)
        col_values <- read_df_col(file, col_path)
        df <- cbind(df, col = col_values)
    }

    colnames(df) <- col_order

    return (df)
}


get_obs <- function(file){
    path = "/obs"
    return (read_df(file, path))
}


get_var <- function(file, layer = NaN){
    path = "/var"
    if (layer == "raw") {path <- paste0("/raw", path)}
    return (read_df(file, path))
}


adapt_naming <- function(name){
    # Update the name of the embedding to be compatible with Seurat
    # Remove 'X_' prefix + replace all '_' to "."
    # X_pca -> pca; X_fake_obsm -> fake.obsm
    if (substr(name, 1, 2) == "X_") {name <- substr(name, 3, nchar(name))}
    name <- gsub("_", ".", name)
    return (name)
}


get_obsm <- function(file){
    path <- "/obsm"
    assert(read_encoding_type(file, path), "dict", "The encoding type of AnnData obsm section must be <dict>")
    assert(read_encoding_version(file, path), "0.1.0", "The encoding version of AnnData <dict> must be 0.1.0")

    obsm_group = H5Gopen(file, path)
    on.exit(H5Gclose(obsm_group))

    obsm_list <- h5ls(obsm_group, recursive = FALSE)
    obsm <- list()

    for (emb_name in obsm_list[['name']]){
        print(emb_name)
        
        emb_path = paste0(path, "/", emb_name)
        type = read_encoding_type(file, emb_path)
        if (type == 'array') {
            matrix <- read_dense_matrix(file, emb_path, FALSE)
        } else if (type == "dataframe") {
            # TODO: do we need to implement it?
            print("DataFrame in obsm is not supported yet")
        }else {
            stop(paste0("Unknown encoding type", emb_path))
        }
        adapted_name <- adapt_naming(emb_name)
        obsm[[adapted_name]] <- matrix
    }
    return (obsm)
}


get_raw_X_var <- function(file){
    path <- "/raw"

    raw_group <- NULL

    # TODO: perhaps it is better to explicitly get the list of entiities in the file and check if the raw exists
    tryCatch({
        raw_group = H5Gopen(file, path)
        on.exit(H5Gclose(raw_group))
    }, error = function(e){
        print("No raw data found")
    })

    # If 'raw_group' is NULL, return early
    if (is.null(raw_group)) {
        return(NaN)
    }

    raw.X <- get_X(file, paste0(path, "/X"))
    raw.var <- get_var(file, layer = "raw")
    res <- list(X = raw.X, var = raw.var)
    return (res)
}


read_mapping <- function(file, path, transpose) {
    # Function assumes that the mapping exists in the file
    # Check it before the call!

    assert(read_encoding_type(file, path), "dict", paste0("The encoding type of AnnData mapping section ", path ," must be <dict>"))
    assert(read_encoding_version(file, path), "0.1.0", paste0("The encoding version of AnnData <dict> ", path ," must be 0.1.0"))

    group <- H5Gopen(file, path)
    on.exit(H5Gclose(group))

    # recursively read all element in the group

    elements <- h5ls(group, recursive = FALSE)[['name']]
    result <- list()

    for (element in elements){
        element_path <- paste0(path, "/", element)
        element_type <- read_encoding_type(file, element_path)
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
            values <- h5read(file, element_path)
        } else {
            stop(paste0("Unsupported encoding type", element_type))
        }

        
        new_name <- adapt_naming(element)
        result[[ new_name ]] <- values
    }

    return (result)
}

get_layers <- function(file) {
    path <- "/layers"     
    layers_group <- NULL

    # TODO: perhaps it is better to explicitly get the list of entiities in the file and check if the layers exists
    tryCatch({
        layers_group = H5Gopen(file, path)
        on.exit(H5Gclose(layers_group))
    }, error = function(e){
        print("No layers data found")
    })

    if (is.null(layers_group)) {
        return(NaN)
    }

    assert(read_encoding_type(file, path), "dict", "The encoding type of AnnData layers section must be <dict>")
    assert(read_encoding_version(file, path), "0.1.0", "The encoding version of AnnData <dict> must be 0.1.0")

    result <- read_mapping(file, path, transpose = TRUE)
    return (result)
}


get_uns <- function(file) {
    path <- "/uns"

    assert(read_encoding_type(file, path), "dict", "The encoding type of AnnData uns section must be <dict>")
    assert(read_encoding_version(file, path), "0.1.0", "The encoding version of AnnData <dict> must be 0.1.0")

    result <- read_mapping(file, path, transpose = FALSE)
    return (result)
}


read_h5ad <- function(path) {
    f <- H5Fopen(path)
    on.exit(H5Fclose(f))
    
    assert(read_encoding_version(f, "/"), '0.1.0', "The encoding version of AnnData file must be 0.1.0")
    assert(read_encoding_type(f, "/"), 'anndata', "The encoding type of AnnData file must be anndata")
    
    x <- get_X(f, "/X")
    obs <- get_obs(f)
    var <- get_var(f)
    rownames(x) <- rownames(var)
    colnames(x) <- rownames(obs)
    obsm <- get_obsm(f)
    raw <- get_raw_X_var(f)
    layers <- get_layers(f)
    uns <- get_uns(f)

    anndata <- list(X = x, obs = obs, var = var, obsm = obsm, raw = raw, layers = layers, uns = uns)
    return (anndata)
}


