source("src/h5adRead.R")
library(Seurat)


cap_h5ad_to_seurat <- function(h5ad_path){
    # create v5 assays
    options(Seurat.object.assay.version = "v5")

    adata <- read_h5ad(h5ad_path)
    
    if (is.null(adata$raw)) {
        # No raw layer
        main_assay <- CreateAssay5Object(counts = adata$X)
        main_assay <- AddMetaData(main_assay, adata$var)  # use raw as it is wider
    } else {
        # Raw layer exists
        main_assay <- CreateAssay5Object(counts = adata$raw$X, data = adata$X)
        main_assay <- AddMetaData(main_assay, adata$raw$var)  # use raw as it is wider
    }
    
    seurat_obj <- CreateSeuratObject(main_assay)
    seurat_obj <- AddMetaData(seurat_obj, adata$obs)
    Misc(seurat_obj, "uns") <- adata$uns

    for (emb in names(adata$obsm)) {
        matrix <- adata$obsm[[emb]]
        colnames(matrix) <-  paste0(emb, seq_len(ncol(adata$obsm[[emb]])))
        rownames(matrix) <- rownames(adata$obs)
        seurat_obj[[emb]] <- CreateDimReducObject(embeddings = matrix, key = paste0(emb, "_"), assay = "RNA")
    }
    return(seurat_obj)
}


h5ad2seuratObject <- function() {
    file <- "wks/mur.h5ad"
    srt <- cap_h5ad_to_seurat(file)
    return(srt)
}

h5ad2rds <- function(h5ad_path) {
    rds_path <- gsub(".h5ad", ".rds", h5ad_path)
    srt <- cap_h5ad_to_seurat(h5ad_path)
    saveRDS(srt, rds_path)
}
