source("src/h5adRead.R")
library(Seurat)


cap_h5ad_to_seurat <- function(h5ad_path){
    log_debug("Start cap_h5ad_to_seurat...")
    # create v5 assays
    options(Seurat.object.assay.version = "v5")

    adata <- read_h5ad(h5ad_path)
    
    if (is.null(adata$raw)) {
        # No raw layer
        log_debug("No raw layer found in h5ad file! Use assay@layers$counts=adata.X")
        main_assay <- CreateAssay5Object(counts = adata$X)
        log_debug("Add var as assay@meta.data")
        main_assay <- AddMetaData(main_assay, adata$var)  # use raw as it is wider
    } else {
        # Raw layer exists
        log_debug("Raw layer found in h5ad file! Use assay@layers$counts=adata.raw.X, assay@layers$data=adata.X")
        main_assay <- CreateAssay5Object(counts = adata$raw$X, data = adata$X)
        log_debug("Add raw.var as assay@meta.data")
        main_assay <- AddMetaData(main_assay, adata$raw$var)  # use raw as it is wider
    }
    log_debug("Create Seurat ojbect from Assay5")
    seurat_obj <- CreateSeuratObject(main_assay)
    log_debug("Add obs section as @meta.data")
    seurat_obj <- AddMetaData(seurat_obj, adata$obs)
    log_debug("Add uns section as seurat@misc$uns")
    Misc(seurat_obj, "uns") <- adata$uns

    log_debug("Add obsm section as seurat@reductions")
    for (emb in names(adata$obsm)) {
        matrix <- adata$obsm[[emb]]
        colnames(matrix) <-  paste0(emb, seq_len(ncol(adata$obsm[[emb]])))
        rownames(matrix) <- rownames(adata$obs)
        seurat_obj[[emb]] <- CreateDimReducObject(embeddings = matrix, key = paste0(emb, "_"), assay = "RNA")
        log_debug(paste0("Added embeddings: ", emb))
    }
    log_debug("Finish cap_h5ad_to_seurat!")
    return(seurat_obj)
}


h5ad2seuratObject <- function() {
    file <- "wks/mur.h5ad"
    srt <- cap_h5ad_to_seurat(file)
    return(srt)
}


h5ad2rds <- function(h5ad_path) {
    log_info(paste0("Start converting h5ad: ", h5ad_path, " to Seurat v5 RDS..."))
    rds_path <- gsub(".h5ad", ".rds", h5ad_path)
    srt <- cap_h5ad_to_seurat(h5ad_path)
    log_info(paste0("Convertion done! File saved to: ", rds_path))
    saveRDS(srt, rds_path)
}
