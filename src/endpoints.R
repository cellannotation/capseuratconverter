source("src/context.R")
source("src/h5ad2rds.R")
library("httr")


h5ad2rds_endpoint <- function(h5ad_url, rds_push_url, push_to_bucket = TRUE) {
    # Input parameters: h5ad_url, rds_push_url, web_hook_success, web_hook_failure
        
    # Download h5ad file
    # TODO: add check for downloaded files already
    
    h5ad_name <- basename(h5ad_url)
    h5ad_path <- paste0(DATA_DIR, "/", h5ad_name)
    download.file(h5ad_url, h5ad_path)

    # Run conversion
    rds_path <- h5ad2rds(h5ad_path)

    # Push RDS to bucket
    # if (push_to_bucket) {
    #     httr$
    # }
}



