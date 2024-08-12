source("src/context.R")
source("src/h5ad2rds.R")
library("httr")
library(logger)


h5ad2rds_endpoint <- function(h5ad_url, rds_push_url, push_to_bucket = TRUE) {
    # Input parameters: h5ad_url, rds_push_url, web_hook_success, web_hook_failure
    # Download h5ad file
    # TODO: add check for downloaded files already
    log_debug("Start h5ad2rds_endpoint...")
    h5ad_name <- basename(h5ad_url)
    h5ad_path <- paste0(DATA_DIR, "/", h5ad_name)
    # download.file(h5ad_url, h5ad_path) TODO: uncomment this line
    log_debug(paste0("Downloaded h5ad file to: ", h5ad_path))

    # Run conversion
    rds_path <- h5ad2rds(h5ad_path)

    # Push RDS to bucket
    resp <- NULL
    if (push_to_bucket) {
        log_debug(paste0("Pushing RDS to bucket: ", rds_push_url))
        resp <- PUT(
            url=rds_push_url, 
            config=add_headers(`content-type` = "application/octet-stream"),
            body=upload_file(rds_path)
            )
        log_debug(paste0("Pushed RDS to bucket with response: ", resp))
    }
    return (resp)
}



