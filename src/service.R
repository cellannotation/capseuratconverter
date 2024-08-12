# source("src/h5ad2rds.R")
source("src/endpoints.R")
library(RestRserve)
library(httr)


app = Application$new()


app$add_get(
  path = "/echo", 
  FUN = function(.req, .res) {
    resp <- paste0("Got message: ", .req$parameters_query$message)
    .res$set_content_type("application/json")
    .res$set_body(list(message = resp))
  })


app$add_post(
  path = "/addone", 
  FUN = function(.req, .res) {
    result = list(x = .req$body$x + 1L)
    .res$set_content_type("application/json")
    .res$set_body(result)
  })

app$add_get(
  path = "/h5ad2rds",
  FUN = function(.req, .res) {
    h5ad_url <- .req$parameters_query$h5ad_url
    rds_push_url <- .req$parameters_query$rds_push_url
    web_hook_failure <- .req$parameters_query$web_hook_failure
    web_hook_success <- .req$parameters_query$web_hook_success

    tryCatch({
      response <- h5ad2rds_endpoint(h5ad_url, rds_push_url, push_to_bucket = TRUE)

      if (is.null(response) || response$status_code != 200) {
        stop("Failed to push RDS to bucket with status code: ", response$status_code)
        }

      # POST(
      #   url=web_hook_success, 
      #   body=list(message = "RDS was pushed to bucket successfully"), 
      #   encode = "json"
      # )
      
      }, error = function(e) {
        # POST(url=web_hook_failure, body=list(message = e$message), encode = "json")
        log_error(paste0("Failed to convert h5ad to RDS: ", e$message))
      }
    )
  }
)


backend = BackendRserve$new()
backend$start(app, http_port = 8080)
