#' Title
#'
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
chea3_demo<-function(genes){

  library(httr)
  library(jsonlite)

  url = "https://maayanlab.cloud/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = "myQuery", gene_set = genes)

  #POST to ChEA3 server
  response = POST(url = url, body = payload, encode = encode)
  json = content(response, "text")

  #results as list of R dataframes
  results = fromJSON(json)
  return(results)
}
