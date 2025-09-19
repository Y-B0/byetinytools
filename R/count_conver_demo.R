#' Title
#'
#' @param expr the expression matrix, rowname is the gene name and coulmn is the sample name
#' @param specie select specie from human, mouse and rat
#' @param genetype which kinds of gene name in expr select from hgnc_symbol, mgi_symbol and ensembl_gene_id.
#' @param from the count type of expr. select from count, fpkm, rpkm, tpm
#'
#' @return
#' @export
#'
#' @examples
count_conver_demo <- function(expr, specie = c("human", "mouse", "rat"), genetype = c("hgnc_symbol", "mgi_symbol", "ensembl_gene_id"), from) {
  library(biomaRt)

  # Match arguments
  specie <- match.arg(specie)
  genetype <- match.arg(genetype)

  # Choose appropriate biomart dataset
  if (specie == "human") {
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", "start_position", "end_position"),
                       filters = genetype, values = rownames(expr), mart = mart)
  } else if (specie == "mouse") {
    mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    gene_info <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position"),
                       filters = genetype, values = rownames(expr), mart = mart)
  } else if (specie == "rat") {
    mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
    gene_info <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "start_position", "end_position"),
                       filters = genetype, values = rownames(expr), mart = mart)
  }

  # Match genes
  gene_info <- gene_info[match(intersect(rownames(expr), gene_info[, genetype]), gene_info[, genetype]), ]
  expr <- expr[match(intersect(rownames(expr), gene_info[, genetype]), rownames(expr)), ]

  # Calculate gene length in kilobases (kb)
  gene_length <- (gene_info$end_position - gene_info$start_position + 1) / 1000

  # Define conversion functions
  rpkm <- function(values, gene_length) {
    total_counts <- sum(values)
    rpkm <- (values / gene_length) / (total_counts / 1e6)
    return(rpkm)
  }

  rpkm_to_fpkm <- function(values) {
    return(values / 2)
  }

  fpkm_to_rpkm <- function(values) {
    return(values * 2)
  }

  rpkm_to_tpm <- function(values) {
    total_rpkm <- sum(values)
    tpm <- (values / total_rpkm) * 1e6
    return(tpm)
  }

  count <- function(values, gene_lengths) {
    total_rpkm <- sum(values)
    total_counts_estimated <- total_rpkm
    counts <- values * gene_lengths * (total_counts_estimated / 1e6)
    return(counts)
  }

  # Convert based on 'from' argument
  if (from == "count") {
    rpkm_matrix <- apply(expr, 2, function(count) rpkm(count, gene_length))
    fpkm_matrix <- apply(rpkm_matrix, 2, rpkm_to_fpkm)
    tpm_matrix <- apply(rpkm_matrix, 2, rpkm_to_tpm)
    result <- list(count = expr, rpkm = rpkm_matrix, fpkm = fpkm_matrix, tpm = tpm_matrix)
  } else if (from == "rpkm") {
    fpkm_matrix <- apply(expr, 2, rpkm_to_fpkm)
    tpm_matrix <- apply(expr, 2, rpkm_to_tpm)
    count_matrix <- apply(expr, 2, function(rpkm) count(rpkm, gene_length))
    result <- list(rpkm = expr, fpkm = fpkm_matrix, tpm = tpm_matrix, count = round(count_matrix))
  } else if (from == "fpkm") {
    rpkm_matrix <- apply(expr, 2, fpkm_to_rpkm)
    tpm_matrix <- apply(rpkm_matrix, 2, rpkm_to_tpm)
    count_matrix <- apply(rpkm_matrix, 2, function(rpkm) count(rpkm, gene_length))
    result <- list(fpkm = expr, rpkm = rpkm_matrix, tpm = tpm_matrix, count = round(count_matrix))
  } else {
    stop("Invalid 'from' value. Must be one of 'count', 'rpkm', or 'fpkm'.")
  }

  result<-lapply(result, as.data.frame)

  return(result)
}
