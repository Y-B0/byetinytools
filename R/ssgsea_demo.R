#' ssgsea and gsva analysis
#'
#' @param exp a expression data with gene symbol as rowname and sample name as colname
#' @param pathway the gene sets list, if pathway is NULL it will default as kegg for human
#' @param method method can choose from gsva or ssgsea
#' @param species species are name in msigdbr (msigdbr)
#' @param category category are name in msigdbr
#' @param subcategory subcategory are name in msigdbr
#'
#' @return a matrix of gsea or gsve result
#' @export
#'
#' @examples
ssgsea_deom<-function(exp,pathway=NULL,method="gsva",species="Homo sapiens",category="C2",subcategory="KEGG"){
  library(msigdbr)
  library(GSVA)
  ### exp is a expression data with gene symbol as rowname and sample name as colname
  ### pathway is the gene sets list, if pathway is NULL it will default as kegg for human
  ### species, category and subcategory are name in msigdbr (msigdbr)
  ### method can choose from ssgsea or gsva

  if (is.null(pathway)) {
    pathway <- msigdbr(species = species, category = "C2", subcategory = "KEGG")
    pathway <- pathway %>% split(x = .$gene_symbol, f = .$gs_name)
    names(pathway) <- gsub("KEGG_","",names(pathway))%>%gsub("_"," ",.)
  }else {
    pathway=pathway
  }
  pathway <- lapply(pathway, unique)
  gs.exp <- gsva(as.matrix(exp),pathway, method = method)

}
