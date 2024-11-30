#' Title
#'
#' @param exp the expression data contain one col indicate the gene name
#' @param symbol the colname of the gene name col
#' @param method how to merge the same gene count ,default as mean
#'
#' @return
#' @export
#'
#' @examples
gene_merge_demo<-function(exp,symbol,method=mean){
  exp = aggregate(exp, by=list(exp[,symbol]),method)
  exp = exp[!is.na(exp$Group.1) & exp$Group.1!="",]
  rownames(exp) = exp$Group.1
  exp<-exp[,!(names(exp) %in% c("Group.1",symbol))]
  return(exp)
}
