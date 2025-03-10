#' deg analysis
#'
#' @param exp expression data (if merge==T,exp need a colum to contain gene symbol).
#' @param group group is the vector which indicate the sample's group
#' @param compared compared used to identify which group vs which group (usually <conventional> - <control>). also use for multi group condition.
#' @param normalize if the array expression data need to normalized between different sample
#' @param log2 if the expression data need to be convert to log2
#' @param merge if expression data contain same gene symol, you can use it to merge them
#' @param symbol accordance with "merge==T", specific the colum name of symbol
#' @param rna.count used to execute rna-seq
#' @param p.name the colname of pvalue
#' @param fc.name the colname of fold change
#' @param p.value the threshold of pvalue
#' @param fc.value the threshold of fold change
#' @param file.name if need to write a txt file to save the deg result, you can assign a file name
#'
#' @return deg result merged with normalized expression data
#' @export
#'
#' @examples
limma_demo<-function(exp,group,compared,normalize=F,log2=F,merge=F,symbol=NULL,rna.count=F,p.name = "P.Value", fc.name = "logFC", p.value = 0.05, fc.value = 0.585,add_expr=T,file.name=NULL){
  library(limma)
  library(edgeR)
  library(byetinytools)
  library(DESeq2)

  print(compared)

  if (merge==T) {
    exp<-gene_merge_demo(exp,symbol)
  }

  if (rna.count==T) {
    design <- model.matrix(~0 + factor(group))
    exp <- DGEList(counts = exp, group = factor(group))
    exp <- calcNormFactors(exp)
    exp<- voom(exp,design)$E
  }

  if (normalize==T) {
    exp<-normalizeBetweenArrays(exp)
  }
  if (log2==T) {
    exp<-log2(exp+1)
  }

  design <- model.matrix(~0 + factor(group))
  colnames(design) <- levels(factor(group))
  rownames(design) <- colnames(exp)
  #exp <- exp[, rownames(design)]

  fit <- lmFit(exp, design)
  contrast.matrix <- makeContrasts(contrasts=compared, levels = colnames(design))

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = compared, n = Inf, sort.by = "none")

  x$sig[(x[, p.name] > p.value | x[, p.name] == "NA") | (x[, fc.name] < fc.value) & x[, fc.name] > -fc.value] <- "Stable"
  x$sig[x[, p.name] <= p.value & x[, fc.name] >= fc.value] <- "Up"
  x$sig[x[, p.name] <= p.value & x[, fc.name] <= -fc.value] <- "Down"

  if (add_expr) {
    x<-cbind(x,exp)
  }else{
    x<-x
  }

  if (!is.null(file.name)) {
    write.table(data.frame(Symbol=rownames(x),x,exp),file = file.name,sep = "\t",quote = F,row.names = F,col.names = T)
  }
  return(x)

}
