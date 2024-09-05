#' deg analysis
#'
#' @param exp expression data (if merge==T,exp need a colum to contain gene symbol).
#' @param group group can be multi group, and first col is sample name; second col is condition; coef, and multi used to select multi group, keep in mind, control sample priority.
#' @param compared compared used to identify which group vs which group (usually <conventional> - <control>). also use for multi group condition.
#' @param normalize if the array expression data need to normalized between different sample
#' @param log2 if the expression data need to be convert to log2
#' @param merge if expression data contain same gene symol, you can use it to merge them
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
limma_demo<-function(exp,group,compared,normalize=F,log2=F,merge=F,rna.count=F,p.name = "P.Value", fc.name = "logFC", p.value = 0.05, fc.value = 0.585,file.name=NULL){
  library(limma)

  print(compared)

  if (merge==T) {
    exp = aggregate(exp, by=list(rownames(exp)),mean)
    exp = exp[!is.na(exp$Group.1) & exp$Group.1!="",]
    rownames(exp) = exp$Group.1
    exp<-exp[,-1]
  }

  if (rna.count==T) {
    design <- model.matrix(~0 + factor(group[, 2]))
    exp<- voom(exp,design,normalize="quantile")$E
  }

  if (normalize==T) {
    exp<-normalizeBetweenArrays(exp)
  }
  if (log2==T) {
    exp<-log2(exp+1)
  }

  design <- model.matrix(~0 + factor(group[, 2]))
  colnames(design) <- unique(group[, 2])
  rownames(design) <- group[, 1]
  exp <- exp[, rownames(design)]

  fit <- lmFit(exp, design)
  contrast.matrix <- makeContrasts(contrasts=compared, levels = colnames(design))

  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  x <- topTable(fit2, coef = compared, n = Inf, sort.by = "none")

  x$sig[(x[, p.name] > p.value | x[, p.name] == "NA") | (x[, fc.name] < fc.value) & x[, fc.name] > -fc.value] <- "Stable"
  x$sig[x[, p.name] <= p.value & x[, fc.name] >= fc.value] <- "Up"
  x$sig[x[, p.name] <= p.value & x[, fc.name] <= -fc.value] <- "Down"

  if (!is.null(file.name)) {
    write.table(data.frame(Symbol=rownames(output),x,exp),file = file.name,sep = "\t",quote = F,row.names = F,col.names = T)
  }
  return(cbind(x,exp))

}
