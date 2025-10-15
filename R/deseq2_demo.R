#' Title
#'
#' @param exp expression data (if merge==T,exp need a colum to contain gene symbol).
#' @param group group is the vector which indicate the sample's group
#' @param compared compared used to identify which group vs which group (usually <conventional> - <control>). also use for multi group condition.
#' @param merge if expression data contain same gene symol, you can use it to merge them
#' @param symbol accordance with "merge==T", specific the colum name of symbol
#' @param p.name the colname of pvalue
#' @param fc.name the colname of fold change
#' @param p.value the threshold of pvalue
#' @param fc.value the threshold of fold change
#' @param file.name if need to write a tdata:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCCxt file to save the deg result, you can assign a file name
#'
#' @return
#' @export
#'
#' @examples
deseq2_demo<-function(exp, group, compared,merge=F, p.name = "pvalue", symbol=NULL, fc.name = "log2FoldChange", p.value = 0.05, fc.value = 0.585,add_expr=T,file.name=NULL){
  library(DESeq2)
  library(stringr)
  library(magrittr)
  if (merge==T) {
    exp = aggregate(exp, by=list(exp[,symbol]),mean)
    exp = exp[!is.na(exp$Group.1) & exp$Group.1!="",]
    rownames(exp) = exp$Group.1
    exp<-exp[,!(names(exp) %in% c("Group.1",symbol))]
  }

  compared<-as.vector(str_split(compared,"-",simplify = T))%>%str_trim()

  dds <- DESeqDataSetFromMatrix(exp, data.frame(group), design= ~ group )
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  dds <- DESeq(dds)
  x <- results(dds,contrast = c("group",compared[[1]],compared[[2]]))%>%as.data.frame()
  count<-counts(dds,normalized=T)

  x$sig[(x[, p.name] > p.value | x[, p.name] == "NA") | (x[, fc.name] < fc.value) & x[, fc.name] > -fc.value] <- "Stable"
  x$sig[x[, p.name] <= p.value & x[, fc.name] >= fc.value] <- "Up"
  x$sig[x[, p.name] <= p.value & x[, fc.name] <= -fc.value] <- "Down"

  if (add_expr) {
    x<-cbind(x,count)
  }else{
    x<-x
  }

  if (!is.null(file.name)) {
    write.table(data.frame(Symbol=rownames(x),x,exp),file = file.name,sep = "\t",quote = F,row.names = F,col.names = T)
  }

  return(x)
}

