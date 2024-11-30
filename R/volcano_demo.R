#' volcano plot
#'
#' @param x deg result
#' @param p.name the colname of pvalue
#' @param fc.name the colname of fold change
#' @param p.value the threshold of pvalue
#' @param fc.value the threshold of fold change
#' @param plot.name if want to save the volcano plot file please assign it
#' @param gene.repel if want to point some special gene，please input a vector of gene symbol
#'
#' @return deg result and volcano plot
#' @export
#'
#' @examples
volcano_demo<-function(x, p.name = "P.Value", fc.name = "logFC", p.value = 0.05, fc.value = 0.585,plot.name=NULL,gene.repel=NULL,width=7,height=4.5) {
  library(ggplot2)
  library(ggrepel)

  p <- ggplot() +
    theme_bw() +
    xlim(-10, 10) +
    geom_point(aes(x = x[, fc.name], y = -1 * log10(x[, p.name]), color = x$sig)) +
    theme(text = element_text(size = 20)) +
    labs(x = "log2(FoldChange)", y = paste("-log10(", p.name, ")", sep = "")) +
    scale_color_manual(name = "", values = c("#0072B5", "grey", "#BC3C28"))

  if (!is.null(gene.repel)) {
    data<- subset(x, rownames(x) %in% gene.repel)
    p<-p+geom_text_repel(
      aes(x = data[[fc.name]], y = -1 * log10(data[[p.name]]),label = rownames(data)),
      size = 3,
      box.padding = unit(0.5, "lines"),
      point.padding = unit(0.5, "lines"),
      max.overlaps = Inf,
      segment.color = "grey50"
    )
  }

  if (!is.null(plot.name)) {
    ggsave(plot.name, p, "pdf", width = width, height = height)
  }

  return(list(data = x, plot = p))
}
