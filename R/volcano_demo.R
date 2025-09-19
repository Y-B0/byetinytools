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
volcano_demo<-function(x, p.name = "P.Value", fc.name = "logFC",curve=T, p.value = 0.05, fc.value = 0.585,plot.name=NULL,gene.repel=NULL,width=7,height=4.5) {
  library(ggplot2)
  library(ggrepel)


  p <- ggplot() +
    theme_bw() +
    xlim(-10, 10) + ylim(0, max(-log10(x[[p.name]]), na.rm = TRUE) * 1.2) +
    geom_point(aes(x = x[, fc.name], y = -1 * log10(x[, p.name]), color = x$sig),alpha=0.8) +
    theme(text = element_text(size = 20)) +
    labs(x = "log2(FoldChange)", y = paste("-log10(", p.name, ")", sep = "")) +
    scale_color_manual(name = "", values = c("#81b8df", "grey", "#fe817d"))

  if (curve==T) {
    build_curve_data <- function(max_x, log2FC_threshold, pval_threshold) {
      x_vals <- seq(0.0001, max_x, by = 0.0001)
      y_vals <- 1 / x_vals - log10(pval_threshold)
      # 构造左右对称的曲线数据（平移至 log2FC_threshold 位置）
      curve_df <- rbind(
        data.frame(x = x_vals + log2FC_threshold, y = y_vals),
        data.frame(x = -x_vals - log2FC_threshold, y = y_vals)
      )
      return(curve_df)
    }
    curve_data <- build_curve_data(10, fc.value, p.value)
    p <- ggplot() +
      theme_bw() +
      xlim(-10, 10) + ylim(0, max(-log10(x[[p.name]]), na.rm = TRUE) * 1.2) +
      geom_point(aes(x = x[, fc.name], y = -1 * log10(x[, p.name]), color = x$sig),alpha=0.8) +
      theme(text = element_text(size = 20)) +
      labs(x = "log2(FoldChange)", y = paste("-log10(", p.name, ")", sep = "")) +
      scale_color_manual(name = "", values = c("#81b8df", "grey", "#fe817d"))+
      geom_line(data = curve_data, aes(x = x, y = y),
                color = "black", linetype = "dashed", size = 0.7)
  }


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
