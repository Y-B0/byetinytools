#' survival analysis
#'
#' @param expr the data frame contain the expression of variates group and time, status. keep in mind, the colname of time and status must be 'time' and 'status'
#' @param group the colname that need to be compared group in expr
#'
#' @return survival plot
#' @export
#'
#' @examples
surv_demo <- function(expr,group){
  library("survival")
  library("survminer")

  model<-as.formula(paste('Surv(time, status)~', group))
  fit <- surv_fit(model, data = expr)

  p <- ggsurvplot(fit,data = expr,
                  pval = TRUE, conf.int = TRUE,
                  risk.table = TRUE, # Add risk table
                  risk.table.col = group, # Change risk table color by groups
                  surv.median.line = "hv", # Specify median survival
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c("#E7B800", "#2E9FDF"),font.x = 16,
                  font.y = 16,
                  font.title = 18,
                  font.legend = 14,
                  font.xtick = 14,
                  font.ytick = 14)

  return(p)
}
