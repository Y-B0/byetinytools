#' go analysis
#'
#' @param genesymbol the enrichment gene list
#' @param ntop choose how many category should exhibit in barchart
#' @param plot if need to plot barchart
#' @param species choose species
#'
#' @return a list contain a datafram and a ggplot object
#' @export
#'
#' @examples
go_demo<-function(genesymbol,ntop=10,plot=T,plot.name=NULL,file.name=NULL,color="#4DBBD5FF",plot.width=7,plot.height=7,species=c("org.Hs.eg.db","org.Mm.eg.db"),ont="ALL",pvalueCutoff=0.05,qvalueCutoff=0.2,keytype="SYMBOL"){
  library(clusterProfiler)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(DOSE)
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
  library(dplyr)
  library(magrittr)
  library(enrichplot)
  try({
    go <- enrichGO(gene = genesymbol,
                   OrgDb = eval(parse(text=species)),
                   keyType = keytype,
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   readable = T,pvalueCutoff=pvalueCutoff,qvalueCutoff=qvalueCutoff)
    go.res <- data.frame(go)
    if (!is.null(file.name)) {
      write.table(go.res,,file = file.name,sep = "\t",quote = F,row.names = F,col.names = T)
    }


    if (plot==T) {
      if (ont!="ALL") {
        go.res<-go.res[go.res$ONTOLOGY%in%ont,]
      }
      go.df <- go.res %>%
        group_by(ONTOLOGY) %>%
        slice_head(n = 10) %>%
        ungroup()
      go.df<-na.omit(go.df)
      go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))

      go_bar<-ggplot() + geom_col(data = go.df, aes(x = -log10(pvalue), y = Description), fill = color, width = 0.5)+
        scale_x_continuous(expand = c(0,0)) + geom_text(data = go.df, aes(x = 0.1, y = Description, label = Description), size = 4.5, hjust = 0) +
        geom_text(data = go.df, aes(x = 0.1, y = Description, label = geneID),color = color, size = 4, hjust = 0, vjust = 2.4) +
        labs(x = expression(-Log[10]*"P"), y = '') + theme_classic() +
        theme(legend.position ='none', plot.title = element_text(size = 16), axis.title = element_text(size = 16),
              axis.text = element_text(size = 14), axis.ticks.y = element_blank(), axis.text.y = element_blank())


      if (!is.null(plot.name)) {
        ggsave(go_bar,filename = plot.name,width = plot.width,height = plot.height)
      }
    }
    return(list(go=go,plot=go_bar))
  })
}
