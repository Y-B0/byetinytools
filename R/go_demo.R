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
go_demo<-function(genesymbol,ntop=10,plot=T,species=c("org.Hs.eg.db","org.Mm.eg.db"),ont="ALL",pvalueCutoff=0.05,qvalueCutoff=0.2){
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
  try({
    id_list <- mapIds(eval(parse(text=species)),genesymbol,"ENTREZID","SYMBOL")
    id_list <- na.omit(id_list)

    go <- enrichGO(gene = id_list,
                   OrgDb = eval(parse(text=species)),
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   readable = T,pvalueCutoff=pvalueCutoff,qvalueCutoff=qvalueCutoff
    )
    go.res <- data.frame(go)
    #write.csv(go.res,"Table_GO_result.csv",quote = F)


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
      go_bar <- ggplot(data = go.df,
                       aes(x = Description, y = Count,fill = ONTOLOGY))+
        geom_bar(stat = "identity",width = 0.9)+
        coord_flip()+theme_bw()+
        scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+
        labs(x = "GO terms",y = "GeneNumber",title = "Barplot of Enriched GO Terms")+
        theme(text=element_text(size=18), axis.text.y = element_text(size = 10),
              plot.title = element_text(size = 14,hjust = 0.5,face = "bold"),
              plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))+scale_fill_aaas()
      #ggsave(go_bar,filename = "GO_Barplot.pdf",width = 9,height = 7)
    }
    return(list(go=go,plot=go_bar))
  })
}
