#' kegg analysis
#'
#' used to execute kegg enrichment
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
kegg_demo<-function(genesymbol,ntop=10,plot=T,species="org.Hs.eg.db"){
  library(clusterProfiler)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(DOSE)
  library(ggplot2)
  library(ggrepel)
  library(ggsci)

  sp<-ifelse(species!="org.Hs.eg.db",ifelse(species=="org.Mm.eg.db","mmu","rno"),"hsa")

  try({
    id_list <- mapIds(eval(parse(text = species)),genesymbol,"ENTREZID","SYMBOL")
    id_list <- na.omit(id_list)

    kegg <- enrichKEGG(id_list, keyType = 'kegg', pvalueCutoff = 0.05, pAdjustMethod = "BH",
                       minGSSize = 5, maxGSSize = 500, organism = sp, use_internal_data = FALSE)
    kegg<-DOSE::setReadable(kegg, OrgDb=species, keyType = 'ENTREZID')

    p<-ggplot(kegg[1:ntop,], aes(x=GeneRatio, y=Description,size=Count, color=pvalue)) + geom_point() +
      scale_colour_gradient(low="green",high="red") + labs(color=expression(padj),size="Gene number", x="GeneRatio",y="Pathway name",title="KEGG Pathway enrichment")
    return(list(kegg=kegg,plot=p))
  })
}
