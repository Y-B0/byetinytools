#' Title
#'
#' @param symbol gene symbol vector for the enrichment, must matched with rank
#' @param rank rank vector for the symbol, only matched with symbol, can be logfc or some sort index for symbol
#' @param geneset can specific some special gene sets
#' @param go.ont go.not can be "ALL", "BP", "MF" and "CC"
#' @param GO whether need go gene set
#' @param species.go specify the specie
#' @param KEGG whether need kegg gene set
#' @param species.kegg specify the specie
#' @param n how many pathway need to be plot
#' @param pathway_ID specify some interest geneset
#'
#' @return
#' @export
#'
#' @examples
gsea_demo<-function(symbol,rank,geneset=NULL,go.ont="ALL",GO=TRUE,species=c("org.Hs.eg.db","org.Mm.eg.db"),KEGG=TRUE,n=10,pathway_ID=NULL,pvalueCutoff = 0.05){

  library(ggplot2)
  library(stringr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(Hmisc)
  library(enrichplot)
  library(RColorBrewer)
  library(org.Mm.eg.db)

  names(rank)<-symbol
  geneList<-sort(rank,decreasing = T)

  if (GO) {
    gse.GO <- try({
      dat<-gseGO(
        geneList,
        ont = go.ont,
        OrgDb = species,
        keyType = "SYMBOL",
        pvalueCutoff = pvalueCutoff,
        pAdjustMethod = "BH")
      dat@result<-dat@result[order(dat@result$NES),]
      dat
    })
  }

  if (KEGG) {
    gse.KEGG <- try({
      sp<-ifelse(species!="org.Hs.eg.db",ifelse(species=="org.Mm.eg.db","mmu","rno"),"hsa")
      id_list <- mapIds(eval(parse(text = species)),names(geneList),"ENTREZID","SYMBOL")%>%unlist()
      names(geneList)<-as.character(id_list)
      geneList <- na.omit(geneList)
      dat<-gseKEGG(geneList,
              organism = sp,pAdjustMethod = "BH",pvalueCutoff = pvalueCutoff)
      dat@result<-dat@result[order(dat@result$NES),]
      dat
    })
  }

  if (!is.null(geneset)) {
    gse.geneset <- try({
      dat<-GSEA(geneList,
           TERM2GENE = geneset,
           pvalueCutoff = 1)
      dat@result<-dat@result[order(dat@result$NES),]
      dat
    })
  } else {
    gse.geneset<-NULL
  }

  if (!is.null(pathway_ID)) {
    p_go<-try({gseaplot2(gse.GO,geneSetID = pathway_ID,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_kegg<-try({gseaplot2(gse.KEGG,geneSetID = pathway_ID,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_geneset<-try({gseaplot2(gse.geneset,geneSetID = pathway_ID,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
  } else {
    p_go<-try({gseaplot2(gse.GO,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_kegg<-try({gseaplot2(gse.KEGG,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_geneset<-try({gseaplot2(gse.geneset,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
  }

  try(return(list(gse.GO=gse.GO,gse.KEGG=gse.KEGG,gse.geneset=gse.geneset,p_go=p_go,p_kegg=p_kegg,p_geneset=p_geneset)))

}
