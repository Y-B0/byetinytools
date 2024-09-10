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
gsea_demo<-function(symbol,rank,geneset=NULL,go.ont="ALL",GO=TRUE,species.go="org.Hs.eg.db",KEGG=TRUE,species.kegg="hsa",n=10,pathway_ID=NULL){

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
      gseGO(
        geneList,
        ont = go.ont,
        OrgDb = species.go,
        keyType = "SYMBOL",
        pvalueCutoff = 1,
        pAdjustMethod = "BH")
    })
  }

  if (KEGG) {
    gse.KEGG <- try({
      gseKEGG(geneList_entrez,
              organism = species.kegg,
              keyType = "SYMBOL",
              pvalueCutoff = 1)
    })
  }

  if (!is.null(geneset)) {
    gse.geneset <- try({
      GSEA(geneList,
           TERM2GENE = geneset,
           pvalueCutoff = 1)
    })
  }

  if (!is.null(pathway_ID)) {
    p_go<-try({gseaplot2(gse.GO,geneSetID = pathway_ID,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_kegg<-try({gseaplot2(gse.KEGG,geneSetID = pathway_ID,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_geneset<-try({gseaplot2(gse.KEGG,geneSetID = pathway_ID,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
  } else {
    p_go<-try({gseaplot2(gse.GO,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_kegg<-try({gseaplot2(gse.KEGG,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
    p_geneset<-try({gseaplot2(gse.KEGG,geneSetID = 1:n,pvalue_table = F,base_size = 13,color = brewer.pal(10,'Paired'))})
  }

  try(return(list(gse.GO=gse.GO,gse.KEGG=gse.KEGG,gse.geneset=gse.geneset,p_go=p_go,p_kegg=p_kegg,p_geneset=p_geneset)))
}
