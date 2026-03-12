scrna_step<-function(){

  #load enviroment
  env_load_s1<-function(){
    library(Seurat)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(clustree)
    library(patchwork)
    library(future)
    library(DoubletFinder)
    library(scMayoMap)
    library(starTracer)
    library(dplyr)
    library(future.apply)
    library(recall)
    library(parallel)
  }

  #file read in
  file_read_s2<-function(path,specise="human",min.cells=3,min.features=200,Rfile=NULL,h5=F,...){

    ##default the input path contain multi sample folder, each sample folder have named accordant with sample info file, every sample folder contain the three basic files,

    if (h5==T) {
      dir_name <<- list.files(path,pattern = "*.h5*",recursive = T,full.names = F)
    }else{
      dir_name <<- list.dirs(path,recursive = F)
    }

    scRNAlist <- list()
    for (i in dir_name) {
      if (h5==T) {
        counts <- Read10X_h5(filename = paste(path,i,sep = "/"),...)
      } else{
        counts <- Read10X(data.dir = i,...)
      }
      scRNAlist[[i]] <- CreateSeuratObject(counts = counts,
                                           project = str_split(i,"\\\\|/",simplify = T)%>%.[,ncol(.)],
                                           min.cells = min.cells,
                                           min.features = min.features)
    }
    names(scRNAlist)<-str_split(dir_name,"\\\\|/",simplify = T)%>%.[,ncol(.)]

    ## feature culculate
    if(specise=="human"){
      pattern<-"MT-"
      HB_genes_raw <- c("HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ")
    } else if(specise=="mouse"){
      pattern<-"mt-"
      HB_genes_raw <- c("Hba1|Hba2|Hbb|Hbd|Hbe1|Hbg1|Hbg2|Hbm|Hbq1|Hbz")
    }
    for (i in 1:length(scRNAlist)) {
      sc <- scRNAlist[[i]]
      sc[['MT_percent']] <- PercentageFeatureSet(sc,pattern = pattern)
      HB_m <- match(grep(HB_genes_raw,rownames(sc@assays$RNA),ignore.case = T,value = T),rownames(sc@assays$RNA))
      HB_genes <- rownames(sc@assays$RNA)[HB_m]
      HB_genes <- HB_genes[!is.na(HB_genes)]
      sc[['HB_percent']] <- PercentageFeatureSet(sc,features = HB_genes)
      scRNAlist[[i]] <- sc
      rm(sc)
    }
    print(paste("Totally ",length(scRNAlist)," sample loaded",sep = ""))

    raw_data<-scRNAlist
    if (!is.null(Rfile)) {
      saveRDS(raw_data,file = Rfile)
    }
    return(raw_data)
  }

  #feature exhibition
  feature_show_s3<-function(scRNAlist,features=c('nFeature_RNA','nCount_RNA','MT_percent','HB_percent'),merge=F,stack = F,...){
    violin<-list()
    for (i in 1:length(scRNAlist)){
      violin[[i]] <- VlnPlot(scRNAlist[[i]],
                             features = features,
                             pt.size = 0.01,
                             ncol = length(features),stack = stack,...)
    }
    if (!merge) {
      return(violin)
    } else {
      p<-purrr::reduce(violin, `/`)
      return(p)
    }

  }

  #feature filer
  cell_filter_s4<-function(scRNAlist,specise="human",nFeature_RNA_min=200,nFeature_RNA_max=4000,mt_percent_max=10,HB_percent_max=5,nCount_RNA_min=1000,
                           nCount_RNA_max=20000,specfic_gene=T,other_specific=NULL,innergene_rate=0,Rfile=NULL){

    filter_data <- lapply(scRNAlist, function(x) {
      subset(x,
             nFeature_RNA > nFeature_RNA_min & nFeature_RNA < nFeature_RNA_max &
               MT_percent < mt_percent_max & HB_percent < HB_percent_max &
               nCount_RNA > nCount_RNA_min & nCount_RNA < nCount_RNA_max)
    })

    if (specfic_gene==T) {

      if(specise=="human"){
        MT_gene <- "MT-"
        HB_gene <- c("HBA1|HBA2|HBB|HBD|HBE1|HBG1|HBG2|HBM|HBQ1|HBZ")
      } else if(specise=="mouse"){
        MT_gene <- "mt-"
        HB_gene <- c("Hba1|Hba2|Hbb|Hbd|Hbe1|Hbg1|Hbg2|Hbm|Hbq1|Hbz")
      }

      filter_data <- lapply(filter_data, function(x) {
        x<-subset(x, features = setdiff(rownames(x), c(other_specific,HB_gene,grep("^MT-",rownames(x),value = T))))
        expr_genes <- rowSums(as.matrix(x[["RNA"]]$counts) > 0) / ncol(x) > innergene_rate
        x <- x[expr_genes, ]
        x
      })

    }


    if (!is.null(Rfile)) {
      saveRDS(filter_data,file = Rfile)
    }
    print(paste("Filter complete"))

    return(filter_data)
  }

  #merge sample
  merge_data_s5<-function(scRNAlist,merge_count=F,Rfile=NULL){

    merged_data <- merge(scRNAlist[[1]], y = scRNAlist[-1])

    if (merge_count==T) {
      merged_data[["RNA"]]<-JoinLayers(merged_data[["RNA"]])
    }
    if (!is.null(Rfile)) {
      saveRDS(merged_data,file = Rfile)
    }

    print("Merge complete")
    return(merged_data)
  }

  #SCT normlizeation and intergrate and batch remove
  sct_data_s6<-function(scRNAdata, all.gene = FALSE,cores=1 , Rfile = NULL, nfeatures = 3000,
                        reduction = c("cca", "rpca", "jpca", "rlsi"), future.seed = T, k.weight = 100,...){
    sct_data<-mclapply(filter_data, function(x){
      SCTransform(
        x,
        method = "glmGamPoi",
        return.only.var.genes = T,
        vars.to.regress = c("MT_percent","HB_percent"),
        conserve.memory = TRUE
      )
    },mc.cores = cores)
    gc()

    features <- SelectSCTIntegrationFeatures(
      object.list = sct_data,
      nfeatures = nfeatures
    )

    sct_data <- PrepSCTIntegration(
      object.list = sct_data,
      anchor.features = features
    )
    gc()

    anchors <- FindIntegrationAnchors(
      object.list = sct_data,
      normalization.method = "SCT",
      anchor.features = features,
      reduction = reduction
    )
    gc()

    sct_data <- IntegrateData(
      anchors,
      normalization.method = "SCT"
    )

    if (!is.null(Rfile)) {
      saveRDS(sct_data, file = Rfile)
    }
    print("normlization complete")
    return(sct_data)
  }

  #pca
  pca_reduction_s7<-function(scRNAdata,n.dims=50,Rfile=NULL){
    scRNAdata<-RunPCA(scRNAdata,npcs=n.dims)

    pca_data<-scRNAdata
    if (!is.null(Rfile)) {
      saveRDS(pca_data,file = Rfile)
    }
    print("pca complete")
    print(ElbowPlot(scRNAdata,ndims = 50))
    return(pca_data)
  }


  #reduction
  umap_reduction_s9<-function (scRNAdata, n.dims = 20, assay = "SCT", reduction = c("pca", "harmony"), learning_rate_tsne = 200,n_iter_tsne = 1000,perplexity_tsne = 30,
                               n.neighbors_umap = 30L, n.components = 3L, min.dist_umap = 0.3,Rfile = NULL,
                               ...)
  {
    scRNAdata <- RunUMAP(scRNAdata, dims = 1:n.dims, reduction = reduction, min.dist = min.dist_umap,
                         n.neighbors = n.neighbors_umap, n.components = n.components,
                         ...)
    scRNAdata <- RunTSNE(scRNAdata, dims = 1:n.dims, reduction = reduction, learning_rate = learning_rate_tsne,n_iter = n_iter_tsne,
                         dim.embed=n.components, perplexity = perplexity_tsne,
                         ...)
    print("Reduction complete")
    if (!is.null(Rfile)) {
      saveRDS(scRNAdata, file = Rfile)
    }
    return(scRNAdata)
  }

  umap_cluster_s10<-function(scRNAdata,n.dims=20,resolution,clustermethod=c("normal","recall"),assay="SCT",reduction=c("pca","harmony"),Rfile=NULL,algorithm = "louvain",core=4,...){
    #if resolution_start length >1 then use old method to clustree, otherwise use findclustersrecall
    print(paste("dims: ",n.dims))
    print(paste("resolution: ",resolution))

    scRNAdata <- FindNeighbors(scRNAdata, reduction = reduction, dims = 1:n.dims)

    if (clustermethod == "normal") {
      reduct_data <- FindClusters(scRNAdata, resolution = resolution)
      if (length(resolution) > 1) {
        clustree(reduct_data)
      }
    }
    else if (clustermethod == "recall") {
      reduct_data <- FindClustersRecall(scRNAdata, dims = 1:ndims,
                                        resolution_start = resolution, algorithm = "louvain",
                                        assay = "SCT", cores = core)
    }

    if (!is.null(Rfile)) {
      saveRDS(reduct_data,file = Rfile)
    }

    return(reduct_data)
  }

  #reduction plot
  clust_plot_s10<-function(scRNAdata,ident="seurat_clusters",...){
    Idents(scRNAdata)<-ident
    umap_integrated_1 <- DimPlot(scRNAdata,reduction = 'umap',group.by = 'orig.ident')
    umap_integrated_2 <- DimPlot(scRNAdata,reduction = 'tsne', label = T,...)
    umap_integrated_3 <- DimPlot(scRNAdata,reduction = 'umap', label = T,...)
    umap_integrated_4 <- DimPlot(scRNAdata,reduction = 'pca', label = T,...)
    integrated_plot <- list(batch_check=umap_integrated_1,pca=umap_integrated_4,tsne=umap_integrated_2,umap=umap_integrated_3)
    return(integrated_plot)
  }

  #find markers
  marker_find_s11<-function(scRNAdata,node=4,only.pos = T,assay = "SCT",slot="scale.data",test.use = "wilcox",Rfile=NULL,...){

    if (node>1) {
      options(future.globals.maxSize = node * 1024^2)
      plan("multisession", workers = node)
      all.markers <- FindAllMarkers(scRNAdata,...)
      plan("sequential")
    } else {
      all.markers <- FindAllMarkers(scRNAdata,only.pos = only.pos,assay = assay,slot=slot,test.use = test.use,...)
    }
    gc()
    if (!is.null(Rfile)) {
      saveRDS(all.markers,file = Rfile)
    }
    print("Find marker complete")
    return(all.markers)
  }

  #cluster annotation

  cluster_anno_s12<-function(de_data,scRNAdata,rowname="celltype",database=NULL,Rfile=NULL,return_marker=T,...){

    if (is.null(database)) {
      anno_data <- scMayoMap(data = de_data, ...)
    } else {
      database$value<-1
      colnames(database)<-c("gene","celltype","value")
      db <- tidyr::spread(demodata, key = c('celltype'), value = 'value')
      anno_data <- scMayoMap(data = de_data, database = databse,...)
    }

    tmp<-anno_data$markers %>%
      group_by(cluster) %>%
      slice_max(order_by = score, n = 1)
    markers<-tmp
    markers$gene<-str_split(tmp$genes,",")
    markers$cluster<-as.numeric(markers$cluster)
    markers<-markers[order(markers$cluster),]
    anno_data$markerlist<-markers
    tmp<-tmp[match(scRNAdata$seurat_clusters,tmp$cluster),]
    tmp$celltype<-factor(tmp$celltype,as.character(unique(tmp$celltype))[order(as.numeric(as.character(unique(tmp$cluster))))])
    scRNAdata[[rowname]]<-tmp$celltype


    if (return_marker==T) {
      data<-list(scRNAdata=scRNAdata,marker=anno_data)
    }else{
      data<-scRNAdata
    }
    if (!is.null(Rfile)) {
      saveRDS(data,file = Rfile)
    }
    print("annotation complete")
    return(data)
  }


  return(list(env_load_s1=env_load_s1,file_read_s2=file_read_s2,feature_show_s3=feature_show_s3,cell_filter_s4=cell_filter_s4,merge_data_s5=merge_data_s5,
              sct_data_s6=sct_data_s6,norm_data_s6=norm_data_s6,sct_data_new_s6=sct_data_new_s6,pca_reduction_s7=pca_reduction_s7,batch_rm_s8=batch_rm_s8,umap_reduction_s9=umap_reduction_s9,umap_cluster_s10=umap_cluster_s10,clust_plot_s10=clust_plot_s10,
              marker_find_s11=marker_find_s11,cluster_anno_s12=cluster_anno_s12))
}
