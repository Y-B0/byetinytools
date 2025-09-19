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

  #SCT normlize
  sct_data_s6<-function(scRNAdata,copy2RNA=T,all.gene=TRUE,Rfile=NULL,node=10,...){
    options(future.globals.maxSize = Inf)
    sct_data <- SCTransform(scRNAdata,method = "glmGamPoi",return.only.var.genes=!all.gene,vars.to.regress = c('MT_percent',"HB_percent"),...)
    sct_data <- PrepSCTFindMarkers(sct_data)
    if (copy2RNA==T) {
      sct_data[['RNA']] <- JoinLayers(sct_data[['RNA']])
      sct_data@assays$RNA$data<-sct_data@assays$SCT$data
      sct_data@assays$RNA$scale.data<-sct_data@assays$SCT$scale.data
    }
    if (!is.null(Rfile)) {
      saveRDS(sct_data,file = Rfile)
    }
    print("normlization complete")

    return(sct_data)
  }

  sct_data_new_s6<-function(scRNAdata,all.gene=TRUE,Rfile=NULL,node=10,nfeatures=3000,...){
    options(future.globals.maxSize = Inf)
    sct_data<-lapply(scRNAdata,SCTransform,method = "glmGamPoi",return.only.var.genes=!all.gene,vars.to.regress = c('MT_percent',"HB_percent"),...)
    features <- SelectIntegrationFeatures(sct_data, nfeatures = nfeatures)
    sct_data <- PrepSCTIntegration(sct_data, anchor.features = features)
    anchors <- FindIntegrationAnchors(object.list = sct_data, normalization.method = "SCT",
                                      anchor.features = features)
    sct_data <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
    if (!is.null(Rfile)) {
      saveRDS(sct_data,file = Rfile)
    }
    print("normlization complete")

    return(sct_data)
  }


  norm_data_s6<-function(scRNAdata,vars.to.regress=c('MT_percent',"HB_percent"),Rfile=NULL){

    scRNAdata <- NormalizeData(scRNAdata)
    scRNAdata <- FindVariableFeatures(scRNAdata,nFeature_RNA = nFeature_RNA)
    scRNAdata <- ScaleData(scRNAdata,vars.to.regress = vars.to.regress)

    if (!is.null(Rfile)) {
      saveRDS(sct_data,file = Rfile)
    }
    print("normlization complete")

    return(scRNAdata)
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

  #batch effect remove
  batch_rm_s8<-function(scRNAdata,Rfile=NULL){
    batchrm_data <- IntegrateLayers(object = scRNAdata,
                                 method = HarmonyIntegration,  # Method can be changed as needed
                                 orig.reduction = 'pca',  # Must use PCA here
                                 new.reduction = 'harmony')
    pca_deviation <- ElbowPlot(batchrm_data, ndims = 50)
    print(pca_deviation)

    if (!is.null(Rfile)) {
      saveRDS(batchrm_data,file = Rfile)
    }
    print("remove batch effect complete")

    return(batchrm_data)
  }

  #reduction
  umap_reduction_s9<-function(scRNAdata,n.dims=20,resolution=seq(0,1,0.1),assay="SCT",reduction=c("pca","harmony"),Rfile=NULL,...){

    scRNAdata <- FindNeighbors(scRNAdata,reduction = reduction,dims = 1:n.dims)
    scRNAdata <- FindClusters(scRNAdata,resolution = resolution)
    scRNAdata <- RunUMAP(scRNAdata,dims = 1:n.dims,reduction = reduction,...)
    scRNAdata <- RunTSNE(scRNAdata,dims = 1:n.dims,reduction = reduction,...)
    print("Reduction complete")
    print(paste("dims: ",n.dims))
    print(paste("resolution: ",resolution))

    reduct_data<-scRNAdata
    if (!is.null(Rfile)) {
      saveRDS(reduct_data,file = Rfile)
    }
    if (length(resolution)>1) {
      clustree(reduct_data)
    }
    return(reduct_data)
  }

  #reduction plot
  clust_plot_s10<-function(scRNAdata,...){
    umap_integrated_1 <- DimPlot(scRNAdata,reduction = 'umap',group.by = 'orig.ident')
    umap_integrated_2 <- DimPlot(scRNAdata,reduction = 'tsne', label = T,...)
    umap_integrated_3 <- DimPlot(scRNAdata,reduction = 'umap', label = T,...)
    umap_integrated_4 <- DimPlot(scRNAdata,reduction = 'pca', label = T,...)
    integrated_plot <- umap_integrated_1+umap_integrated_2+umap_integrated_3
    return(integrated_plot)
  }

  #find markers
  marker_find_s11<-function(scRNAdata,node=4,only.pos = T,assay = "SCT",slot="scale.data",test.use = "wilcox",Rfile=NULL,multcore=F,...){

    if (multcore==T) {
      options(future.globals.maxSize = node * 1024^2)
      plan("multisession", workers = node)
    }

    all.markers <- FindAllMarkers(scRNAdata,...)

    if (multcore==T) {
      plan("sequential")
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

    tmp<-group_by(anno_data$markers,cluster)%>%arrange(desc(score),.by_group = T)%>%slice_head(.,n=1)
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
              sct_data_s6=sct_data_s6,norm_data_s6=norm_data_s6,sct_data_new_s6=sct_data_new_s6,pca_reduction_s7=pca_reduction_s7,batch_rm_s8=batch_rm_s8,umap_reduction_s9=umap_reduction_s9,clust_plot_s10=clust_plot_s10,
              marker_find_s11=marker_find_s11,cluster_anno_s12=cluster_anno_s12))
}
