#' Title
#'
#' @param path the file path of single cell rna seq, every sample use the sample name as folder name
#' @param specise indicate the specise, human or mouse
#' @param min.cells minimal cell number of sample
#' @param min.features minimal feature number of per cell
#' @param n.dims the dims of pca, default as 20
#' @param resolution the resolution of umap, default as 0.6. if you want to specific it according to the cuttree, you can use interactive_mode to assign it at specific step.
#' @param assay default matrix to use. default as SCT, you can alter it to RNA.
#' @param interactive_mode if you want to assign cell filter parament or resolution you can alter it as TURE, or it will go through the whole pipline using default parament for all steps.
#'
#' @return
#' @export
#'
#' @examples
scRNA_demo <- function(path, specise = "human", min.cells = 3, min.features = 200,
                         n.dims = 20, resolution = seq(0.1, 1.0, 0.1), assay = "SCT",
                         interactive_mode = FALSE) {

  #enviroment load
  {
    library(Seurat)
    library(stringr)
    library(purrr)
    library(ggplot2)
    library(clustree)
    library(patchwork)
  }

  #count load
  file_read<-function(path,specise="human",min.cells=min.cells,min.features=min.features){

    ##default the input path contain multi sample folder, each sample folder have named accordant with sample info file, every sample folder contain the three basic files,
    dir_name <<- list.dirs(path,recursive = F)
    scRNAlist <- list()

    for (i in dir_name) {
      counts <- Read10X(data.dir = i)
      scRNAlist[[i]] <- CreateSeuratObject(counts = counts,
                                           project = str_split(i,"\\\\|/",simplify = T)%>%.[,ncol(.)],
                                           min.cells = min.cells,
                                           min.features = min.features)
    }
    names(scRNAlist)<-str_split(dir_name,"\\\\|/",simplify = T)%>%.[,ncol(.)]

    ## feature culculate
    if(specise=="human"){
      pattern<-"^MT-"
      HB_genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
    } else if(specise=="mouse"){
      pattern<-"^mt-"
      HB_genes <- c("Hba1","Hba2","Hbb","Hbd","Hbe1","Hbg1","Hbg2","Hbm","Hbq1","Hbz")
    }
    for (i in 1:length(scRNAlist)) {
      sc <- scRNAlist[[i]]
      sc[['mt_percent']] <- PercentageFeatureSet(sc,pattern = pattern)
      HB_m <- match(HB_genes,rownames(sc@assays$RNA))
      HB_genes <- rownames(sc@assays$RNA)[HB_m]
      HB_genes <- HB_genes[!is.na(HB_genes)]
      sc[['HB_percent']] <- PercentageFeatureSet(sc,features = HB_genes)
      scRNAlist[[i]] <- sc
      rm(sc)
    }
    print(paste("Totally ",length(scRNAlist)," sample loaded",sep = ""))
    return(scRNAlist)

  }

  #parament select
  feature_show<-function(scRNAlist){
    violin_before<-list()
    for (i in 1:length(scRNAlist)){
      violin_before[[i]] <- VlnPlot(scRNAlist[[i]],
                                    features = c('nFeature_RNA','nCount_RNA','mt_percent','HB_percent'),
                                    pt.size = 0.01,
                                    ncol = 4)+ theme(axis.text.x = element_text(angle = 0))+xlab("")
    }
    p<-reduce(violin_before, `/`)
    return(p)
  }

  #cell filter
  cell_filter<-function(scRNAlist,interactive_mode=interactive_mode){
    get_input <- function(prompt, default) {
      input <- readline(prompt)
      if (input == "") default else as.numeric(input)
    }

    print("cell filter parament select")
    nFeature_RNA_min <- ifelse(interactive_mode == F,200,get_input("nFeature_RNA 最小值 (默认 200): ", 200))
    nFeature_RNA_max <- ifelse(interactive_mode == F,4000,get_input("nFeature_RNA 最大值 (默认 4000): ", 4000))
    mt_percent_max <- ifelse(interactive_mode == F,10,get_input("mt_percent 最大值 (默认 10): ", 10))
    HB_percent_max <- ifelse(interactive_mode == F,5,get_input("HB_percent 最大值 (默认 5): ", 5))
    nCount_RNA_min <- ifelse(interactive_mode == F,500,get_input("nCount_RNA 最小值 (默认 500): ", 1000))
    nCount_RNA_max <- ifelse(interactive_mode == F,20000,get_input("nCount_RNA 最大值 (默认 20000): ", 20000))

    scRNAlist <- lapply(scRNAlist, function(x) {
      subset(x,
             nFeature_RNA > nFeature_RNA_min & nFeature_RNA < nFeature_RNA_max &
               mt_percent < mt_percent_max & HB_percent < HB_percent_max &
               nCount_RNA > nCount_RNA_min & nCount_RNA < nCount_RNA_max)
    })
    print(paste("Filter complete"))
    return(scRNAlist)
  }

  #cell reduction
  reduction<-function(scRNAlist,n.dims=20,resolution=0.6,assay="SCT"){
    assay <- if (assay != "SCT") {
      DefaultAssay(scRNAlist) <- "SCT"
    } else {
      "SCT"
    }
    scRNAlist <- FindNeighbors(scRNAlist,reduction = 'harmony',dims = 1:n.dims)
    scRNAlist <- FindClusters(scRNAlist,resolution = resolution)
    scRNAlist <- RunUMAP(scRNAlist,dims = 1:n.dims,reduction = 'harmony')
    scRNAlist <- RunTSNE(scRNAlist,dims = 1:n.dims,reduction = 'harmony')
    print("Reduction complete")
    print(paste("dims: ",n.dims))
    print(paste("resolution: ",resolution))
    return(scRNAlist)
  }

  #clust plot
  clust_plot<-function(scRNAlist){
    umap_integrated_1 <- DimPlot(scRNAlist,reduction = 'umap',group.by = 'orig.ident')
    umap_integrated_2 <- DimPlot(scRNAlist,reduction = 'tsne', label = T)
    umap_integrated_3 <- DimPlot(scRNAlist,reduction = 'umap', label = T)
    umap_integrated_4 <- DimPlot(scRNAlist,reduction = 'pca', label = T)
    integrated_plot <- umap_integrated_1+umap_integrated_2+umap_integrated_3
    return(integrated_plot)
  }

  # Helper function: Asks the user if they want to proceed, returns default if not in interactive mode
  ask_user <- function(prompt, default = TRUE) {
    if (!interactive_mode) return(default)  # If not in interactive mode, return default
    user_input <- readline(paste0(prompt, " (default ", ifelse(default, "Y", "N"), "): "))
    if (user_input == "") return(default)
    return(tolower(user_input) == "y")
  }

  # Step 1: Data Loading
  print(path)
  raw_data <- file_read(path, specise = specise, min.cells = min.cells, min.features = min.features)
  print(raw_data)

  # Step 2: Plotting Raw Features
  feature_raw <- feature_show(raw_data)
  print(feature_raw)

  # Step 3: Data Filtering
  if (interactive_mode == F) {
    qc_data <- cell_filter(raw_data,interactive_mode = interactive_mode)
    print(qc_data)
    # Step 4: Plotting Filtered Features
    feature_filter <- feature_show(qc_data)
    print(feature_filter)
  } else {
    error<-T
    while (error==T) {
      tryCatch({
        repeat {
          if (ask_user("run the data filtering step?", default = TRUE)){
            qc_data <- cell_filter(raw_data,interactive_mode = interactive_mode)
            print(qc_data)
            # Step 4: Plotting Filtered Features
            feature_filter <- feature_show(qc_data)
            print(feature_filter)
            break
          }
          else {
            error=FALSE
            break
          }
        }
      },error=function(e){
        e<-TRUE
        print("Error execute, reinput")
        return(e)
      })
    }
  }

  # Step 5: Sample Merging
  merge_data <- merge(x = qc_data[[1]], y = qc_data[-1],
                      add.cell.ids = str_split(dir_name, "\\\\|/", simplify = TRUE) %>% .[, ncol(.)])
  print(merge_data)
  print("Merge complete")

  # Step 6: Normalization
  merge_data <- SCTransform(merge_data, vars.to.regress = c('mt_percent', "HB_percent"))
  print("Normalization completed")

  # Step 7: PCA
  merge_data <- RunPCA(merge_data)
  print("PCA completed.")

  # Step 8: Batch Effect Correction
  norm_data <- IntegrateLayers(object = merge_data,
                               method = HarmonyIntegration,  # Method can be changed as needed
                               orig.reduction = 'pca',  # Must use PCA here
                               new.reduction = 'harmony')
  norm_data[['RNA']] <- JoinLayers(norm_data[['RNA']])
  pca_deviation <- ElbowPlot(norm_data, ndims = 50)
  print(pca_deviation)

  # Step 9: Dimensionality Reduction and Clustering
  if (length(resolution)!=1) {
    tmp <- reduction(norm_data, n.dims = n.dims, resolution = resolution, assay = assay)
    resolution_select <- clustree(tmp)
    print(resolution_select)
  }
  if (interactive_mode == F) {
    if (length(resolution)!=1) {
      resolution <- 0.6
    }
    norm_data <- reduction(norm_data, resolution = resolution)
    clust_perform <- clust_plot(norm_data)
    print(clust_perform)
  } else {
    error<-T
    while (error==T) {
      tryCatch({
        repeat {
          if (ask_user("Rerun clustering and dimensionality reduction step?", default = FALSE)) {
            resolution <- as.numeric(readline("Input a resolution parament: "))
            print(resolution)
            norm_data <- reduction(norm_data, resolution = resolution)
            # Step 10: Final Clustering Plot
            clust_perform <- clust_plot(norm_data)
            print(clust_perform)
          } else {
            error=FALSE
            break
          }
        }
      },error=function(e){
        e<-TRUE
        print("Resolution error, reinput")
        return(e)
      })
    }

  }
  sample_batch <- DimPlot(norm_data, reduction = 'umap', group.by = 'orig.ident') +
    ggtitle('harmony')
  print(sample_batch)

  # Return processed data and plots
  return(list(data = list(raw_data = raw_data, qc_data = qc_data, norm_data = norm_data),
              plot = list(feature_raw = feature_raw, feature_filter = feature_filter,
                          pca_deviation = pca_deviation, resolution_select = resolution_select,
                          sample_batch = sample_batch, clust_perform = clust_perform)))
}


