.libPaths(c("/home/sjt/R/x86_64-pc-linux-gnu-library/4.2","/usr/local/lib64/R/library"))
library(Seurat)
library(magrittr)
library(reticulate)
library(Cottrazm)
library(mailR)

process_data <- function(){
  ## input ----
  args <- commandArgs(trailingOnly = TRUE)
  folder <- args[1]
  stFile <- args[2]
  scrnaMatFile <- args[3]
  scrnaMetadataFile <- args[4]
  TumorCluster_ori <- args[5]
  splitElements_tumor  <-  strsplit(TumorCluster_ori, "[, /]+")[[1]]
  TumorCluster <- trimws(splitElements_tumor)
  
  EpithelialCluster_ori <- args[6]
  splitElements_epi  <-  strsplit(EpithelialCluster_ori, "[, /]+")[[1]]
  EpithelialCluster <- trimws(splitElements_epi)
  
  StrmalCluster_ori <- args[7]
  splitElements_stromal  <-  strsplit(StrmalCluster_ori, "[, /]+")[[1]]
  StrmalCluster <- trimws(splitElements_stromal)
  
  email <- args[8]
  
  Sys.unsetenv("RETICULATE_PYTHON")
  use_condaenv("TumorBoundary", required = TRUE,conda = "/work/sjt/software/miniconda3/condabin/conda")
  use_python("/work/sjt/software/miniconda3/envs/TumorBoundary/bin/python", required = TRUE)
  print(paste0("folder=",folder))
  print(paste0("StrmalCluster=",StrmalCluster))
  ## pipeline ----
  ## 0. LOAD ----
  source('/work/sjt/SpatialTME/ST_pipeline/Function/ReadST_hires.R')
  source('/work/sjt/SpatialTME/ST_pipeline/Function/get_sig_exp_largeMat.R')
  MainDir <- folder
  # 获取解压后的目录名
  unzipped_dirs <- list.dirs(MainDir, full.names = TRUE, recursive = FALSE)
  # 假设您想要的目录是最新创建的一个
  InDir <- paste0(unzipped_dirs[length(unzipped_dirs)],"/")
  
  OutDir <- paste0(MainDir,"/outs/")
  dir.create(OutDir)
  Sample <- 'User'
  res = 1.5
  
  print(paste0("InDir=",InDir))
  
  aa_try <- try(Xdata <- Seurat::Read10X(data.dir = paste(InDir,
                                                          "filtered_feature_bc_matrix", sep = "")), silent = T)
  if (is(aa_try, "try-error")) {
    Xdata <- Seurat::Read10X_h5(filename = paste(InDir, "filtered_feature_bc_matrix.h5",
                                                 sep = ""))
  }else {
    Xdata <- Xdata
  }
  XF <- CreateSeuratObject(counts = Xdata, project = Sample,
                           min.spots = 0, assay = "Spatial")
  if (file.exists(paste0(InDir,'spatial','/tissue_lowres_image.png'))) {
    Ximage <- Read10X_Image(image.dir = paste(InDir, "spatial",
                                              sep = ""))
  }else{
    Ximage <- Read10X_Image_edit(image.dir = paste0(InDir,"/spatial"),image.name = "tissue_hires_image.png")
  }
  
  Seurat::DefaultAssay(Ximage) <- "Spatial"
  Ximage <- Ximage[colnames(XF)]
  XF[["image"]] <- Ximage
  TumorST <- XF
  TumorST[["Mito.percent"]] <- PercentageFeatureSet(TumorST,
                                                    pattern = "^MT-")
  rm(XF,Xdata,Ximage)
  
  ## 1. Morph adjust ----
  print(paste0("startSME"))
  
  if ("filtered_feature_bc_matrix.h5"%in%list.files(InDir)&&"tissue_lowres_image.png"%in%list.files(paste0(InDir,"spatial"))) {
    source_python(system.file("python/Rusedtile.py", package = "Cottrazm"))
  }else if ("filtered_feature_bc_matrix.h5"%in%list.files(InDir)==TRUE&"tissue_lowres_image.png"%in%list.files(paste0(InDir,"spatial"))==FALSE) {
    source_python("/work/sjt/SpatialTME/ST_pipeline/Function/Rusedtile_hires_h5file.py")
  }else if("filtered_feature_bc_matrix.h5"%in%list.files(InDir)==FALSE&"tissue_lowres_image.png"%in%list.files(paste0(InDir,"spatial"))==TRUE){
    source_python("/work/sjt/SpatialTME/ST_pipeline/Function/Rusedtile_lowres_mtx.py")
  }else if("filtered_feature_bc_matrix.h5"%in%list.files(InDir)==FALSE&"tissue_hires_image.png"%in%list.files(paste0(InDir,"spatial"))==TRUE){
    source_python("/work/sjt/SpatialTME/ST_pipeline/Function/Rusedtile_hires_mtx.py")
  }
  
  Adjusted_expr_mtx <- ME_normalize(inDir = InDir, outDir = OutDir,
                                    sample = Sample)
  aa_try <- try(rownames(Adjusted_expr_mtx) <- colnames(TumorST@assays$Spatial@counts),
                silent = TRUE)
  if (is(aa_try, "try-error")) {
    library(Matrix)
    Adjusted_expr_mtx <- Matrix::readMM(paste(OutDir, Sample,
                                              "_raw_SME_normalizeA.mtx", sep = ""))
    rownames(Adjusted_expr_mtx) <- colnames(TumorST@assays$Spatial@counts)
    colnames(Adjusted_expr_mtx) <- rownames(TumorST@assays$Spatial@counts)
  }else {
    rownames(Adjusted_expr_mtx) <- colnames(TumorST@assays$Spatial@counts)
    colnames(Adjusted_expr_mtx) <- rownames(TumorST@assays$Spatial@counts)
  }
  Adjusted_expr_mtxF <- t(as.matrix(Adjusted_expr_mtx))
  MorphMatirxSeurat <- CreateSeuratObject(counts = as.matrix(Adjusted_expr_mtxF))
  MorphMatirxSeurat <- subset(MorphMatirxSeurat, cells = rownames(TumorST@meta.data))
  TumorST@assays$Morph <- MorphMatirxSeurat@assays$RNA
  TumorST <- NormalizeData(TumorST, assay = "Morph")
  TumorST <- FindVariableFeatures(object = TumorST, mean.function = ExpMean,
                                  dispersion.function = LogVMR, assay = "Morph")
  TumorST <- ScaleData(object = TumorST, assay = "Morph")
  TumorST <- RunPCA(object = TumorST, npcs = 50, verbose = FALSE,
                    assay = "Morph")
  TumorST <- FindNeighbors(TumorST, reduction = "pca", dims = 1:50,
                           assay = "Morph")
  TumorST <- RunUMAP(object = TumorST, dims = 1:50, assay = "Morph")
  TumorST <- FindClusters(TumorST, resolution = res, algorithm = 1,
                          graph.name = "Morph_snn")
  TumorST@meta.data$seurat_clusters <- TumorST@meta.data[,
                                                         paste("Morph_snn_res.", res, sep = "")]
  
  rm(Adjusted_expr_mtx,Adjusted_expr_mtxF,MorphMatirxSeurat)
  ## 2. Deconvolution ----
  ### 2.1 scRNA process ----
  print(paste0("scRNA process"))
  
  if (grepl(".csv",scrnaMatFile)) {
    sc_mat <- read.csv(file = paste0(MainDir,"/",scrnaMatFile),row.names = 1)
  }else if(grepl(".rds",scrnaMatFile)){
    sc_mat <- readRDS(file = paste0(MainDir,"/",scrnaMatFile))
  }
  sc_metadata <- read.csv(paste0(MainDir,"/",scrnaMetadataFile),row.names = 1)
  colnames(sc_metadata) <- 'DeconTypes'
  
  sc_obj <- CreateSeuratObject(counts = sc_mat,
                               min.cells = 0,
                               min.features = 0,
                               project = "seurat_object")
  sc_obj <- AddMetaData(sc_obj,sc_metadata)
  Idents(sc_obj) <- 'DeconTypes'
  clustermarkers <-
    Seurat::FindAllMarkers(object = sc_obj,
                           logfc.threshold = 0.25,
                           only.pos = TRUE)
  clustermarkers_list <- split(clustermarkers, clustermarkers$cluster)
  clustermarkers_list <-
    lapply(names(clustermarkers_list), function(cluster) {
      sub_markers <- clustermarkers_list[[cluster]]$gene
    })
  names(clustermarkers_list) <-
    names(split(clustermarkers, clustermarkers$cluster))
  #sig_exp
  sig_exp <-
    get_sig_exp_largeMat(
      se.obj = sc_obj,
      DefineTypes = "DeconTypes",
      sig_scran = unique(unlist(clustermarkers_list))
    )
  
  rm(sc_obj,clustermarkers)
  ### 2.2 Enrichment analysis of ST data ----
  TumorST <- NormalizeData(TumorST, assay = "Spatial")
  TumorST@meta.data$Location <- 'total'
  TumorST@meta.data$Decon_topics <-
    paste(TumorST@meta.data$Location,
          TumorST@meta.data$seurat_clusters,
          sep = "_")
  expr_values = as.matrix(TumorST@assays$Spatial@data) #ST expr log
  nolog_expr = 2 ^ (expr_values) - 1 #ST expr nolog
  meta_data <-
    TumorST@meta.data[, c("nCount_Spatial", "Decon_topics","Location")]
  for (cluster in names(clustermarkers_list)) {
    cluster_markers = clustermarkers_list[[cluster]][1:25]
    cluster_score <-
      apply(TumorST@assays$Spatial@data[rownames(TumorST@assays$Spatial@data) %in% cluster_markers, ], 2, mean)
    meta_data <- cbind(meta_data, cluster_score)
  }
  colnames(meta_data) <-
    c("nCount_Spatial",
      "Decon_topics",
      'Location',
      names(clustermarkers_list))
  intersect_gene = intersect(rownames(sig_exp), rownames(nolog_expr))
  filter_sig = sig_exp[intersect_gene, ]
  filter_expr = nolog_expr[intersect_gene, ]
  filter_log_expr = expr_values[intersect_gene, ]
  enrich_matrix <-
    get_enrich_matrix(filter_sig = filter_sig,
                      clustermarkers_list = clustermarkers_list)
  enrich_result <-
    enrich_analysis(filter_log_expr = filter_log_expr,
                    enrich_matrix = enrich_matrix)
  rm(expr_values,nolog_expr,intersect_gene,filter_log_expr)
  ### 2.3 Spot deconvolution ----
  print(paste0("Decon process"))
  
  DeconData <- SpatialDecon(
    enrich_matrix = enrich_matrix,
    enrich_result = enrich_result,
    filter_expr = filter_expr,
    filter_sig = filter_sig,
    clustermarkers_list = clustermarkers_list,
    meta_data = meta_data,
    malignant_cluster = TumorCluster,
    tissue_cluster = EpithelialCluster,
    stromal_cluster = StrmalCluster
  )
  
  ### 2.4 Decon Filter 0.05 ----
  DeconData <- as.matrix(DeconData)
  DeconData[DeconData<0.05] <- 0
  DeconData <- as.data.frame(DeconData)
  DeconData[,2:ncol(DeconData)] <- apply(DeconData[,2:ncol(DeconData)],2,as.numeric)
  
  TumorST <- AddMetaData(TumorST,DeconData%>%as.data.frame%>%tibble::column_to_rownames('cell_ID'))
  TumorST.meta <- TumorST@meta.data[,setdiff(colnames(TumorST@meta.data),c("Location","Decon_topics","Morph_snn_res.1.5"))]
  colnames(TumorST.meta) <- gsub("seurat_clusters","MorphAdjustedCluster",colnames(TumorST.meta))
  write.csv(TumorST.meta,paste0(OutDir,'ST_results.csv'),row.names = TRUE,quote = FALSE)
  
  ## 3. Reconstruction ----
  print(paste0("Recon process"))
  
  TumorSTRecon <- SpatialRecon(TumorST = TumorST,
                               sig_exp = sig_exp,
                               clustermarkers_list = clustermarkers_list,
                               DeconData = DeconData,
                               Location = c("total"))
  Recon.mtx <- TumorSTRecon@assays$RNA@data
  Recon.mtx.df <- as.data.frame(Recon.mtx)
  saveRDS(Recon.mtx,paste0(OutDir,"Subspot_matrix.rds.gz"),compress = "gzip")
  
  ## 4. Email ----
  # send.mail(from = "yelab2020@gmail.com",
  #           to = email,
  #           subject = "SpatialTME: Your ST results",
  #           body = "Hi, thanks for using SpatialTME. Here are your ST results.",
  #           smtp = list(host.name = "smtp.gmail.com", port = 587,
  #                       user.name = "yelab2020@gmail.com", passwd = "", ssl = TRUE),
  #           authenticate = TRUE,
  #           send = TRUE,
  #           attach.files = c(paste0(OutDir,'ST_results.csv'), paste0(OutDir,"Subspot_matrix.rds.gz")))
  
  send.mail(from = "shijintong127@163.com",
            to = 'shijintong@sjtu.edu.cn',
            subject = "SpatialTME: Your ST results",
            body = "Hi, thanks for using SpatialTME. Here are your ST results.",
            smtp = list(host.name = "smtp.163.com", port = 587,
                        user.name = "shijintong127@163.com", passwd = "MWRQNHJIJOQWPYBT", ssl = TRUE),
            authenticate = TRUE,
            send = TRUE,
            attach.files = c(paste0(OutDir,'ST_results.csv'), paste0(OutDir,"Subspot_matrix.rds.gz")))
  ## 5. Delete ----
  unlink(MainDir, recursive = TRUE)
  print(paste0("Success"))
}

# 使用tryCatch来捕获任何错误
tryCatch({
  process_data()
}, error = function(err) {
  # 如果出错，发送错误信息
  args <- commandArgs(trailingOnly = TRUE)
  folder <- args[1]
  stFile <- args[2]
  scrnaMatFile <- args[3]
  scrnaMetadataFile <- args[4]
  TumorCluster_ori <- args[5]
  splitElements_tumor  <-  strsplit(TumorCluster_ori, "[, /]+")[[1]]
  TumorCluster <- trimws(splitElements_tumor)
  
  EpithelialCluster_ori <- args[6]
  splitElements_epi  <-  strsplit(EpithelialCluster_ori, "[, /]+")[[1]]
  EpithelialCluster <- trimws(splitElements_epi)
  
  StrmalCluster_ori <- args[7]
  splitElements_stromal  <-  strsplit(StrmalCluster_ori, "[, /]+")[[1]]
  StrmalCluster <- trimws(splitElements_stromal)
  
  email <- args[8]
  
  send.mail(from = "shijintong127@163.com",
            to = email,
            subject = "SpatialTME: Error in processing your data",
            body = paste0("An error occurred:", err$message),
            smtp = list(host.name = "smtp.163.com", port = 587,
                        user.name = "shijintong127@163.com", passwd = "MWRQNHJIJOQWPYBT", ssl = TRUE),
            authenticate = TRUE,
            send = TRUE)
  # 删除文件夹

  unlink(folder, recursive = TRUE)
})
