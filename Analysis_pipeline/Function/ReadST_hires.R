library(Seurat)
library(png)
library(jsonlite)
library(ggplot2)
Read10X_Image_edit <- function (image.dir, 
                                filter.matrix = TRUE, image.name) {
  image <- readPNG(source = file.path(image.dir, image.name))
  scale.factors <- fromJSON(txt = file.path(image.dir, "scalefactors_json.json"))
  tissue.positions <- read.csv(file = file.path(image.dir, 
                                                "tissue_positions_list.csv"), col.names = c("barcodes", 
                                                                                            "tissue", "row", "col", "imagerow", "imagecol"), header = FALSE, 
                               as.is = TRUE, row.names = 1)
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(x = tissue.positions$tissue == 
                                                 1), , drop = FALSE]
  }
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * 
    scale.factors$tissue_hires_scalef
  spot.radius <- unnormalized.radius/max(dim(x = image))
  return(new(Class = "VisiumV1", image = image, scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                                                             fiducial = scale.factors$fiducial_diameter_fullres, hires = scale.factors$tissue_hires_scalef, 
                                                                             scale.factors$tissue_hires_scalef), coordinates = tissue.positions, 
             spot.radius = spot.radius))
}


STPreProcess_edit <- function (InDir = InDir, Sample = Sample, OutDir = NULL) {
  if (is.null(OutDir) == TRUE) {
    OutDir = paste(getwd(), "/", Sample, "/", sep = "")
    dir.create(OutDir)
  }
  aa_try <- try(Xdata <- Seurat::Read10X(data.dir = paste(InDir, 
                                                          "filtered_feature_bc_matrix", sep = ""),gene.column = 1), silent = TRUE)
  if (is(aa_try, "try-error")) {
    Xdata <- Seurat::Read10X_h5(filename = paste(InDir, "filtered_feature_bc_matrix.h5", 
                                                 sep = ""))
  }
  else {
    Xdata <- Xdata
  }
  XF <- CreateSeuratObject(counts = Xdata, project = Sample, 
                           min.spots = 0, assay = "Spatial")
  Ximage <- Read10X_Image_edit(image.dir = paste(InDir, "spatial",sep = ""),image.name = "tissue_hires_image.png")
  Seurat::DefaultAssay(Ximage) <- "Spatial"
  Ximage <- Ximage[colnames(XF)]
  XF[["image"]] <- Ximage
  TumorST <- XF
  dir.create(paste(OutDir, "QC", sep = ""))
  TumorST[["Mito.percent"]] <- PercentageFeatureSet(TumorST, 
                                                    pattern = "^MT-")
  pdf(paste(OutDir, "QC/Vlnplot.pdf", sep = ""), width = 6, 
      height = 4)
  p <- VlnPlot(TumorST, features = c("nFeature_Spatial", "nCount_Spatial" 
  ), pt.size = 0, combine = F)
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoLegend() + theme(axis.title.x = element_blank(), 
                                          axis.text.x = element_text(angle = 0))
  }
  p <- cowplot::plot_grid(plotlist = p, ncol = 3)
  print(p)
  dev.off()
  pdf(paste(OutDir, "QC/featurplot.pdf", sep = ""), width = 7, 
      height = 7)
  p <- SpatialFeaturePlot(TumorST, features = c("nFeature_Spatial", 
                                                "nCount_Spatial", "Mito.percent"), combine = F)
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + theme(axis.title.x = element_blank(), 
                             axis.text.x = element_text(angle = 0))
  }
  print(cowplot::plot_grid(plotlist = p, ncol = 3))
  dev.off()
  QCData <- TumorST@meta.data[, c("nCount_Spatial", "nFeature_Spatial", 
                                  "Mito.percent")]
  openxlsx::write.xlsx(QCData, paste(OutDir, "QC/QCData.xlsx", 
                                     sep = ""), overwrite = TRUE)
  return(TumorST)
}