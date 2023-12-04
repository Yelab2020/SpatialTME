checkGeneList <- function(Genes, expr_data) {
  Genes <- toupper(Genes)  # All the genes input are converted to uppercase.
  genes_input <- unlist(strsplit(Genes, ",| |\n|\t"))
  genes_input <- genes_input[genes_input != ""]
  all_genes_avaliable_flag = 1
  right_gene <- c()
  wrong_gene <- c()
  obj_GeneList <- row.names(expr_data)
  
  for (gene in genes_input) {
    if (!(gene %in% obj_GeneList)) {
      all_genes_avaliable_flag = 0
      wrong_gene <- c(wrong_gene, gene)
    }else{
      right_gene <- c(right_gene, gene)
    }
  }
  GeneListInfo <- list(
    genes_input = toupper(genes_input),
    all_genes_avaliable_flag = all_genes_avaliable_flag,
    gene_number = length(right_gene),
    right_gene = right_gene,
    wrong_gene = wrong_gene
  )
  return(GeneListInfo)
}

## plot function ----
geom_spatial <-  function(
    mapping = NULL,
    data = NULL,
    image = image,
    image.alpha = image.alpha,
    crop = crop,
    stat = "identity",
    position = "identity",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE,
    ...
) {
  layer(
    geom = GeomSpatial,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, image = image, image.alpha = image.alpha, crop = crop, ...)
  )
}
SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
col_Decon <-c("EndothelialCells"="#ffffad","Fibroblasts"='#ff9833',
              "BCells"='#9c4da5','PlasmaCells'='#C9B3E3','Neutrophil'='#84E4BD','B/Plasma'='#9c4da5',
              "CD4TCells"='#F070C8',"CD8TCells"='#F5A1BD','NKCells'='#63c3a5','TCells'='#f070c8', 'MyeloidCells'='#396db5',
              "Macrophages"='#396db5',"DCs"='#84b2d6',
              "EpithelialCells"="#33a02c","TumorCells"="#e31a1c","CNSCells"="#33a02c","Melanocytes"='#c1e6a1',
              "Hepatocytes"="#33a02c","Melanocyte"='#c1e6a1',"AdiposeCells"='#FFE4C4')
col_location <- c("Mal"="#CB181D",'Bdy'= "#1f78b4",'nMal'="#fdb462")
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3',
               '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282',
               '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8',
               '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755',
               '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3',
               '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD',
               '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
               '#968175')
.cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE",
                   "#882E72", "#B17BA6", "#FF7F00", "#FDB462", "#E7298A",
                   "#E78AC3", "#33A02C", "#B2DF8A", "#55B1B1", "#8DD3C7",
                   "#A6761D", "#E6AB02", "#7570B3", "#BEAED4", "#666666",
                   "#999999", "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3",
                   "#808000", "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d",
                   "#ffff00")

c102 <- c("#FFFF00","#1CE6FF","#FF34FF","#FF4A46","#008941", #1
          "#006FA6","#A30059","#FFDBE5","#7A4900","#0000A6", #2
          "#63FFAC","#B79762","#004D43","#8FB0FF","#997D87", #3
          "#5A0007","#809693","#6A3A4C","#1B4400","#4FC601", #4
          "#3B5DFF","#4A3B53","#FF2F80","#61615A","#BA0900", #5
          "#6B7900","#00C2A0","#FFAA92","#FF90C9","#B903AA", #6
          "#D16100","#DDEFFF","#000035","#7B4F4B","#A1C299", #7
          "#300018","#0AA6D8","#013349","#00846F","#372101", #8
          "#FFB500","#C2FFED","#A079BF","#CC0744","#C0B9B2", #9
          "#C2FF99","#001E09","#00489C","#6F0062","#0CBD66", #10
          "#EEC3FF","#456D75","#B77B68","#7A87A1","#788D66", #11
          "#885578","#FAD09F","#FF8A9A","#D157A0","#BEC459", #12
          "#456648","#0086ED","#886F4C","#34362D","#B4A8BD", #13
          "#00A6AA","#452C2C","#636375","#A3C8C9","#FF913F", #14
          "#938A81","#575329","#00FECF","#B05B6F","#8CD0FF", #15
          "#3B9700","#04F757","#C8A1A1","#1E6E00","#7900D7", #16
          "#A77500","#6367A9","#A05837","#6B002C","#772600", #17
          "#D790FF","#9B9700","#549E79","#FFF69F","#201625", #18
          "#72418F","#BC23FF","#99ADC0","#3A2465","#922329", #19
          "#5B4534","#FDE8DC","#404E55","#0089A3","#CB7E98", #20
          "#A4E804","#324E72")                               #21

## BoundaryPlot  ----
DoBoundaryPlot <- function(structure.type,metadata,image.use,alpha,crop = TRUE,pt.size.factor,legend.size.factor){
  if (structure.type=='Morphological adjusted cluster') {
    group.by = 'Morph_clusters'
  }else if(structure.type=="Mal-Bdy-nMal axis"){
    group.by = 'Location'
  }
  images <-"image"
  data <- metadata[,group.by]%>%as.data.frame()
  rownames(data) <- rownames(metadata)
  colnames(data) <- group.by
  for (group in group.by) {
    if (!is.factor(x = data[, group])) {
      data[, group] <- factor(x = data[, group])
    }
  }
  features <- colnames(x = data)
  coordinates <- GetTissueCoordinates(object = image.use)
  if (structure.type=='Mal-Bdy-nMal axis') {
    plot <- SingleSpatialPlot(data = cbind(coordinates,data[rownames(x = coordinates), features,drop = FALSE]),
                              image = image.use, image.alpha = 1, 
                              col.by = features, cols = col_location, 
                              alpha.by = NULL, 
                              pt.alpha = alpha[1],
                              geom = "spatial", cells.highlight = NULL, cols.highlight = NULL, 
                              pt.size.factor = pt.size.factor, stroke = 0.25, 
                              crop = crop)+
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"))+
      guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor)))
  }else{
    plot <- SingleSpatialPlot(data = cbind(coordinates,data[rownames(x = coordinates), features,drop = FALSE]),
                              image = image.use, image.alpha = 1, 
                              col.by = features, 
                              alpha.by = NULL, 
                              pt.alpha = alpha[1],
                              geom = "spatial", cells.highlight = NULL, cols.highlight = NULL, 
                              pt.size.factor = pt.size.factor, stroke = 0.25, 
                              crop = crop)+
      scale_fill_manual(values = c(.cluster_cols,my36colors))+
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"))+
      guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor)))
  }
  
  return(plot)
}

## GeneExrPlot----
DoGeneExprPlot <- function(gene,metadata,expr_data,ncol,image.use,seperate.plot,combine=TRUE,pt.size.factor=1.6,stroke=0.25,crop=TRUE){
  if(length(gene)==1){
    data <- expr_data[gene,]%>%as.data.frame()
    colnames(data) <- gene
    features <- gene
    plots <- vector(mode = "list", length = length(x = features))
    coordinates <- GetTissueCoordinates(object = image.use)
    plot <- SingleSpatialPlot(data = cbind(coordinates, data[rownames(x = coordinates), features, drop = FALSE]),
                              image = image.use, image.alpha = 1, 
                              col.by = features, cols = NULL,
                              alpha.by = NULL,
                              pt.alpha = NULL,
                              geom = "spatial",
                              cells.highlight = NULL, cols.highlight = c("#DE2D26", "grey50"), 
                              pt.size.factor = pt.size.factor,
                              stroke = stroke, 
                              crop = crop) +
      scale_fill_gradientn(name = features, colours = SpatialColors(n = 100)) +
      theme(legend.position = "top") + scale_alpha(range = alpha) + guides(alpha = FALSE)
    plots<- plot
  }else{ #gene>2
    if (seperate.plot=='Seperate') {
      data <- expr_data[gene,]%>%as.matrix()%>%t%>%as.data.frame()
      features <- colnames(data)
      plots <- vector(mode = "list", length = length(x = features))
      coordinates <- GetTissueCoordinates(object = image.use)
      for (j in 1:length(x = features)) {
        plot <- SingleSpatialPlot(data = cbind(coordinates, data[rownames(x = coordinates), features[j], drop = FALSE]),
                                  image = image.use, image.alpha = 1, 
                                  col.by = features[j], cols = NULL,
                                  alpha.by = NULL,
                                  pt.alpha = NULL,
                                  geom = "spatial",
                                  cells.highlight = NULL, cols.highlight = c("#DE2D26", "grey50"), 
                                  pt.size.factor = pt.size.factor,
                                  stroke = stroke, 
                                  crop = crop) +
          scale_fill_gradientn(name = features[j], colours = SpatialColors(n = 100)) +
          theme(legend.position = "top") + scale_alpha(range = alpha) + guides(alpha = FALSE)
        plots[[j]] <- plot
      }
      plots <- wrap_plots(plots = plots, ncol = ncol)
    }else{
      data <- expr_data[gene,]%>%as.matrix()%>%t%>%as.data.frame()
      data[,"mean"] <- rowMeans(data)
      features <- "mean"
      plots <- vector(mode = "list", length = length(x = features))
      coordinates <- GetTissueCoordinates(object = image.use)
      plot <- SingleSpatialPlot(data = cbind(coordinates, data[rownames(x = coordinates), features, drop = FALSE]),
                                image = image.use, image.alpha = 1, 
                                col.by = features, cols = NULL,
                                alpha.by = NULL,
                                pt.alpha = NULL,
                                geom = "spatial",
                                cells.highlight = NULL, cols.highlight = c("#DE2D26", "grey50"), 
                                pt.size.factor = pt.size.factor,
                                stroke = stroke, 
                                crop = crop) +
        scale_fill_gradientn(name = features, colours = SpatialColors(n = 100)) +
        theme(legend.position = "top") + scale_alpha(range = alpha) + guides(alpha = FALSE)
      plots<- plot
    }
  }
  return(plots)
}
## DEValcanoplot ----
DoDEValcanoPlot <- function(DiffGenes,cluster,cut_off_pvalue,cut_off_logFC,n,legend.size.factor){
  LocationDiff <- DiffGenes[[cluster]]
  LocationDiff <- LocationDiff[LocationDiff$Symbol %in% grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL", 
                                                             LocationDiff$Symbol, invert = TRUE, value = TRUE), ]
  dataset <- LocationDiff %>% tibble::rownames_to_column() %>% 
    set_colnames(., c("gene", colnames(LocationDiff)))
  dataset$color <- ifelse(LocationDiff$FDR > cut_off_pvalue | 
                            abs(LocationDiff$Diff) <= cut_off_logFC, "grey", ifelse(LocationDiff$Diff > 
                                                                                      0, cluster, "Other"))
  if (min(dataset$pvalue)==0) {
    sorted_numbers <- sort(dataset$pvalue)
    dataset[which(dataset$pvalue==0),'pvalue'] <- unique(sorted_numbers)[2]
  }
  datasetN <- dataset %>% dplyr::group_by(color) %>% dplyr::top_n(n, 
                                                                  abs(Diff))
  datasetN <- datasetN[datasetN$color != "grey", ]
  dataset$color <- factor(dataset$color,levels = c(cluster,'grey','Other'))
  plot <- ggplot(dataset, aes(x = Diff, y = (-log10(pvalue)), 
                              color = color)) + geom_point(aes(fill = color), size = 1) + 
    scale_color_manual(values = c("red", "grey", "blue")) + 
    ggrepel::geom_text_repel(data = datasetN, aes(label = datasetN$Symbol),
                             max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                             size = 3*legend.size.factor, point.padding = unit(0.35, "lines"),
                             segment.color = "black", show.legend = FALSE,
                             color = "black", fontface = "bold") +
    geom_vline(xintercept = c(-cut_off_logFC, cut_off_logFC), lty = 4, col = "black", lwd = 0.8) + 
    geom_hline(yintercept = cut_off_pvalue, lty = 4, col = "black", 
               lwd = 0.8) + labs(x = "log2(FoldChange)", y = "-log10(pvalue)") + 
    theme(panel.background = element_rect(color = "black", fill = NA), legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = 9*legend.size.factor),legend.key.size = unit(1*legend.size.factor, "lines"),
          axis.text.y = element_text(size = 11*legend.size.factor),axis.text.x = element_text(size = 9*legend.size.factor),
          axis.title = element_text(size = 13*legend.size.factor))+
    guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor)))
  
  return(plot)
}
## DEGOplot ----
DoGOGene <- function(DiffGenes,cluster,cut_off_pvalue,cut_off_logFC){
  Location=cluster
  sub <- DiffGenes[[cluster]]
  sub <- sub[sub$Diff >= cut_off_logFC & sub$FDR <= cut_off_pvalue, ]
  sub <- sub[sub$Symbol %in% grep("^IG[HJKL]|^RNA|^MT-|^RPS|^RPL", 
                                  sub$Symbol, invert = TRUE, value = TRUE), ]
  sub <- sub[!is.na(sub$Diff), ]
  
  LocationDiffFeatures <- sub$Symbol
  return(LocationDiffFeatures)
}
DoGOres <- function(DiffGenes){
  GO <- clusterProfiler::enrichGO(gene = DiffGenes, keyType = "SYMBOL", 
                                  OrgDb = "org.Hs.eg.db", ont = "BP", pAdjustMethod = "fdr", 
                                  pvalueCutoff = 0.2) %>%as.data.frame()
  return(GO)
}

DoGOPlot <- function(GO,n,legend.size.factor){
  if (n<nrow(GO)) {
    plot <- ggplot(GO[1:n,],aes(x = -log10(p.adjust),y = reorder(Description,-log10(p.adjust)),fill=Count))+
      geom_bar(stat = "identity")+
      scale_fill_gradientn(colors = brewer.pal(9,"Reds")[2:8])+
      ylab("")+theme_cowplot()+
      theme(legend.text = element_text(size = 9*legend.size.factor),legend.key.size = unit(1*legend.size.factor, "lines"),
            legend.title = element_text(size = 13*legend.size.factor),
            axis.text.y = element_text(size = 11*legend.size.factor),axis.text.x = element_text(size = 9*legend.size.factor),
            axis.title = element_text(size = 13*legend.size.factor))+scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30) )
  }else{
    plot <- ggplot(GO,aes(x = -log10(p.adjust),y = reorder(Description,-log10(p.adjust)),fill=Count))+
      geom_bar(stat = "identity")+
      scale_fill_gradientn(colors = brewer.pal(9,"Reds")[2:8])+
      ylab("")+theme_cowplot()+
      theme(legend.text = element_text(size = 9*legend.size.factor),legend.key.size = unit(1*legend.size.factor, "lines"),
            legend.title = element_text(size = 13*legend.size.factor),
            axis.text.y = element_text(size = 11*legend.size.factor),axis.text.x = element_text(size = 9*legend.size.factor),
            axis.title = element_text(size = 13*legend.size.factor))+scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 30) )
  }
  return(plot)
}

## DeconPlot-slide----
# DoDeconPlot_slide <- function(metadata,image.use,legend.size.factor){
#   metadata.use <- metadata
#   data <- metadata.use[,c("first_type","Morph_clusters")]%>%as.data.frame()
#   rownames(data) <- rownames(metadata.use)
#   colnames(data) <- c("CellType","Morph_clusters")
#   features <- 'CellType'
#   coordinates <- GetTissueCoordinates(object = image.use)
#   
#   plot <- SingleSpatialPlot(data = cbind(coordinates,data[rownames(x = coordinates), features,drop = FALSE]),
#                             image = image.use, image.alpha = 1, 
#                             col.by = features, cols = col_Decon[unique(data$CellType)], 
#                             alpha.by = NULL, 
#                             geom = "spatial", cells.highlight = NULL, cols.highlight = NULL, 
#                             pt.size.factor = 1, stroke = 0.25)+
#     theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
#           legend.key.size = unit(1*legend.size.factor, "lines"))+
#     guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor)))
#   return(plot)
# }
## DeconBarplot ----
DoDeconBarPlot <- function(structure,metadata,legend.size.factor){
  if (structure=='Morphological adjusted cluster') {
    group <- 'Morph_clusters'
  }else{
    group <- 'Location'
  }
  metadata_use <- metadata[,c(group,colnames(metadata)[6:ncol(metadata)])]
  plot_col <- colnames(metadata_use)[2:ncol(metadata_use)]
  Structure <- do.call(rbind, lapply(unique(metadata_use[,1]), 
                                     function(x) {
                                       metadata_split <- metadata_use[metadata_use[,1] == x, 
                                       ] %>% as.data.frame()
                                       sub <- colSums(metadata_split[, plot_col], na.rm = TRUE) %>% 
                                         data.frame() %>% tibble::rownames_to_column() %>% 
                                         set_colnames(., c("Types", "Sum")) %>% dplyr::mutate(Per = 100 * 
                                                                                                Sum/sum(Sum))
                                     }))
  Structure$Group <- rep(unique(metadata_use[,1]), each = length(plot_col))
  Structure$Types <- factor(Structure$Types, levels = plot_col)
  if (structure=='Morphological adjusted cluster') {
    Barplot <- ggplot(Structure, aes(x = Group, y = Per, fill = Types)) + 
      geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values = col_Decon[plot_col])+
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"),
            axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
      guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))
  }else{
    Barplot <- ggplot(Structure, aes(x = Group, y = Per, fill = Types)) + 
      geom_bar(stat = "identity") + theme_bw() + scale_fill_manual(values = col_Decon[plot_col])+
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"))+
      guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))
  }
  return(Barplot)
}
DoDeconBarPlot_slide <- function(metadata,legend.size.factor){
  metadata_use <- metadata[,c('Morph_clusters','first_type')]
  colnames(metadata_use) <- c('Group','CellType')
  p <- 
    ggplot(metadata_use, aes(x = Group, fill = CellType)) +
    geom_bar(position = "fill") +
    theme_bw() + scale_fill_manual(values = col_Decon[unique(metadata_use$CellType)])+
    theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
          legend.key.size = unit(1*legend.size.factor, "lines"),
          axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
    guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))
  return(p)
}

## CellLocationplot ----
DoCellLocationPlot <- function(structure,metadata,image.use,celltype,location,pt.size.factor=1){
  images <-"image"
  if (structure=='Morphological adjusted cluster') {
    group <- 'Morph_clusters'
  }else{
    group <- 'Location'
  }
  metadata.use <- metadata[,c(group,colnames(metadata)[6:ncol(metadata)])]
  colnames(metadata.use)[1] <- 'Structure'
  ##choose spots
  if (location=='All') {
    data <- metadata[,celltype]%>%as.data.frame()
    rownames(data) <- rownames(metadata)
    colnames(data) <- celltype
    features <- colnames(x = data)
    coordinates <- GetTissueCoordinates(object = image.use)
    
  }else{
    metadata.use <- metadata.use %>% dplyr::filter(Structure==location)
    data <- metadata.use[,celltype]%>%as.data.frame()
    rownames(data) <- rownames(metadata.use)
    colnames(data) <- celltype
    features <- colnames(x = data)
    coordinates <- GetTissueCoordinates(object = image.use)[rownames(metadata.use),]
  }
  
  plot <- SingleSpatialPlot(data = cbind(coordinates, data[rownames(x = coordinates), features, drop = FALSE]),
                            image = image.use, image.alpha = 1, 
                            col.by = features, cols = NULL,
                            pt.alpha = NULL,
                            # alpha=c(0,1),
                            geom = "spatial",
                            cells.highlight = NULL, cols.highlight = c("#DE2D26", "grey50"), 
                            pt.size.factor = pt.size.factor,
                            stroke = 0.25, 
                            crop = FALSE) +
    scale_fill_gradientn(name = features, colours = SpatialColors(n = 100)) +
    theme(legend.position = "top") + scale_alpha(range = c(0,1))
  
  return(plot)
  
}
DoCellLocationPlot_slide <- function(metadata,image.use,celltype,location,pt.size.factor=1,legend.size.factor=1){
    metadata.use <- metadata
    data <- metadata.use[,c("first_type","Morph_clusters")]%>%as.data.frame()
    rownames(data) <- rownames(metadata.use)
    colnames(data) <- c("CellType","Morph_clusters")
    features <- 'CellType'
    if (location=='All'&celltype=='All') {
      coordinates <- GetTissueCoordinates(object = image.use)
    }else if(location=='All'&celltype!='All'){
      data <- data%>%dplyr::filter(CellType==celltype)
      coordinates <- GetTissueCoordinates(object = image.use)[rownames(data),]
    }else if(location!='All'&celltype=='All'){
      data <- data%>%dplyr::filter(Morph_clusters==location)
      coordinates <- GetTissueCoordinates(object = image.use)[rownames(data),]
    }else{
      data <- data%>%dplyr::filter(Morph_clusters==location&CellType==celltype)
      coordinates <- GetTissueCoordinates(object = image.use)[rownames(data),]
    }
    plot <- SingleSpatialPlot(data = cbind(coordinates,data[rownames(x = coordinates), features,drop = FALSE]),
                              image = image.use, image.alpha = 1,
                              col.by = features, cols = col_Decon[unique(data$CellType)],
                              alpha.by = NULL,
                              geom = "spatial", cells.highlight = NULL, cols.highlight = NULL,
                              pt.size.factor = 1, stroke = 0.25)+
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"))+
      guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor)))
    return(plot)
}
## CellLocationBoxplot ----
DoCellLocationBoxPlot <- function(structure,metadata,celltype,legend.size.factor,test){
  if (structure=='Morphological adjusted cluster') {
    group <- 'Morph_clusters'
    metadata.use <- metadata[,c(group,colnames(metadata)[6:ncol(metadata)])]
    colnames(metadata.use)[1] <- 'Structure'
    data <- metadata.use[,c('Structure',celltype)]%>%as.data.frame()
    if (test=='Kruskal-Wallis') {
      plot <- ggplot(data,aes(x = Structure,y = get(celltype),fill=Structure))+
        geom_boxplot()+theme_bw()+
        scale_fill_manual(values = c(.cluster_cols,my36colors))+
        ylab(celltype)+
        theme(legend.title = element_text(size = 11*legend.size.factor),
              legend.text = element_text(size = 9*legend.size.factor),
              axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
        theme(legend.position = 'none')+
        guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))+
        stat_compare_means(method = "kruskal.test")
    }else if(test=='Anova'){
      plot <- ggplot(data,aes(x = Structure,y = get(celltype),fill=Structure))+
        geom_boxplot()+theme_bw()+
        scale_fill_manual(values = c(.cluster_cols,my36colors))+
        ylab(celltype)+
        theme(legend.title = element_text(size = 11*legend.size.factor),
              legend.text = element_text(size = 9*legend.size.factor),
              axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))+
        theme(legend.position = 'none')+
        guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))+
        stat_compare_means(method = "anova")
    }
    
  }else{
    group <- 'Location'
    metadata.use <- metadata[,c(group,colnames(metadata)[6:ncol(metadata)])]
    colnames(metadata.use)[1] <- 'Structure'
    data <- metadata.use[,c('Structure',celltype)]%>%as.data.frame()
    if (test=='Kruskal-Wallis') {
      plot <- ggplot(data,aes(x = Structure,y = get(celltype),fill=Structure))+
        geom_boxplot()+theme_bw()+
        scale_fill_manual(values = c(.cluster_cols,my36colors))+
        ylab(celltype)+
        theme(legend.title = element_text(size = 11*legend.size.factor),
              legend.text = element_text(size = 9*legend.size.factor))+
        theme(legend.position = 'none')+
        guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))+
        stat_compare_means(method = "kruskal.test")
    }else if (test=='Anova'){
      plot <- ggplot(data,aes(x = Structure,y = get(celltype),fill=Structure))+
        geom_boxplot()+theme_bw()+
        scale_fill_manual(values = c(.cluster_cols,my36colors))+
        ylab(celltype)+
        theme(legend.title = element_text(size = 11*legend.size.factor),
              legend.text = element_text(size = 9*legend.size.factor))+
        theme(legend.position = 'none')+
        guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))+
        stat_compare_means(method = "anova")
    }else if(test=='t.test'){
      plot <- ggplot(data,aes(x = Structure,y = get(celltype),fill=Structure))+
        geom_boxplot()+theme_bw()+
        scale_fill_manual(values = col_location)+
        ylab(celltype)+
        theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor))+
        guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))+
        theme(legend.position = 'none')+
        stat_compare_means(method='t.test',comparisons = list(c('Bdy','Mal'),c('Mal','nMal'),c('nMal','Bdy')))
    }else if(test=='Wilcoxon'){
      plot <- ggplot(data,aes(x = Structure,y = get(celltype),fill=Structure))+
        geom_boxplot()+theme_bw()+
        scale_fill_manual(values = col_location)+
        ylab(celltype)+
        theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor))+
        guides(fill=guide_legend(override.aes = list(size=0.5*legend.size.factor)))+
        theme(legend.position = 'none')+
        stat_compare_means(comparisons = list(c('Bdy','Mal'),c('Mal','nMal'),c('nMal','Bdy')))
    }
  }
  return(plot)
}

## ori
# data1 <- object[[group.by]]
# for (group in group.by) {
#   if (!is.factor(x = data[, group])) {
#     data[, group] <- factor(x = data[, group])
#   }
# }
# features <- colnames(x = data)
# colnames(x = data) <- features
# rownames(x = data) <- colnames(x = object)
# image.use <- object[[images[[image.idx]]]] #input
# coordinates <- GetTissueCoordinates(object = image.use)
# plot <- SingleSpatialPlot(data = cbind(coordinates,data[rownames(x = coordinates), features,drop = FALSE]),
#                           image = image.use, image.alpha = image.alpha, 
#                           col.by = features, cols = NULL, 
#                           alpha.by = NULL, 
#                           pt.alpha = alpha[1],
#                           geom = "spatial", cells.highlight = NULL, cols.highlight = NULL, 
#                           pt.size.factor = pt.size.factor, stroke = stroke, 
#                           crop = crop)+NoLegend()

## CellChatPlot ----
DoCellChatPlot <- function(source,target,weight,groupSize){
  
  
  if (source=='All'&target=='All') {
    plot <- netVisual_circle(weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge = FALSE,
                             title.name = "Interaction weights")
    
  }else if(target=='All'){
    plot <- netVisual_circle(weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge = FALSE,
                             title.name = "Interaction weights",sources.use = source)
  }else if(source=='All'){
    plot <- netVisual_circle(weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge = FALSE,
                             title.name = "Interaction weights",targets.use = target)
  }else{
    plot <- netVisual_circle(weight,vertex.weight = groupSize,weight.scale = TRUE,label.edge = FALSE,
                             title.name = "Interaction weights",targets.use = target,sources.use = source)
  }
  
  return(plot)
}

DoCellChatSignalingPlot <- function(cellchat,source,target,signaling){
  
  if (source=='All'&target=='All') {
    plot <- netVisual_bubble(cellchat,signaling = signaling)
    
    
  }else if(target=='All'){
    plot <- netVisual_bubble(cellchat,sources.use = source, signaling = signaling)
    
  }else if(source=='All'){
    plot <- netVisual_bubble(cellchat,targets.use = target, sources.use= c(setdiff(levels(cellchat@idents),target)),signaling = signaling)
    
  }else{
    plot <- netVisual_bubble(cellchat,targets.use = target, sources.use = source, signaling = signaling)
    
  }
  return(plot)
}
## ReconPlot ----
DoReconPlot <- function(meta_data,legend.size.factor,group_by){
  x.axis <- "UMAP_1"
  y.axis <- "UMAP_2"
  plot.data <- meta_data[,c(group_by,x.axis,y.axis)]
  ggplot(plot.data, aes_string(x = x.axis, y = y.axis)) +
    geom_point(aes_string(color = group_by), size = 0.8, alpha = 0.8) +
    labs(x = "UMAP_1", y = "UMAP_2") +
    theme_cowplot() +
    scale_color_manual(values = c(.cluster_cols,my36colors,c102)) +
    theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
          legend.key.size = unit(1*legend.size.factor, "lines"))+
    guides(color=guide_legend(override.aes = list(size=1.5*legend.size.factor))) -> p
  return(p)
}

## ReconGenePlot ----
DoReconGenePlot <- function(GenesInfo,expr_data,meta_data,legend.size.factor,facet,ncol){
  x.axis <- "UMAP_1"
  y.axis <- "UMAP_2"
  Expression_PlotData <- as.matrix(expr_data[GenesInfo$right_gene,rownames(meta_data)])
  if (GenesInfo$gene_number > 1 & facet=="Geometric mean") {
    Mean_expression <- apply(Expression_PlotData + 1, 2,
                             psych::geometric.mean) - 1
    plot.data <- data.frame(
      Expression = Mean_expression, meta_data,
      t(Expression_PlotData), check.names = FALSE
    )
    ggplot(plot.data %>% arrange(Expression), aes_string(x = x.axis, y = y.axis)) +
      geom_point(aes(color = Expression), alpha = 0.8, size = 0.8) +
      labs(x = "", y = "") +
      scale_colour_gradientn("", colors = SpatialColors(100)) +
      ggtitle('Mean')+
      theme_cowplot() +
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"))+
      guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor))) -> p
    
  }else if(GenesInfo$gene_number == 1){
    plot.data <- data.frame(
      Expression = Expression_PlotData, meta_data, 
      input_gene = Expression_PlotData)
    names(plot.data)[names(plot.data) == "input_gene"] = GenesInfo$right_gene
    plot.data <- melt(plot.data[,c(GenesInfo$right_gene, x.axis, y.axis)], id.vars = c(x.axis, y.axis))
    ggplot(plot.data %>% arrange(value), aes_string(x = x.axis, y = y.axis)) +
      geom_point(aes(color = value), alpha = 0.8, size = 0.8) +
      labs(x = "UMAP_1", y = "UMAP_2") +
      scale_colour_gradientn("", colors = SpatialColors(100)) +
      ggtitle(GenesInfo$right_gene)+
      theme_cowplot() +
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"))+
      guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor))) -> p
  }else if(GenesInfo$gene_number > 1 & facet=="Seperate"){
    plot.data <- data.frame(
      Expression = t(Expression_PlotData), meta_data)
    names(plot.data)[1:GenesInfo$gene_number] = GenesInfo$right_gene
    plot.data <- melt(plot.data[,c(GenesInfo$right_gene, x.axis, y.axis)], id.vars = c(x.axis, y.axis))
    ggplot(plot.data %>% arrange(value), aes_string(x = x.axis, y = y.axis)) +
      geom_point(aes(color = value), alpha = 0.8, size = 0.8) +
      labs(x = "UMAP_1", y = "UMAP_2") +
      facet_wrap(facets = 'variable',ncol = ncol)+
      scale_colour_gradientn("", colors = SpatialColors(100)) +
      theme_cowplot() +
      theme(legend.title = element_text(size = 11*legend.size.factor),legend.text = element_text(size = 9*legend.size.factor),
            legend.key.size = unit(1*legend.size.factor, "lines"),strip.text = element_text(size = 11*legend.size.factor))+
      guides(fill=guide_legend(override.aes = list(size=1.5*legend.size.factor))) -> p
  }else{
    print('error')
  }
  return(p)
}
