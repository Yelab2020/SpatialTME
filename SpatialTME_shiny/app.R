# Load packages----
library(shiny)
library(shiny.router)
library(shinydashboard)
library(shinyWidgets)
library(magrittr)
library(Seurat)
library(shinyjs)
library(shinyjqui)
library(figpatch)
library(CellChat)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
options(shiny.router.debug = T)
source('./UI.R')
source('./Plot_utility.R')
source('./buttonwithBusyIndicator.R')

# Close the window and shut down app----
# jscode <- "shinyjs.closeWindow = function() { window.close(); }"

# Load data ----
# load("./data/Initialize_expression.rda", envir = .GlobalEnv)   #count矩阵
load("./data/Saved_genes_panel.rda", envir = .GlobalEnv)    #不同cluster的ref-gene

# Both sample pages.
other_page <- div(
  titlePanel("Load")
)

# Creates router. We provide routing path and UI for this page.
router <- make_router(
  route("other", other_page)
)

# ui ----
ui <- fluidPage(
  tags$head(
    tags$style(
      HTML("
        .select-label {
          font-size: 20px; /* 设置字体大小 */
        }
      ")
    )
  ),
  useShinyjs(),
  uiOutput("tabs"),
  absolutePanel(
    id = "Load_button",
    # class = "btn btn-default",
    style = "top: 5px; right: 15px;",
    withBusyIndicatorUI(actionButton(inputId = "Load_dataset",label = "Load data",width = "100px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
  ),
  # router$ui,
  #actionButton("change_bdy", "Change query path-visium"),
  # 
  #actionButton("change_total", "Change query path-slide")

)

# server----
server <- function(input, output, session) {
  router$server(input, output, session)
  ## InputValue ----
  InputValue <- reactiveValues(
    dataset_choose=NULL,
    sample_choose=NULL,
    platform=NULL,
    geneinput = NULL,
    geneinfo = NULL,
    expr_data = NULL,
    image_data = NULL,
    meta_data = NULL,
    ## Structure
    PlotPar_facet = NULL,
    PlotPar_ncol = NULL,
    PlotPar_spotsize = NULL,
    BoundaryPlot_type=NULL,
    BoundaryPlot_spotsize=NULL,
    BoundaryPlot_alpha=NULL,
    BoundaryPlot_legendsize=NULL,
    ## DE
    DE_type=NULL,
    DE_DEG_structure=NULL,
    DE_GO_structure=NULL,
    DE_DEG_cluster=NULL,
    DE_GO_cluster=NULL,
    DE_DEG_pvalue=NULL,
    DE_DEG_FC=NULL,
    DE_DEG_N=NULL,
    DE_DEG_legendsize=NULL,
    DE_GO_pvalue=NULL,
    DE_GO_FC=NULL,
    DE_GO_N=NULL,
    DE_GO_legendsize=NULL,
    GOgene=NULL,
    DiffGenes_bdy=NULL,
    DiffGenes_cluster=NULL,
    GOres=NULL,
    ## Decon
    Decon_type=NULL,
    Decon_structure=NULL,
    DeconPie_spot=NULL,
    DeconBar_legendsize=NULL,
    CellLocation_structure=NULL,
    CellLocation_cell=NULL,
    CellLocation_spot=NULL,
    CellLocation_spotsize=NULL,
    CellLocationBox_legendsize=NULL,
    ## Cell location-slide
    CellLocation_spot_slide=NULL,
    CellLocation_cell_slide=NULL,
    CellLocation_spotsize_slide=NULL,
    CellLocationBox_legendsize_slide=NULL,
    ## Interaction
    CellChat_source=NULL,
    CellChat_target=NULL,
    CellLocation_signaling=NULL,
    cellchat_less=NULL,
    cellchat_df=NULL,
    recon_expr_data=NULL,
    recon_meta_data=NULL,
    Recon_geneinput = NULL,
    Recon_geneinfo = NULL,
    Recon_sample =NULL,
    Recon_groupby =NULL,
    Recon_legendsize = NULL,
    recon_expr_data_use=NULL,
    recon_meta_data_use=NULL,
    ReconGene_facet=NULL,
    ReconGene_ncol=NULL,
    CellLocation_test=NULL
  )
  
  ## load data ----
  observe({
    route_info <- get_query_param()
    selectedTab <- gsub("_"," ",route_info$selectedTab[[1]])
    platform <- route_info$platform[[1]]
    # 根据selectedTab的值设置选中的模块
    if (!is.null(selectedTab)) {
      updateTabsetPanel(session, "maintabset", selected = selectedTab)
    }
    # print(selectedTab)
    output$tabs <- renderUI({
      if(platform == "Visium") {
        tabsetPanel(
          id = "maintabset", type = "tabs",
          tabPanel("Spatial Structure", wellPanel(
            fluidRow(
              column(3, BoundaryInput()),
              column(9, BoundaryPlot())
            )
          )),
          tabPanel("Spatial Gene Expression", wellPanel(
            fluidRow(
              column(3, GeneInput()),
              column(9, GenePlot())
            )
          )),
          tabPanel("Spatial DE Analysis", wellPanel(
            fluidRow(
              column(3, DEInput()),
              column(9, DEPlot())
            )
          )),
          tabPanel("Spatial Cell Pattern", wellPanel(
            fluidRow(
              column(3, DeconInput()),
              column(5, DeconPlot()),
              column(4, DeconStatPlot())
            )
          )),
          tabPanel("Sub-spot GEP", wellPanel(
            fluidRow(
              column(3, ReconInput()),
              column(9, ReconPlot())
            )
          )),
          tabPanel("Cell Interaction", wellPanel(
            fluidRow(
              column(3, CellChatInput()),
              column(4, CellChatPlot()),
              column(5, CellChatSignalingPlot())
            )
          ))
        )
      } else {
        tabsetPanel(
          id = "maintabset", type = "tabs",
          tabPanel("Spatial Structure", wellPanel(
            fluidRow(
              column(3, BoundaryInput()),
              column(9, BoundaryPlot())
            )
          )),
          tabPanel("Spatial Gene Expression", wellPanel(
            fluidRow(
              column(3, GeneInput()),
              column(9, GenePlot())
            )
          )),
          tabPanel("Spatial DE Analysis", wellPanel(
            fluidRow(
              column(3, DEInput()),
              column(9, DEPlot())
            )
          )),
          tabPanel("Spatial Cell Pattern", wellPanel(
            fluidRow(
              column(3, DeconInput_slide()),
              column(5, DeconPlot_slide()),
              column(4, DeconStatPlot_slide())
            )
          )),
          tabPanel("Cell Interaction", wellPanel(
            fluidRow(
              column(3, CellChatInput()),
              column(4, CellChatPlot()),
              column(5, CellChatSignalingPlot())
            )
          ))
        )
      }
    })
    shinyjs::delay(700, {  # 延迟1秒
      shinyjs::hide("GeneExpr_download1")
      shinyjs::hide("Boundary_download1")
      shinyjs::hide("DeconPie_download1")
      shinyjs::hide("DeconPie_download2")
      shinyjs::hide("DeconBar_download1")
      shinyjs::hide("CellLocation_download1")
      shinyjs::hide("CellLocationBox_download1")
      shinyjs::hide("CellChat_download1")
      shinyjs::hide("CellChatSignaling_download1")
      shinyjs::hide("CellChatSignaling_download2")
      shinyjs::hide("DEValcano_download1")
      shinyjs::hide("DEValcano_download2")
      shinyjs::hide("DEGO_download1")
      shinyjs::hide("Recon_download1")
      shinyjs::hide("ReconGene_download1")
      shinyjs::hide("CellLocation_slide_download1")
      shinyjs::hide("CellLocation_slide_download2")
      shinyjs::hide("CellLocationBox_slide_download1")
    })
  })
  observeEvent(input$Load_dataset,{
      shinyjs::hide("GeneExpr_download1")
      shinyjs::hide("Boundary_download1")
      shinyjs::hide("DeconPie_download1")
      shinyjs::hide("DeconPie_download2")
      shinyjs::hide("DeconBar_download1")
      shinyjs::hide("CellLocation_download1")
      shinyjs::hide("CellLocationBox_download1")
      shinyjs::hide("CellChat_download1")
      shinyjs::hide("CellChatSignaling_download1")
      shinyjs::hide("CellChatSignaling_download2")
      shinyjs::hide("DEValcano_download1")
      shinyjs::hide("DEValcano_download2")
      shinyjs::hide("DEGO_download1")
      shinyjs::hide("Recon_download1")
      shinyjs::hide("ReconGene_download1")
      shinyjs::hide("CellLocation_slide_download1")
      shinyjs::hide("CellLocation_slide_download2")
      shinyjs::hide("CellLocationBox_slide_download1")
    withBusyIndicatorServer("Load_dataset", {
      route_info <- get_query_param()
      platform <- route_info$platform[[1]]
      dataset_choose <- route_info$dataset[[1]]
      sample_choose <- route_info$sample[[1]]
      InputValue$dataset_choose <- dataset_choose
      InputValue$sample_choose <- sample_choose
      InputValue$platform <- platform
      expr_data <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_expression.rds.gz"))
      InputValue$expr_data <- expr_data
      meta_data <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_metadata.rds.gz"))
      InputValue$meta_data <- meta_data
      cellchat_less <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_cellchatless.rds.gz"))
      InputValue$cellchat_less <- cellchat_less
      image_data <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_image.rds.gz"))
      InputValue$image_data <- image_data
      if (platform=='Visium') {
        recon_expr_data <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_Recon_expression.rds.gz"))
        recon_meta_data <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_Recon_metadata.rds.gz"))
        InputValue$recon_expr_data <- recon_expr_data
        InputValue$recon_meta_data <- recon_meta_data
        if ('total' %in% unique(meta_data$Location)) {
          DiffGenes_cluster <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_DiffGenes_cluster.rds.gz"))
          InputValue$DiffGenes_cluster <- DiffGenes_cluster
        }else{
          DiffGenes_bdy <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_DiffGenes_bdy.rds.gz"))
          DiffGenes_cluster <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_DiffGenes_cluster.rds.gz"))
          InputValue$DiffGenes_bdy <- DiffGenes_bdy
          InputValue$DiffGenes_cluster <- DiffGenes_cluster
        }
        updateSelectInput(session, "CellLocation_cell", choices = colnames(InputValue$meta_data[,6:ncol(InputValue$meta_data)]))
        if ('total' %in% unique(InputValue$meta_data$Location)) {
          updateSelectInput(session, "BoundaryPlot_type",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
          updateSelectInput(session, "DE_DEG_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
          updateSelectInput(session, "DE_GO_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
          updateSelectInput(session, "Decon_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
          updateSelectInput(session, "CellLocation_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
          updateSelectInput(session, "Recon_groupby",label = "Group by",choices = c("Cell type" = 'Subtypes',"Morphological adjusted clusters" = 'Morph_clusters','Morphological adjusted clusters_Cell type'='Morph_clusters_celltype'))
          updateSelectInput(session, inputId = "DE_GO_cluster",label = 'Location',choices = levels(InputValue$meta_data$Morph_clusters))
          updateSelectInput(session, inputId="DE_DEG_cluster",label = 'Location', choices = levels(InputValue$meta_data$Morph_clusters))
          updateSelectInput(session, inputId = "CellLocation_spot",label = 'Spots for plot',choices = c("All spots"='All',levels(InputValue$meta_data$Morph_clusters)))
          updateSelectInput(session, inputId = "DeconPie_spot",label = 'Spots for plot',choices = c("All spots"='All'))
        }else{
          updateSelectInput(session, "BoundaryPlot_type",label = "Structure type",choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis"),selected = "Mal-Bdy-nMal axis")
          updateSelectInput(session, "DE_DEG_structure",label = "Structure type",choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis"),selected = "Mal-Bdy-nMal axis")
          updateSelectInput(session, "DE_GO_structure",label = "Structure type",choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis"),selected = "Mal-Bdy-nMal axis")
          updateSelectInput(session, "Decon_structure",label = "Structure type",choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis"),selected = "Mal-Bdy-nMal axis")
          updateSelectInput(session, "CellLocation_structure",label = "Structure type",choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis"),selected = "Mal-Bdy-nMal axis")
          updateSelectInput(session, "Recon_groupby",label = "Group by",c("Cell type" = 'Subtypes',"Mal-Bdy-nMal Axis" = 'Location',"Morphological adjusted clusters" = 'Morph_clusters','Mal-Bdy-nMal Axis_Cell type'='Location_celltype','Morphological adjusted clusters_Cell type'='Morph_clusters_celltype'))
        }
      }else if(platform=='Slide_seq'){
        DiffGenes_cluster <- readRDS(paste0("./data/",platform,"/",dataset_choose,"/",sample_choose,"_DiffGenes_cluster.rds.gz"))
        InputValue$DiffGenes_cluster <- DiffGenes_cluster
        updateSelectInput(session, "CellLocation_cell_slide", choices = c('All cells'='All',unique(InputValue$meta_data$first_type)))
        updateSelectInput(session, "CellLocation_spot_slide", choices = c('All spots'='All',levels(InputValue$meta_data$Morph_clusters)))
        updateSelectInput(session, "BoundaryPlot_type",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
        updateSelectInput(session, "DE_DEG_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
        updateSelectInput(session, "DE_GO_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
        updateSelectInput(session, "Decon_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
        updateSelectInput(session, "CellLocation_structure",label = "Structure type",choices = c("Morphological adjusted cluster"),selected = "Morphological adjusted cluster")
        updateSelectInput(session, inputId = "DE_GO_cluster",label = 'Location',choices = levels(InputValue$meta_data$Morph_clusters))
        updateSelectInput(session, inputId="DE_DEG_cluster",label = 'Location', choices = levels(InputValue$meta_data$Morph_clusters))
        updateSelectInput(session, inputId = "CellLocation_spot",label = 'Spots for plot',choices = c("All spots"='All',levels(InputValue$meta_data$Morph_clusters)))
        updateSelectInput(session, inputId = "DeconPie_spot",label = 'Spots for plot',choices = c("All spots"='All'))
      }
      updateSelectInput(session, "CellChat_source", choices = c("","All",levels(cellchat_less@idents)))
      updateSelectInput(session, "CellChat_target", choices = c("","All",levels(cellchat_less@idents)))
      updateSelectInput(session, "CellChat_signaling", choices = 'choose source and target first',selected = NULL)
    })
  })
  
  ## Spatial Structure ----
  observeEvent(input$Boundary_submit,{
    withBusyIndicatorServer("Boundary_submit", {
      InputValue$BoundaryPlot_spotsize <- input$BoundaryPlot_spotsize
      InputValue$BoundaryPlot_alpha <- input$BoundaryPlot_alpha
      InputValue$BoundaryPlot_legendsize <- input$BoundaryPlot_legendsize
      InputValue$BoundaryPlot_type <- input$BoundaryPlot_type
      if (!is.null(InputValue$expr_data)&&!is.null(InputValue$image_data)) {
        output$Boundary_plot <- renderPlot({
          DoBoundaryPlot(
            structure.type=InputValue$BoundaryPlot_type,
            metadata=InputValue$meta_data,
            image.use=InputValue$image_data,
            pt.size.factor = InputValue$BoundaryPlot_spotsize,
            alpha = InputValue$BoundaryPlot_alpha,
            legend.size.factor=InputValue$BoundaryPlot_legendsize
          )
        })
        shinyjs::show("Boundary_download1", anim = TRUE, animType = "fade")
      }else{
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
      
    })
  })
  ## download 
  output$Boundary_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_BoundaryPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$Boundary_plot_size$width / 100,
          height = input$Boundary_plot_size$height / 100 )
      plot(
        DoBoundaryPlot(
          structure.type=InputValue$BoundaryPlot_type,
          metadata=InputValue$meta_data,
          image.use=InputValue$image_data,
          pt.size.factor = InputValue$BoundaryPlot_spotsize,
          alpha = InputValue$BoundaryPlot_alpha,
          legend.size.factor=InputValue$BoundaryPlot_legendsize
        )
      )
      dev.off()
    }
  )
  
  ## gene expr plot ----
  observe({
    updateTextAreaInput(
      session = session,
      inputId = "GeneInput_text",
      value = input$GeneInput_saved
    )
  })
  
  observeEvent(input$GeneInput_text, {
    # 获取输入基因的数量
    geneinput <- toupper(input$GeneInput_text)
    geneinfo <- checkGeneList(isolate(geneinput),InputValue$expr_data)
    if (geneinfo$gene_number >= 1) {
      updateSelectInput(session, "PlotPar_facet", selected = "Geometric mean")
    }else{
      updateSelectInput(session, "PlotPar_facet", selected = "Seperate")
    }
  })
  
  observeEvent(input$GeneExpr_clear, {
    updateTextAreaInput(session, "GeneInput_text", value = "")
    updateTextAreaInput(session, "GeneInput_saved", value = "")
  })
  
  observeEvent(input$GeneExpr_submit,{
    withBusyIndicatorServer("GeneExpr_submit", {
      InputValue$PlotPar_facet <- input$PlotPar_facet
      InputValue$PlotPar_ncol <- input$PlotPar_ncol
      InputValue$PlotPar_spotsize <- input$PlotPar_spotsize
      ## check load ----
      if (!is.null(InputValue$expr_data)&&!is.null(InputValue$image_data)) {
        ## check geneinput ----
        InputValue$geneinput <- toupper(input$GeneInput_text)
        InputValue$geneinfo <- checkGeneList(isolate(InputValue$geneinput),InputValue$expr_data)
        
        if (InputValue$geneinfo$gene_number >= 1) {
          if (!InputValue$geneinfo$all_genes_avaliable_flag) {
            message = paste(
              paste(InputValue$geneinfo$wrong_gene, collapse = ","),
              "not found in the dataset!"
            )
            sendSweetAlert(session = session,title = message,type = "warning")
          }
          if (InputValue$geneinfo$gene_number>10&InputValue$PlotPar_facet=='Seperate') {
            message = paste(
              'Please do not draw more than 10 genes separately!'
            )
            sendSweetAlert(session = session,title = message,type = "warning")
          }else{
            output$GeneExpr_plot <- renderPlot({
              DoGeneExprPlot(
                gene=InputValue$geneinfo$right_gene,
                expr_data=InputValue$expr_data,
                image.use=InputValue$image_data,
                seperate.plot = InputValue$PlotPar_facet,
                ncol = InputValue$PlotPar_ncol,
                pt.size.factor = InputValue$PlotPar_spotsize
              )
            })
            shinyjs::show("GeneExpr_download1", anim = TRUE, animType = "fade")
          }
        } else{
          InputValue$geneinfo <- NULL
          message = ("No input genes could be used!")
          sendSweetAlert(session = session,title = message,type = "error")
        }
        
      }else{ ## unload ----
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  ## download 
  output$GeneExpr_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_EmbeddingPlot.Gene",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$GeneExpr_plot_size$width / 100,
          height = input$GeneExpr_plot_size$height / 100 )
      plot(
        DoGeneExprPlot(
          gene=InputValue$geneinfo$right_gene,
          expr_data=InputValue$expr_data,
          image.use=InputValue$image_data,
          seperate.plot = InputValue$PlotPar_facet,
          ncol = InputValue$PlotPar_ncol,
          pt.size.factor = InputValue$PlotPar_spotsize
        )
      )
      dev.off()
    }
  )
  ## DEG ----
  ##> update cluster  ----
  observeEvent(input$DE_GO_structure, {
    req(input$DE_GO_structure)
    if (!is.null(InputValue$meta_data)) {
      if (input$DE_GO_structure=='Morphological adjusted cluster') {
        updateSelectInput(session, inputId = "DE_GO_cluster",label = 'Location',choices = levels(InputValue$meta_data$Morph_clusters))
      }else{
        updateSelectInput(session, inputId="DE_GO_cluster",label = 'Location', choices = c("Malignant spots" = 'Mal', "Boundary spots" = 'Bdy', "non-Malignant spots" = 'nMal'))
      }
    }else{
      updateSelectInput(session, inputId="DE_GO_cluster",label = 'Location', choices = c('Load data firstly'))
    }
  })
  observeEvent(input$DE_DEG_structure, {
    req(input$DE_DEG_structure)
    
    if (!is.null(InputValue$meta_data)) {
      if (input$DE_DEG_structure=='Morphological adjusted cluster') {
        updateSelectInput(session, inputId="DE_DEG_cluster",label = 'Location', choices = levels(InputValue$meta_data$Morph_clusters))
      }else{
        updateSelectInput(session, inputId="DE_DEG_cluster",label = 'Location', choices = c("Malignant spots" = 'Mal', "Boundary spots" = 'Bdy', "non-Malignant spots" = 'nMal'))
      }
    }else{
      updateSelectInput(session, inputId="DE_DEG_cluster",label = 'Location', choices = c('Load data firstly'))
    }
  })
  ##> valcano plot ----
  observeEvent(input$DE_DEGsubmit,{
    withBusyIndicatorServer("DE_DEGsubmit", {
      InputValue$DE_DEG_structure <- input$DE_DEG_structure
      InputValue$DE_DEG_cluster <- input$DE_DEG_cluster
      InputValue$DE_DEG_pvalue <- input$DE_DEG_pvalue
      InputValue$DE_DEG_FC <- input$DE_DEG_FC
      InputValue$DE_DEG_N <- input$DE_DEG_N
      InputValue$DE_DEG_legendsize <- input$DE_DEG_legendsize
      
      if (!is.null(InputValue$expr_data)&&!is.null(InputValue$image_data)) {
        if (InputValue$DE_DEG_structure=='Mal-Bdy-nMal axis') {
          output$DEValcano_plot <- renderPlot({
            DoDEValcanoPlot(
              DiffGenes=InputValue$DiffGenes_bdy,
              cluster = InputValue$DE_DEG_cluster,
              cut_off_pvalue = InputValue$DE_DEG_pvalue,
              cut_off_logFC = InputValue$DE_DEG_FC,
              n=InputValue$DE_DEG_N,
              legend.size.factor=InputValue$DE_DEG_legendsize
            )
          })
          shinyjs::show("DEValcano_download1", anim = TRUE, animType = "fade")
          shinyjs::show("DEValcano_download2", anim = TRUE, animType = "fade")
        }else{
          output$DEValcano_plot <- renderPlot({
            DoDEValcanoPlot(
              DiffGenes=InputValue$DiffGenes_cluster,
              cluster = InputValue$DE_DEG_cluster,
              cut_off_pvalue = InputValue$DE_DEG_pvalue,
              cut_off_logFC = InputValue$DE_DEG_FC,
              n=InputValue$DE_DEG_N,
              legend.size.factor=InputValue$DE_DEG_legendsize
            )
          })
          shinyjs::show("DEValcano_download1", anim = TRUE, animType = "fade")
          shinyjs::show("DEValcano_download2", anim = TRUE, animType = "fade")
        }
      }else{
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  ##> GO plot ----
  ### analysis
  observeEvent(input$DE_GOanalysis,{
    withBusyIndicatorServer("DE_GOanalysis", {
      InputValue$DE_GO_structure <- input$DE_GO_structure
      InputValue$DE_GO_cluster <- input$DE_GO_cluster
      InputValue$DE_GO_pvalue <- input$DE_GO_pvalue
      InputValue$DE_GO_FC <- input$DE_GO_FC
      
      if (!is.null(InputValue$expr_data)&&!is.null(InputValue$meta_data)) {
        if (InputValue$DE_GO_structure=='Mal-Bdy-nMal axis') {
          InputValue$GOgene <- DoGOGene(
            DiffGenes=InputValue$DiffGenes_bdy,
            cluster=InputValue$DE_GO_cluster,
            cut_off_pvalue = InputValue$DE_GO_pvalue,
            cut_off_logFC = InputValue$DE_GO_FC
          )
          if (length(InputValue$GOgene)==0) {
            message = paste("No genes input!")
            sendSweetAlert(session = session,title = message,type = "warning")
          }else{
            InputValue$GOres <- DoGOres(DiffGenes=InputValue$GOgene)
          }
          
        }else{
          InputValue$GOgene <- DoGOGene(
            DiffGenes=InputValue$DiffGenes_cluster,
            cluster=InputValue$DE_GO_cluster,
            cut_off_pvalue = InputValue$DE_GO_pvalue,
            cut_off_logFC = InputValue$DE_GO_FC
          )
          if (length(InputValue$GOgene)==0) {
            message = paste("No genes input!")
            sendSweetAlert(session = session,title = message,type = "warning")
          }else{
            InputValue$GOres <- DoGOres(DiffGenes=InputValue$GOgene)
          }
        }
      }else{
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  ## plot
  observeEvent(input$DE_GOplot,{
    withBusyIndicatorServer("DE_GOplot", {
      InputValue$DE_GO_N <- input$DE_GO_N
      InputValue$DE_GO_legendsize <- input$DE_GO_legendsize
      
      if (!is.null(InputValue$GOres)) {
        print('yes')
        output$DEGO_plot <- renderPlot({
          DoGOPlot(
            GO = InputValue$GOres,
            n=InputValue$DE_GO_N,
            legend.size.factor=InputValue$DE_GO_legendsize
          )
        })
        shinyjs::show("DEGO_download1", anim = TRUE, animType = "fade")
      }else{
        message = paste("Please analyze firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  ##> download ----
  output$DEValcano_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"-",InputValue$DE_cluster,"_DEValcanoPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      if (InputValue$DE_DEG_structure=='Mal-Bdy-nMal axis'){
        pdf(file,
            width = input$DEValcano_plot_size$width / 100,
            height = input$DEValcano_plot_size$height / 100 )
        plot(
          DoDEValcanoPlot(
            DiffGenes=InputValue$DiffGenes_bdy,
            cluster = InputValue$DE_DEG_cluster,
            cut_off_pvalue = InputValue$DE_DEG_pvalue,
            cut_off_logFC = InputValue$DE_DEG_FC,
            n=InputValue$DE_DEG_N,
            legend.size.factor=InputValue$DE_DEG_legendsize
          )
        )
        dev.off()
      }else{
        pdf(file,
            width = input$DEValcano_plot_size$width / 100,
            height = input$DEValcano_plot_size$height / 100 )
        plot(
          DoDEValcanoPlot(
            DiffGenes=InputValue$DiffGenes_cluster,
            cluster = InputValue$DE_DEG_cluster,
            cut_off_pvalue = InputValue$DE_DEG_pvalue,
            cut_off_logFC = InputValue$DE_DEG_FC,
            n=InputValue$DE_DEG_N,
            legend.size.factor=InputValue$DE_DEG_legendsize
          )
        )
        dev.off()
      }
    }
  )
  output$DEValcano_download2 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"-",InputValue$DE_cluster,"_DEG",stringi::stri_rand_strings(1, 10),".csv")},
    content = function(file){
      set.seed(1234)
      if (InputValue$DE_DEG_structure=='Mal-Bdy-nMal axis'){
        data <- InputValue$DiffGenes_bdy[[InputValue$DE_DEG_cluster]]
        write.csv(data, file, col.names = T, row.names = F, quote = F)
      }else{
        data <- InputValue$DiffGenes_cluster[[InputValue$DE_DEG_cluster]]
        write.csv(data, file, col.names = T, row.names = F, quote = F)
      }
    }
  )
  output$DEGO_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"-",InputValue$DE_cluster,"_DEGOPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$DEGO_plot_size$width / 100,
          height = input$DEGO_plot_size$height / 100 )
      plot(
        DoGOPlot(
          GO = InputValue$GOres,
          n=InputValue$DE_GO_N,
          legend.size.factor=InputValue$DE_GO_legendsize
        )
      )
      dev.off()
    }
  )

  ## Cell Distribution ----
  ##> update cluster  ----
  observeEvent(input$CellLocation_structure, {
    req(input$CellLocation_structure)
    if (!is.null(InputValue$meta_data)) {
      if (input$CellLocation_structure=='Morphological adjusted cluster') {
        updateSelectInput(session, inputId = "CellLocation_spot",label = 'Spots for plot',choices = c("All spots"='All',levels(InputValue$meta_data$Morph_clusters)))
        updateSelectInput(session, inputId = "CellLocation_test",label = 'Statistical analysis',choices = c('Kruskal-Wallis','Anova'))
      }else{
        updateSelectInput(session, inputId="CellLocation_spot",label = 'Spots for plot', choices = c("All spots"='All',"Malignant spots" = 'Mal', "Boundary spots" = 'Bdy', "non-Malignant spots" = 'nMal'))
        updateSelectInput(session, inputId = "CellLocation_test",label = 'Statistical analysis',choices = c("t.test","Wilcoxon",'Kruskal-Wallis','Anova'))
      }
    }else{
      updateSelectInput(session, inputId="CellLocation_spot",label = 'Spots for plot', choices = c('Load data firstly'))
      updateSelectInput(session, inputId = "CellLocation_test",label = 'Statistical analysis',choices = c('Kruskal-Wallis','Anova'))
    }
  })
  observeEvent(input$Decon_structure, {
    req(input$Decon_structure)
    if (!is.null(InputValue$meta_data)) {
      if (input$Decon_structure=='Morphological adjusted cluster') {
        updateSelectInput(session, inputId = "DeconPie_spot",label = 'Spots for plot',choices = c("All spots" = 'All'))
      }else{
        updateSelectInput(session, inputId="DeconPie_spot",label = 'Spots for plot', choices = c("All spots" = 'All',"Malignant spots" = 'Mal', "Boundary spots" = 'Bdy', "non-Malignant spots" = 'nMal'))
      }
    }else{
      updateSelectInput(session, inputId="DeconPie_spot",label = 'Spots for plot', choices = c('Load data firstly'))
    }
  })
  
  ##> Decon plot ----
  observeEvent(input$Decon_submit,{
    withBusyIndicatorServer("Decon_submit", {
      InputValue$DeconPie_spot <- input$DeconPie_spot
      InputValue$DeconBar_legendsize <- input$DeconBar_legendsize
      InputValue$Decon_structure <- input$Decon_structure
      
      if (InputValue$platform=='Visium') {
        if (!is.null(InputValue$meta_data)) {
          output$DeconPie_plot <- renderImage({
            filename <- normalizePath(file.path('./data',
                                                paste0(InputValue$platform,"/",InputValue$dataset_choose,"/",InputValue$sample_choose, '_DeconPieplot_',InputValue$DeconPie_spot,'.png')))
            list(src = filename)
          }, deleteFile = FALSE)
          output$DeconBar_plot <- renderPlot({
            DoDeconBarPlot(
              structure=InputValue$Decon_structure,
              metadata=InputValue$meta_data,
              legend.size.factor=InputValue$DeconBar_legendsize
            )
          })
          shinyjs::show("DeconPie_download1", anim = TRUE, animType = "fade")
          shinyjs::show("DeconPie_download2", anim = TRUE, animType = "fade")
          shinyjs::show("DeconBar_download1", anim = TRUE, animType = "fade")
        }else{
          message = paste("Please load data firstly!")
          sendSweetAlert(session = session,title = message,type = "warning")
        }
      }else{
        if (!is.null(InputValue$meta_data)) {
          output$DeconPie_plot <- renderPlot({
            DoDeconPlot_slide(
              metadata=InputValue$meta_data,
              image.use=InputValue$image_data,
              legend.size.factor=InputValue$DeconBar_legendsize
            )
          })
          output$DeconBar_plot <- renderPlot({
            DoDeconBarPlot_slide(
              metadata=InputValue$meta_data,
              legend.size.factor=InputValue$DeconBar_legendsize
            )
          })
          shinyjs::show("DeconPie_download1", anim = TRUE, animType = "fade")
          shinyjs::show("DeconPie_download2", anim = TRUE, animType = "fade")
          shinyjs::show("DeconBar_download1", anim = TRUE, animType = "fade")
        }else{
          message = paste("Please load data firstly!")
          sendSweetAlert(session = session,title = message,type = "warning")
        }
      }
    })
  })
  ## download 
  output$DeconBar_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_DeconBarPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
        pdf(file,
            width = input$DeconBar_plot_size$width / 100,
            height = input$DeconBar_plot_size$height / 100 )
        plot(
          DoDeconBarPlot(
            structure=InputValue$Decon_structure,
            metadata=InputValue$meta_data,
            legend.size.factor=InputValue$DeconBar_legendsize
          )
        )
        dev.off()
    }
  )

  output$DeconPie_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_DeconPiePlot",stringi::stri_rand_strings(1, 10),".png")},
    content = function(file) {
      set.seed(1234)
        filename <- normalizePath(file.path('./data',
                                            paste0(InputValue$platform,"/",InputValue$dataset_choose,"/",InputValue$sample_choose, '_DeconPieplot_',InputValue$DeconPie_spot,'.png')))
        p <- fig(filename)
        png(file,width = 500,height = 500,units="px",res=500)
        print(p)
        dev.off()
    }
  )
  
  output$DeconPie_download2 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_DeconData",stringi::stri_rand_strings(1, 10),".csv")},
    content = function(file) {
      set.seed(1234)
        data <- InputValue$meta_data[,4:ncol(InputValue$meta_data)]
        write.csv(data, file, col.names = T, row.names = T, quote = F)

    }
  )
  ##> Cell location plot ----
  observeEvent(input$CellLocation_submit,{
    withBusyIndicatorServer("CellLocation_submit", {
      InputValue$CellLocation_structure <- input$CellLocation_structure
      InputValue$CellLocation_cell <- input$CellLocation_cell
      InputValue$CellLocation_spot <- input$CellLocation_spot
      InputValue$CellLocation_spotsize <- input$CellLocation_spotsize
      InputValue$CellLocationBox_legendsize <- input$CellLocationBox_legendsize
      InputValue$CellLocation_test <- input$CellLocation_test
      
      if (!is.null(InputValue$expr_data)&&!is.null(InputValue$image_data)) {
          output$CellLocation_plot <- renderPlot({
            DoCellLocationPlot(
              structure=InputValue$CellLocation_structure,
              metadata=InputValue$meta_data,
              image.use=InputValue$image_data,
              celltype=InputValue$CellLocation_cell,
              location=InputValue$CellLocation_spot,
              pt.size.factor = InputValue$CellLocation_spotsize
            )
          })
          output$CellLocationBox_plot <- renderPlot({
            DoCellLocationBoxPlot(
              structure=InputValue$CellLocation_structure,
              metadata=InputValue$meta_data,
              celltype=InputValue$CellLocation_cell,
              legend.size.factor=InputValue$CellLocationBox_legendsize,
              test=InputValue$CellLocation_test
            )
          })
          shinyjs::show("CellLocation_download1", anim = TRUE, animType = "fade")
          shinyjs::show("CellLocationBox_download1", anim = TRUE, animType = "fade")
      }else{
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  ## download 
  output$CellLocation_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_",InputValue$CellLocation_cell,"_",InputValue$CellLocation_spot,"_CellLocationPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$CellLocation_plot_size$width / 100,
          height = input$CellLocation_plot_size$height / 100 )
      plot(
        DoCellLocationPlot(
          structure=InputValue$CellLocation_structure,
          metadata=InputValue$meta_data,
          image.use=InputValue$image_data,
          celltype=InputValue$CellLocation_cell,
          location=InputValue$CellLocation_spot,
          pt.size.factor = InputValue$CellLocation_spotsize
        )
      )
      dev.off()
    }
  )
  output$CellLocationBox_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_",InputValue$CellLocation_cell,"_CellLocationBoxPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$CellLocationBox_plot_size$width / 100,
          height = input$CellLocationBox_plot_size$height / 100 )
      plot(
        DoCellLocationBoxPlot(
          structure=InputValue$CellLocation_structure,
          metadata=InputValue$meta_data,
          celltype=InputValue$CellLocation_cell,
          legend.size.factor=InputValue$CellLocationBox_legendsize,
          test=InputValue$CellLocation_test
        )
      )
      dev.off()
    }
  )
  ##> Cell location-slide----
  observeEvent(input$CellLocation_submit_slide,{
    withBusyIndicatorServer("CellLocation_submit_slide", {
      InputValue$CellLocation_cell_slide <- input$CellLocation_cell_slide
      InputValue$CellLocation_spot_slide <- input$CellLocation_spot_slide
      InputValue$CellLocation_spotsize_slide <- input$CellLocation_spotsize_slide
      InputValue$CellLocationBox_legendsize_slide <- input$CellLocationBox_legendsize_slide
      
      if (!is.null(InputValue$expr_data)&&!is.null(InputValue$image_data)) {
        output$CellLocation_plot_slide <- renderPlot({
          DoCellLocationPlot_slide(
            metadata=InputValue$meta_data,
            image.use=InputValue$image_data,
            celltype=InputValue$CellLocation_cell_slide,
            location=InputValue$CellLocation_spot_slide,
            pt.size.factor = InputValue$CellLocation_spotsize_slide,
            legend.size.factor=InputValue$CellLocationBox_legendsize_slide
          )
        })
        output$CellLocationBox_plot_slide <- renderPlot({
          DoDeconBarPlot_slide(
            metadata=InputValue$meta_data,
            legend.size.factor=InputValue$CellLocationBox_legendsize_slide
          )
        })
        shinyjs::show("CellLocation_slide_download1", anim = TRUE, animType = "fade")
        shinyjs::show("CellLocation_slide_download2", anim = TRUE, animType = "fade")
        shinyjs::show("CellLocationBox_slide_download1", anim = TRUE, animType = "fade")
      }else{
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  ## download 
  output$CellLocation_slide_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_",InputValue$CellLocation_cell,"_",InputValue$CellLocation_spot,"_CellLocationPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$CellLocation_plot_slide_size$width / 100,
          height = input$CellLocation_plot_slide_size$height / 100 )
      plot(
        DoCellLocationPlot_slide(
          metadata=InputValue$meta_data,
          image.use=InputValue$image_data,
          celltype=InputValue$CellLocation_cell_slide,
          location=InputValue$CellLocation_spot_slide,
          pt.size.factor = InputValue$CellLocation_spotsize_slide
        )
      )
      dev.off()
    }
  )
  output$CellLocationBox_slide_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_",InputValue$CellLocation_cell,"_CellLocationBoxPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$CellLocationBox_legendsize_slide_size$width / 100,
          height = input$CellLocationBox_legendsize_slide_size$height / 100 )
      plot(
        DoDeconBarPlot_slide(
          metadata=InputValue$meta_data,
          legend.size.factor=InputValue$CellLocationBox_legendsize_slide
        )
      )
      dev.off()
    }
  )
  output$CellLocation_slide_download2 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_CellTypeAnno",stringi::stri_rand_strings(1, 10),".csv")},
    content = function(file) {
      set.seed(1234)
      data <- InputValue$meta_data[,c(5,16)]
      write.csv(data, file, col.names = T, row.names = T, quote = F)
    }
  )
  ## Recon plot ----
  observeEvent(input$Recon_submit,{
    withBusyIndicatorServer("Recon_submit", {
      # InputValue$Recon_sample <- input$Recon_sample
      InputValue$Recon_groupby <- input$Recon_groupby
      InputValue$Recon_legendsize <- input$Recon_legendsize
      ## check load ----
      if (!is.null(InputValue$recon_expr_data)&&!is.null(InputValue$recon_meta_data)) {
        InputValue$recon_meta_data_use = InputValue$recon_meta_data 
        InputValue$recon_expr_data_use = InputValue$recon_expr_data
        ## dimplot ----
        output$Recon_plot <- renderPlot({
          DoReconPlot(
            meta_data=InputValue$recon_meta_data_use,
            legend.size.factor = InputValue$Recon_legendsize,
            group_by=InputValue$Recon_groupby
          )
        })
        shinyjs::show("Recon_download1", anim = TRUE, animType = "fade")
      }else{ ## unload ----
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  
  ## ReconGene plot ----
  observeEvent(input$ReconGene_submit,{
    withBusyIndicatorServer("ReconGene_submit", {
      InputValue$ReconGene_legendsize <- input$ReconGene_legendsize
      InputValue$ReconGene_facet=input$ReconGene_facet
      InputValue$ReconGene_ncol=input$ReconGene_ncol
      ## check load ----
      if (!is.null(InputValue$recon_expr_data)&&!is.null(InputValue$recon_meta_data)) {
        InputValue$recon_meta_data_use = InputValue$recon_meta_data 
        InputValue$recon_expr_data_use = InputValue$recon_expr_data
        ## check geneinput ----
        InputValue$Recon_geneinput <- toupper(input$Recon_gene)
        InputValue$Recon_geneinfo <- checkGeneList(isolate(InputValue$Recon_geneinput),InputValue$recon_expr_data_use)
        
        if (InputValue$Recon_geneinfo$gene_number >= 1) {
          if (!InputValue$Recon_geneinfo$all_genes_avaliable_flag) {
            message = paste(
              paste(InputValue$Recon_geneinfo$wrong_gene, collapse = ","),
              "not found in the dataset!"
            )
            sendSweetAlert(session = session,title = message,type = "warning")
          }
          print(InputValue$Recon_geneinfo)
          print(InputValue$ReconGene_legendsize)
          print(InputValue$ReconGene_facet)
          print(InputValue$ReconGene_ncol)
          
          output$ReconGene_plot <- renderPlot({
            DoReconGenePlot(
              GenesInfo=InputValue$Recon_geneinfo,
              expr_data=InputValue$recon_expr_data,
              meta_data=InputValue$recon_meta_data,
              legend.size.factor = InputValue$ReconGene_legendsize,
              facet=InputValue$ReconGene_facet,
              ncol=InputValue$ReconGene_ncol
            )
          })
          shinyjs::show("ReconGene_download1", anim = TRUE, animType = "fade")
          
        } else{
          message = ("No input genes could be used!")
          sendSweetAlert(session = session,title = message,type = "error")
        }
        
      }else{ ## unload ----
        message = paste("Please load data firstly!")
        sendSweetAlert(session = session,title = message,type = "warning")
      }
    })
  })
  
  ## download ----
  output$Recon_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_",InputValue$Recon_groupby,"_ReconPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$Recon_plot_size$width / 100,
          height = input$Recon_plot_size$height / 100 )
      plot(
        DoReconPlot(
          meta_data=InputValue$recon_meta_data_use,
          legend.size.factor = InputValue$Recon_legendsize,
          group_by=InputValue$Recon_groupby
        )
      )
      dev.off()
    }
  )
  output$ReconGene_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_",InputValue$Recon_gene,"_FeaturePlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$ReconGene_plot_size$width / 100,
          height = input$ReconGene_plot_size$height / 100 )
      plot(
        output$ReconGene_plot <- renderPlot({
          DoReconGenePlot(
            GenesInfo=InputValue$Recon_geneinfo,
            expr_data=InputValue$recon_expr_data_use,
            meta_data=InputValue$recon_meta_data_use,
            legend.size.factor = InputValue$ReconGene_legendsize,
            facet=InputValue$ReconGene_facet,
            ncol=InputValue$ReconGene_ncol
          )
        })
      )
      dev.off()
    }
  )
  
  
  ## Cell chat plot ----
  ## update signaling value 
  observeEvent(c(input$CellChat_source, input$CellChat_target), {
    req(input$CellChat_source, input$CellChat_target)
    if (!is.null(input$CellChat_source) && !is.null(input$CellChat_target)) {
      cellchat <- InputValue$cellchat_less
      if(input$CellChat_target=='All' && input$CellChat_source=='All'){
        updateSelectInput(session, "CellChat_signaling", choices = cellchat@netP$pathways)
      }else if(input$CellChat_target=='All'){
        cellchat_df <-  subsetCommunication(cellchat,slot.name = 'net', sources.use = input$CellChat_source)
        updateSelectInput(session, "CellChat_signaling", choices = unique(cellchat_df$pathway_name))
        
      }else if(input$CellChat_source=='All'){
        cellchat_df <-  subsetCommunication(cellchat,slot.name = 'net', targets.use = input$CellChat_target)
        updateSelectInput(session, "CellChat_signaling", choices = unique(cellchat_df$pathway_name))
      }else{
        cellchat_df <-  subsetCommunication(cellchat,slot.name = 'net', sources.use = input$CellChat_source,targets.use = input$CellChat_target)
        updateSelectInput(session, "CellChat_signaling", choices = unique(cellchat_df$pathway_name))
      }
    }else{
      updateSelectInput(session, "CellChat_signaling", choices = 'choose source and target first',selected = NULL)
    }
  })
  
  ## plot 
  observeEvent(input$CellChat_submit,{
    InputValue$CellChat_signaling <- input$CellChat_signaling
    InputValue$CellChat_source <- input$CellChat_source
    InputValue$CellChat_target <- input$CellChat_target
    
    if (!is.null(InputValue$cellchat_less)&&!is.null(InputValue$CellChat_source)&&!is.null(InputValue$CellChat_target)&&!is.null(InputValue$CellChat_signaling)) {
      cellchat <- InputValue$cellchat_less
      output$CellChat_plot <- renderPlot({
        DoCellChatPlot(source=InputValue$CellChat_source,target=InputValue$CellChat_target,weight=cellchat@net$weight,groupSize=as.numeric(table(cellchat@idents)))
      })
      output$CellChatSignaling_plot <- renderPlot({ 
        DoCellChatSignalingPlot(cellchat=cellchat,source=InputValue$CellChat_source,target=InputValue$CellChat_target,signaling=InputValue$CellChat_signaling)
      })
      shinyjs::show("CellChat_download1", anim = TRUE, animType = "fade")
      shinyjs::show("CellChatSignaling_download1", anim = TRUE, animType = "fade")
      shinyjs::show("CellChatSignaling_download2", anim = TRUE, animType = "fade")
      
      withBusyIndicatorServer("CellChat_submit", {
        Sys.sleep(1)
      })
      
    }else if(is.null(InputValue$cellchat_less)){
      message = paste("Please load data firstly!")
      sendSweetAlert(session = session,title = message,type = "warning")
    }else if(is.null(InputValue$CellChat_target)){
      message = paste("Please choose target!")
      sendSweetAlert(session = session,title = message,type = "warning")
    }else if (is.null(InputValue$CellChat_source)){
      message = paste("Please choose source")
      sendSweetAlert(session = session,title = message,type = "warning")
    }else{
      message = paste("Please choose signaling!")
      sendSweetAlert(session = session,title = message,type = "warning")
    }
    
  })
  ## download 
  output$CellChat_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_CellChatPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      cellchat <- InputValue$cellchat_less
      if (InputValue$CellChat_source=='All'&InputValue$CellChat_target=='All') {
        pdf(file,
            width = input$CellChat_plot_size$width / 100,
            height = input$CellChat_plot_size$height / 100 )
        netVisual_circle(cellchat@net$weight,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale = TRUE,label.edge = FALSE,
                                 title.name = "Interaction weights")
        dev.off()
        
      }else if(InputValue$CellChat_target=='All'){
        pdf(file,
            width = input$CellChat_plot_size$width / 100,
            height = input$CellChat_plot_size$height / 100 )
        plot <- netVisual_circle(cellchat@net$weight,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale = TRUE,label.edge = FALSE,
                                 title.name = "Interaction weights",sources.use = InputValue$CellChat_source)
        dev.off()
        
        
      }else if(InputValue$CellChat_source=='All'){
        pdf(file,
            width = input$CellChat_plot_size$width / 100,
            height = input$CellChat_plot_size$height / 100 )
        plot <- netVisual_circle(cellchat@net$weight,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale = TRUE,label.edge = FALSE,
                                 title.name = "Interaction weights",targets.use = InputValue$CellChat_target)
        dev.off()
        
      }else{
        pdf(file,
            width = input$CellChat_plot_size$width / 100,
            height = input$CellChat_plot_size$height / 100 )
        plot <- netVisual_circle(cellchat@net$weight,vertex.weight = as.numeric(table(cellchat@idents)),weight.scale = TRUE,label.edge = FALSE,
                                 title.name = "Interaction weights",targets.use = InputValue$CellChat_target,sources.use = InputValue$CellChat_source)
        dev.off()
        
      }
    }
  )
  output$CellChatSignaling_download1 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_CellChatSignalPlot",stringi::stri_rand_strings(1, 10),".pdf")},
    content = function(file) {
      set.seed(1234)
      pdf(file,
          width = input$CellChatSignaling_plot_size$width / 100,
          height = input$CellChatSignaling_plot_size$height / 100 )
      plot(
        DoCellChatSignalingPlot(cellchat=InputValue$cellchat_less,source=InputValue$CellChat_source,target=InputValue$CellChat_target,signaling=InputValue$CellChat_signaling)
      )
      dev.off()
    }
  )
  output$CellChatSignaling_download2 <- downloadHandler(
    filename = function() {paste0(InputValue$dataset_choose,"-",InputValue$sample_choose,"_CellChatData",stringi::stri_rand_strings(1, 10),".csv")},
    content = function(file) {
      set.seed(1234)
      data <-  subsetCommunication(InputValue$cellchat_less,slot.name = 'net')
      write.csv(data, file, col.names = T, row.names = F, quote = F)
    }
  )
  
  ####other ----
  observeEvent(input$change_bdy, {
    # change_page("other?platform=Visium&dataset=BRCA_10x&sample=BRCA_BlockASection2_10x&selectedTab=Spatial_Gene_Expression")
    change_page("other?platform=Slide_seq&dataset=Melanoma_GSE185386&sample=GSM6025946&selectedTab=Spatial_Gene_Expression")
  })
  observeEvent(input$change_total, {
     change_page("other?platform=Visium&dataset=BRCA_10x&sample=BRCA_BlockASection2_10x&selectedTab=Spatial_Gene_Expression")
    #change_page("other?platform=Slide&dataset=Melanoma_GSE185386&sample=GSM6025946&selectedTab=Spatial_Gene_Expression")
  })
}

# Run server in a standard way.
shinyApp(ui, server)