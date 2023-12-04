## Spatial Structure ----
BoundaryInput <- function(){
  fluidPage(
    wellPanel(
      div(
        width='100%',height='500px',
        h4("Plot parameters"),
        selectInput(inputId = "BoundaryPlot_type",label = "Structure type",choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis")),
        numericInput(inputId = "BoundaryPlot_spotsize",label = "Spot size",value = 1.6,min = 0.8,max = 2,step = 0.1),
        numericInput(inputId = "BoundaryPlot_alpha",label = "Alpha",value = 1,min = 0, max = 1,step = 0.1),
        numericInput(inputId = "BoundaryPlot_legendsize",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1),
        fluidRow(
          column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "Boundary_submit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
        )
      )
    )
  )
}

BoundaryPlot <- function() {
  fluidPage(
    wellPanel(
      fluidRow(align='center',
               jqui_resizable(plotOutput('Boundary_plot', width = '600px', height = '500px'))),
      fluidRow(align='right',downloadButton(outputId = 'Boundary_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

## Spatial Expression ----
GeneInput <- function(){
  fluidPage(
    wellPanel(
      div(
        width='100%',height='500px',
        h4("Gene input"),
        textAreaInput(
          inputId = "GeneInput_text",
          label = "Type a gene or geneset",
          value = NULL
        ),
        selectizeInput(
          inputId = "GeneInput_saved",
          label = "Select from the saved geneset:",
          multiple = TRUE,
          choice = Saved_genes_panel,
          selected = NULL
        ),
        # textAreaInput(inputId = "GeneInput_text", label = "Type a gene or geneset", value = NULL),
        h4("Plot parameters"),
        selectInput(inputId = "PlotPar_facet",label = "Multi gene",choices = c("Seperate", "Geometric mean"),selected = "Seperate"),
        numericInput(inputId = "PlotPar_ncol",label = "Col number",value = 2),
        numericInput(inputId = "PlotPar_spotsize",label = "Spot size",value = 1.6,min = 0.8,max = 2,step = 0.1),
        fluidRow(
          column(6,align="center",withBusyIndicatorUI(actionButton(inputId = "GeneExpr_clear",label = "Clear",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
          column(6,align="center",withBusyIndicatorUI(actionButton(inputId = "GeneExpr_submit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
          )
      )
    )
  )
}

GenePlot <- function() {
  fluidPage(
    wellPanel(
      fluidRow(align='center',
               jqui_resizable(plotOutput('GeneExpr_plot', width = '500px', height = '500px'))),
      fluidRow(align='right',
               downloadButton(outputId = 'GeneExpr_download1',label = 'Download figure',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

## DE ----
DEInput <- function(){
  fluidPage(
    wellPanel(
      div(
        width='100%',height='500px',
        selectInput("DE_type", label = tags$span("Analysis type", class = "select-label"), choices = c("DEG", "Functions"), selected = "DEG"),
        h4("Plot parameters"),
        conditionalPanel(
          condition = "input.DE_type == 'DEG'",
          selectInput("DE_DEG_structure", "Structure type", choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis")),
          selectInput(inputId = "DE_DEG_cluster",label = "Location",choices = NULL),
          fluidRow(
            column(6,align="left",numericInput(inputId = "DE_DEG_pvalue",label = "Cut-off: pvalue",value = 0.05,min = 0.01, max = 0.1,step = 0.01)),
            column(6,align="left",numericInput(inputId = "DE_DEG_FC",label = "Cut-off: logFC",value = 0.25,min = 0, max = 100,step = 0.01))
          ),
          fluidRow(
            column(6,align="left",numericInput(inputId = "DE_DEG_N",label = "Labels",value = 10,min = 1, max = 100,step = 1)),
            column(6,align="left",numericInput(inputId = "DE_DEG_legendsize",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1))
          ),
          fluidRow(
            column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "DE_DEGsubmit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
          )
        ),
        conditionalPanel(
          condition = "input.DE_type == 'Functions'",
          selectInput("DE_GO_structure", "Structure type", choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis")),
          selectInput(inputId = "DE_GO_cluster",label = "Location",choices = NULL),
          fluidRow(
            column(6,align="left",numericInput(inputId = "DE_GO_pvalue",label = "Cut-off: pvalue",value = 0.05,min = 0.01, max = 0.1,step = 0.01)),
            column(6,align="left",numericInput(inputId = "DE_GO_FC",label = "Cut-off: logFC",value = 0.25,min = 0, max = 100,step = 0.01))
          ),
          fluidRow(
            column(6,align="left",numericInput(inputId = "DE_GO_N",label = "Pathways",value = 20,min = 1, max = 100,step = 1)),
            column(6,align="left",numericInput(inputId = "DE_GO_legendsize",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1))
          ),
          fluidRow(
            column(6,align="center",withBusyIndicatorUI(actionButton(inputId = "DE_GOanalysis",label = "Analysis",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))),
            column(6,align="center",withBusyIndicatorUI(actionButton(inputId = "DE_GOplot",label = "Plot",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
          )
        )
      )
    )
  )
}

DEPlot <- function() {
  fluidPage(
    conditionalPanel(
      condition = "input.DE_type == 'DEG'",
      fluidRow(align='center',
               jqui_resizable(plotOutput('DEValcano_plot', width = '500px', height = '500px'))),
      fluidRow(align='right',
               downloadButton(outputId = 'DEValcano_download1',label = 'Download figure',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
               downloadButton(outputId = 'DEValcano_download2',label = 'Download table',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    ),
    conditionalPanel(
      condition = "input.DE_type == 'Functions'",
      fluidRow(align='center',
               jqui_resizable(plotOutput('DEGO_plot', width = '500px', height = '700px'))),
      fluidRow(align='right',
               downloadButton(outputId = 'DEGO_download1',label = 'Download figure',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
               downloadButton(outputId = 'DEGO_download2',label = 'Download table',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

## Cell Distribution ----
DeconInput <- function(){
  fluidPage(
    wellPanel(
      div(
        width='100%',height='500px',
        selectInput("Decon_type", label = tags$span("Visualization type", class = "select-label"), choices = c("Spot composition", "Cell location"), selected = "Spot composition"),
        h4("Plot parameters"),
        conditionalPanel(
          condition = "input.Decon_type == 'Spot composition'",
          selectInput("Decon_structure", "Structure type", choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis")),
          selectInput(inputId = "DeconPie_spot",label = "Spots for plot",choices = NULL),
          numericInput(inputId = "DeconBar_legendsize",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1),
          fluidRow(
            column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "Decon_submit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))          
          )
        ),
        conditionalPanel(
          condition = "input.Decon_type == 'Cell location'",
          selectInput("CellLocation_structure", "Structure type", choices = c("Morphological adjusted cluster","Mal-Bdy-nMal axis")),
          selectInput(inputId = "CellLocation_spot",label = "Spots for plot",choices = NULL),
          selectInput(inputId = "CellLocation_cell",label = "Cells for plot",choices = NULL),
          fluidRow(
            column(6,align="left",numericInput(inputId = "CellLocation_spotsize",label = "Spot size",value = 1.6,min = 0.8,max = 2,step = 0.1)),
            column(6,align="left",numericInput(inputId = "CellLocationBox_legendsize",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1))
          ),
          selectInput(inputId = "CellLocation_test",label = "Statistical analysis",choices = NULL),
          fluidRow(
            column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "CellLocation_submit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
          )
        )
      )
    )
  )
}

DeconPlot <- function() {
  fluidPage(
    conditionalPanel(
      condition = "input.Decon_type == 'Spot composition'",
      fluidRow(align='center',
               jqui_resizable(plotOutput('DeconPie_plot', width = '500px', height = '500px'))),
      fluidRow(align='right',
               downloadButton(outputId = 'DeconPie_download1',label = 'Download figure',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
               downloadButton(outputId = 'DeconPie_download2',label = 'Download table',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    ),
    conditionalPanel(
      condition = "input.Decon_type == 'Cell location'",
      fluidRow(align='center',
               jqui_resizable(plotOutput('CellLocation_plot', width = '500px', height = '500px'))),
      fluidRow(align='right',
               downloadButton(outputId = 'CellLocation_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

DeconStatPlot <- function() {
  fluidPage(
    conditionalPanel(
      condition = "input.Decon_type == 'Spot composition'",
      fluidRow(align='center',div(style = "height:500px; display: flex; align-items: center; justify-content: center;",
               jqui_resizable(plotOutput('DeconBar_plot', width = '300px', height = '300px')))),
      fluidRow(align='right',downloadButton(outputId = 'DeconBar_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    ),
    conditionalPanel(
      condition = "input.Decon_type == 'Cell location'",
      fluidRow(align='center',div(style = "height:500px;display: flex; align-items: center; justify-content: center;",
               jqui_resizable(plotOutput('CellLocationBox_plot', width = '300px', height = '300px')))),
      fluidRow(align='right',downloadButton(outputId = 'CellLocationBox_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    ),
    conditionalPanel(
      condition = "input.Decon_type == 'Co-location'",
      fluidRow(align='center',div(style = "height:500px;display: flex; align-items: center; justify-content: center;",
                                  jqui_resizable(plotOutput('CoLocation_plot', width = '300px', height = '300px')))),
      fluidRow(align='right',downloadButton(outputId = 'CoLocation_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

## Cell Distribution-slide ----
DeconInput_slide <- function(){
  fluidPage(
    wellPanel(
      div(
        width='100%',height='500px',
        h4("Plot parameters"),
        selectInput(inputId = "CellLocation_spot_slide",label = "Spots for plot",choices = NULL),
        selectInput(inputId = "CellLocation_cell_slide",label = "Cells for plot",choices = NULL),
        fluidRow(
          column(12,align="left",numericInput(inputId = "CellLocation_spotsize_slide",label = "Spot size",value = 1.6,min = 0.8,max = 2,step = 0.1))
        ),
        fluidRow(
          column(12,align="left",numericInput(inputId = "CellLocationBox_legendsize_slide",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1))
        ),
        fluidRow(
          column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "CellLocation_submit_slide",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
        )
      )
    )
  )
}

DeconPlot_slide <- function() {
  fluidPage(
    fluidRow(align='center',
             jqui_resizable(plotOutput('CellLocation_plot_slide', width = '500px', height = '500px'))),
    fluidRow(align='right',
             downloadButton(outputId = 'CellLocation_slide_download1',label = 'Download figure',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
             downloadButton(outputId = 'CellLocation_slide_download2',label = 'Download table',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
  )
}

DeconStatPlot_slide <- function() {
  fluidPage(
    fluidRow(align='center',div(style = "height:500px;display: flex; align-items: center; justify-content: center;",
                                jqui_resizable(plotOutput('CellLocationBox_plot_slide', width = '300px', height = '300px')))),
    fluidRow(align='right',downloadButton(outputId = 'CellLocationBox_slide_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
  )
}

## Reconstruction ----
ReconInput <- function(){
  fluidPage(
    wellPanel(
      div(
        width='100%',height='500px',
        selectInput("Recon_type", label = tags$span("Visualization type", class = "select-label"), choices = c("UMAP Plot", "Gene Expression Plot"), selected = "UMAP Plot"),
        conditionalPanel(
          condition = "input.Recon_type == 'UMAP Plot'",
          h4("Plot parameters"),
          selectInput("Recon_groupby",label = "Group by",choices = c("Cell type" = 'Subtypes',"Morphological adjusted clusters" = 'Morph_clusters','Morphological adjusted clusters_Cell type'='Morph_clusters_celltype')),
          numericInput(inputId = "Recon_legendsize",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1),
          fluidRow(
            column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "Recon_submit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
          )
        ),
        conditionalPanel(
          condition = "input.Recon_type == 'Gene Expression Plot'",
          h4("Gene input"),
          textAreaInput(inputId = "Recon_gene", label = "Type a gene or geneset", value = NULL),
          h4("Plot parameters"),
          selectInput(inputId = "ReconGene_facet",label = "Multi gene",choices = c("Seperate", "Geometric mean"),selected = "Seperate"),
          numericInput(inputId = "ReconGene_ncol",label = "Col number",value = 2),
          numericInput(inputId = "ReconGene_legendsize",label = "Legend size",value = 1,min = 0.1, max = 2,step = 0.1),
          fluidRow(
            column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "ReconGene_submit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
          )
        )
      )
    )
  )
}

ReconPlot <- function() {
  fluidPage(
    conditionalPanel(
      condition = "input.Recon_type == 'UMAP Plot'",
      fluidRow(align='center',
               jqui_resizable(plotOutput('Recon_plot', width = '500px', height = '500px'))),
      fluidRow(align='right',downloadButton(outputId = 'Recon_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    ),
    conditionalPanel(
      condition = "input.Recon_type == 'Gene Expression Plot'",
      fluidRow(align='center',
               jqui_resizable(plotOutput('ReconGene_plot', width = '500px', height = '500px'))),
      fluidRow(align='right',downloadButton(outputId = 'ReconGene_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

## Cell Interaction ----
CellChatInput <- function(){
  fluidPage(
    wellPanel(
      div(
        width='100%',height='500px',
        h4("Plot parameters"),
        selectInput(inputId = "CellChat_source",label = "Source",choices = NULL),
        selectInput(inputId = "CellChat_target",label = "Target",choices = NULL),
        selectInput(inputId = "CellChat_signaling",label = "Signaling",choices = NULL),
        fluidRow(
          column(12,align="center",withBusyIndicatorUI(actionButton(inputId = "CellChat_submit",label = "Submit",width = "80px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4")))
        )
      )
    )
  )
}

CellChatPlot <- function() {
  fluidPage(
    wellPanel(
      fluidRow(align='center',
               jqui_resizable(plotOutput('CellChat_plot', width = '300px', height = '400px'))),
      fluidRow(align='right',downloadButton(outputId = 'CellChat_download1',label = 'Download',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

CellChatSignalingPlot <- function(){
  fluidPage(
    wellPanel(
      fluidRow(align='center',
               jqui_resizable(plotOutput('CellChatSignaling_plot', width = '400px', height = '400px'))),
      fluidRow(align='right',
               downloadButton(outputId = 'CellChatSignaling_download1',label = 'Download figure',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"),
               downloadButton(outputId = 'CellChatSignaling_download2',label = 'Download table',width = "120px",style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"))
    )
  )
}

