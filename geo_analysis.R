#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#


library(shiny)
library(plyr)
library(tidyverse)
library(readr)
library(cli)
library(grid)
library(stats)
library(GEOquery)
library(shinyWidgets)
library(shinycssloaders)

library(limma)
library(dplyr)
library(stringr)
library(DescTools)
library(ggplot2)
library(gridExtra)
library(umap)
library(maptools)
library(pathfindR)
library(tibble)
library(ropls)
library(Biobase)
library(DT)
library(shinyFiles)
library(data.table)

source("diff_cell.R")

#### UI ####


# Define UI 
ui <- fluidPage(

    # Application title
    titlePanel("GEO Accession Data Analysis"),

    # SIDEBAR 
    sidebarLayout(
        sidebarPanel(
            textInput(inputId = "geoacc", 
                      label = "GEO Accession Code",
                      value = ""),
            
            actionButton(inputId="button", label="Load Data"),
            br(),
            fluidRow(
              br(),
              uiOutput("meta"),
              br(),
              uiOutput("genects")
            )
            # input GEO Accession ID
            # column(width = 7, progressBar(id = "pb1", value = 0, display_pct = T))
        ),
        

        # MAIN PANEL
        mainPanel(
          tabsetPanel(
            tabPanel("Phenotype Data", 
                     DTOutput("geodata"),
                     br(),
                     downloadButton("phenotabsave", "Download Data Table")),
            tabPanel("Top Gene Expression Table", 
                     DTOutput("genetables"),
                     br(),
                     downloadButton("geneexpsave", "Download Data Table")),
            tabPanel("Histogram of P-Values", 
                     plotOutput("genehist"),
                     downloadButton("histsave", "Download Graph")),
            tabPanel("Q-Q Plot", 
                     plotOutput("geneqq"), 
                     br(),
                     downloadButton("qqplotsave", "Download Graph")),
            tabPanel("Volcano Plots", 
                     fluidRow(uiOutput("volcanopval")),
                     plotOutput("genevolcano"),
                     br(),
                     downloadButton("volplotsave", "Download Graph"), 
                     br(),
                     DTOutput("volctop"),
                     br(),
                     downloadButton("voltabsave", "Download Data Table")),
            tabPanel("Mean Difference Plots", 
                     fluidRow(uiOutput("meandiffpval")), 
                     plotOutput("genemd"),
                     br(),
                     downloadButton("mdplotsave", "Download Graph"),
                     br(),
                     DTOutput("mdtop"),
                     br(),
                     downloadButton("mdtabsave", "Download Data Table")),
            tabPanel("Density and Intensity Plot", 
                     plotOutput("genedense"),
                     br(),
                     downloadButton("disave", "Download Graph")),
            tabPanel("UMAP",
                     fluidRow(uiOutput("neighbours")),
                     br(),
                     textOutput("umaperr"),
                     plotOutput("geneumap"),
                     br(),
                     downloadButton("umapsave", "Download Graph")),
            tabPanel("Mean Variance Trend", 
                     plotOutput("genemeanvar"),
                     br(),
                     downloadButton("mvsave", "Download Graph")),
            tabPanel("Top Gene Pathways",
                     fluidRow(uiOutput("geneset"),
                              br(),
                              uiOutput("paths")),
                     br(),
                     textOutput("patherror"),
                     plotOutput("genepath"),
                     br(),
                     downloadButton("pathplotsave", "Download Graph"),
                     br(),
                     DTOutput("genepathdata"),
                     br(),
                     downloadButton("pathtabsave", "Download Data Table")),
            tabPanel("PLS-DA Analysis",
                     br(),
                     fluidRow(uiOutput("plsdabutt")),
                     plotOutput("geneplsda"),
                     br(),
                     downloadButton("plsdasave", "Download Graph"))
          )
           
        )
    )

)


#### SERVER ####

# the dynamic input phenodata should be for factor/character variables for further analysis.

# Define server 
server <- function(input, output, session) {
  Sys.setenv(VROOM_CONNECTION_SIZE=10000001)
  readr::local_edition(1)
  
  output$geodata <- renderDataTable(data.frame(col1=NA))
  

    GSE  <- eventReactive(input$button,{
      withProgress(message="Reading in GEO data", value=0,{
        incProgress(1/2, message = "Reading in GEO data")
      getGEO(input$geoacc, GSEMatrix =TRUE, AnnotGPL=FALSE)
      })
    })

    GSEpheno.clean <- eventReactive(input$button,{
      GSEpheno <- pluck(GSE(), paste0(input$geoacc,"_series_matrix.txt.gz"))@phenoData@data
      diff_cell(GSEpheno)
    })
    GSEassay  <- eventReactive(input$button,{
      GSEassay <- pluck(GSE(), paste0(input$geoacc,"_series_matrix.txt.gz"))@assayData[["exprs"]]
      as.data.frame(GSEassay)
    })
    GSE.t.assay  <- eventReactive(input$button,{
      GSE.t.assay <- t(GSEassay())
      as.data.frame(GSE.t.assay)
    })
    
    
    phenotab <- reactive({
      data.table(GSEpheno.clean(), keep.rownames=TRUE)
    })   
    
    
    output$geodata <- renderDataTable({
      phenotab()
      }, options=list(sDom = 'lfrtip', scrollY="600px", scrollX=TRUE))
    
    
    output$phenotabsave <- downloadHandler(
      filename = function() {
        paste("Phenotype_Data", ".csv", sep="") 
        }, 
      content = function(file) {
        fwrite(phenotab(), file, row.names = FALSE)
      }
    )
  
  observeEvent(input$button, { 
    output$meta <- renderUI({
      fluidPage(selectInput("meta", 
                            label = "Please Select Classification Variable for Analysis", 
                            choices = as.list(GSEpheno.clean() %>% select_if(~is.factor(.) & nlevels(.)>=2) %>% colnames())),
                actionButton(inputId="meta_button", label="Run Analysis"))
    })
  })
  
  observeEvent(input$meta_button, {
    GSEexp <- top_gene_exp(GSEpheno.clean(), GSE.t.assay(), GSEassay(), input$meta)
    GSEtable <- top_gene_table(GSEexp$tT.logfc, GSEexp$tT)
    
    output$genects <- renderUI({
      fluidPage(selectInput("ctschoice", 
                            label = "Please Select Pairwise Analysis to View Top Gene Expression Table and Associated Plots", 
                            choices = as.list(GSEexp$cts)))
    })
    
    
    genexptab <- reactive({
      if(is.null(GSEexp$tT)){
        data.table(GSEtable$tT.subset.logfc, keep.rownames=TRUE)
      } else{
        GSEcts <- GSEexp$cts
        ind <- which(GSEcts==input$ctschoice)
        data.table(GSEtable$tT.subset.logfc[[ind]], keep.rownames=TRUE)
      }
    })
    
    output$genetables <- renderDataTable({
      genexptab()
    }, options=list(sDom="lfrtip",scrollY="600px", scrollX=TRUE))
    
    output$geneexpsave <- downloadHandler(
      filename = function() {
        paste0(input$ctschoice, "_Top_Gene_Expression_Data", ".csv", sep="") 
      }, 
      content = function(file) {
        fwrite(genexptab(), file, row.names = FALSE)
      }
    )
    
    
    
    histplot <- function(){
      gene_hist(GSEpheno.clean(), input$meta, GSEexp$tT.logfc, GSEexp$tT)
    }
    
    output$genehist <- renderPlot({
      histplot()
    })
    
    output$histsave <- downloadHandler(
      filename = function(){
        paste("Histogram", ".png", sep="")
      },
      content = function(file){
        png(file)
        histplot()
        dev.off()
      }
    )
    
    qqplot <- function(){
      gene_qqplot(GSEexp$fit2)
    }
    
    output$geneqq <- renderPlot({
      qqplot()
    })
    
    output$qqplotsave <- downloadHandler(
      filename = function(){
        paste("QQ_Plot", ".png", sep="")
      },
      content = function(file){
        png(file)
        qqplot()
        dev.off()
      }
    )

    diplot <- function(){
      gene_densities(GSEassay(), GSEexp$gs)
    }
    
    output$genedense <- renderPlot({
      diplot()
    })
    
    output$disave <- downloadHandler(
      filename = function(){
        paste("Density_and_Intensity_Plot", ".png", sep="")
      },
      content = function(file){
        png(file)
        diplot()
        dev.off()
      }
    )
    
    mvplot <- function(){
      gene_meanvar(GSEexp$fit2, input$geoacc)
    }
    
    output$genemeanvar <- renderPlot({
      mvplot()
    })
    
    output$mvsave <- downloadHandler(
      filename = function(){
        paste("Mean_Variance_Trend", ".png", sep="")
      },
      content = function(file){
        png(file)
        mvplot()
        dev.off()
      }
    ) 
    
    
    
    output$paths <- renderUI({
      fluidPage(numericInput("nterms", 
                             label = "Number of Top Gene Pathways Displayed", 5, min=5, max=40, step=5)
      )
    })  
    
    output$geneset <- renderUI({
      sets <- c("KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC", "GO-MF", "cell_markers", "mmu_KEGG")
      fluidPage(selectInput("set", 
                  label = "Select Database for Gene Search", 
                  choices = as.list(sets)),
      actionButton(inputId="term_button", label="Generate Top Gene Pathways"))
    })

    
    output$neighbours <- renderUI({
      fluidPage(numericInput("nneighs", label="Number of neighbours for UMAP analysis", 3, min=3, max=40, step=1),
                actionButton(inputId="neighs_button", label="Run UMAP Analysis"))
    })
    
    output$volcanopval <- renderUI({
      fluidPage(numericInput("volpval", label="Choose Significance Threshold", 0.001, min=0.001, max=1, step=0.0005),
                actionButton(inputId="vol_button", label="Generate Volcano Plots"))
    })
    
    output$meandiffpval <- renderUI({
      fluidPage(numericInput("mdpval", label="Choose Significance Threshold", 0.001, min=0.001, max=1, step=0.0005),
                actionButton(inputId="md_button", label="Generate Mean Difference Plots"))
    })
    
    output$plsdabutt <- renderUI({
      fluidPage(actionButton(inputId="plsda_button", label="Run PLS-DA Analysis"))
    })
    
  })
  
  observeEvent(input$neighs_button, {
    GSEexp <- top_gene_exp(GSEpheno.clean(), GSE.t.assay(), GSEassay(), input$meta)
    if(input$nneighs >= nrow(GSEpheno.clean())){
      output$umaperr <- renderText({"Number of neighbours must be strictly smaller than the number of samples"})
    } else {
    output$umaperr <- renderText({NULL})
    
    umapplot <- function(){
      gene_umap(GSE.t.assay(), GSEexp$gs, nneighs=input$nneighs)
    }
    
    output$geneumap <- renderPlot({
      withProgress(message="Generating UMAP", value=0,{
        incProgress(1/2, message = "Generating UMAP")
        umapplot()
      })
    })
    
    
    output$umapsave <- downloadHandler(
      filename = function(){
        paste("UMAP", ".png", sep="")
      },
      content = function(file){
        ggsave(file, plot=umapplot(), width = 10, height = 7, dpi = 300, device='png')
      }
    )
    
    } 
  })
  

  observeEvent(input$vol_button, {
    GSEexp <- top_gene_exp(GSEpheno.clean(), GSE.t.assay(), GSEassay(), input$meta)
    
    volplot <- function(){
      if(is.null(GSEexp$tT)){
        gene_volcano(GSEexp$cts, GSEexp$tT.logfc, input$volpval)
      } else{
        GSEcts <- GSEexp$cts
        ind <- which(GSEcts==input$ctschoice)
        gene_volcano(GSEexp$cts, GSEexp$tT.logfc, input$volpval)[[ind]]
      }
    }
    
    output$genevolcano <- renderPlot({
      volplot()
    })
    
    output$volplotsave <- downloadHandler(
      filename = function(){
        paste(input$ctschoice, "_Volcano_Plot.png", sep="")
      },
      content = function(file){
        ggsave(file, plot=volplot(), width = 9, height = 6, dpi = 300, device='png')
      }
    )
    
    
    voltab <- reactive({
      if(is.null(GSEexp$tT)){
        gene_vol_table(GSEexp$cts, GSEexp$tT.logfc, input$volpval)
      } else{
        GSEcts <- GSEexp$cts
        ind <- which(GSEcts==input$ctschoice)
        gene_vol_table(GSEexp$cts, GSEexp$tT.logfc, input$volpval)[[ind]]
      }      
    })
    
    output$volctop <- renderDataTable({
      voltab()
    })
    
    output$voltabsave <- downloadHandler(
      filename = function() {
        paste0(input$ctschoice, "_Volcano_Data", ".csv", sep="") 
      }, 
      content = function(file) {
        fwrite(voltab(), file, row.names = TRUE)
      }
    )
    
    
  })
  
  observeEvent(input$md_button, {
    GSEexp <- top_gene_exp(GSEpheno.clean(), GSE.t.assay(), GSEassay(), input$meta)
    
    mdplot <- function(){
      if(is.null(GSEexp$tT)){
        gene_md(GSEexp$cts, GSEexp$tT.logfc, input$mdpval)
      } else{
        GSEcts <- GSEexp$cts
        ind <- which(GSEcts==input$ctschoice)
        gene_md(GSEexp$cts, GSEexp$tT.logfc, input$mdpval)[[ind]]
      }
    }      
    
    output$genemd <- renderPlot({
      mdplot()
    })
    
    
    output$mdplotsave <- downloadHandler(
      filename = function(){
        paste(input$ctschoice, "_Mean_Difference_Plot.png", sep="")
      },
      content = function(file){
        ggsave(file, plot=mdplot(), width = 9, height = 6, dpi = 300, device='png')
      }
    )
    
    mdtab <- reactive({
      if(is.null(GSEexp$tT)){
        gene_md_table(GSEexp$cts, GSEexp$tT.logfc, input$mdpval)
      } else{
        GSEcts <- GSEexp$cts
        ind <- which(GSEcts==input$ctschoice)
        gene_md_table(GSEexp$cts, GSEexp$tT.logfc, input$mdpval)[[ind]]
      }      
    })
    
    output$mdtop <- renderDataTable({
      mdtab()
    })
    
    output$mdtabsave <- downloadHandler(
      filename = function() {
        paste0(input$ctschoice, "_Mean_Difference_Data", ".csv", sep="") 
      }, 
      content = function(file) {
        fwrite(mdtab(), file, row.names = TRUE)
      }
    )
    
  })
  
  observeEvent(input$term_button, {
    
    GSEexp <- top_gene_exp(GSEpheno.clean(), GSE.t.assay(), GSEassay(), input$meta)
    withProgress(message="Generating Top Gene Pathways", value=0,{
      incProgress(1/4, message = "Generating Top Gene Pathways")
      
        output_df <- gene_paths(GSEexp$tT.logfc, GSEexp$tT, input$set)
        if(is.data.frame(output_df)){
          
          pathplot <- function(){
            gene_paths_plot(output_df, nterms = input$nterms)
          }
          
          
          pathtab <- reactive({
            data.table(output_df)
          })
 
        } else {
          pathplot <- function(){NULL}
          pathtab <- function(){NULL}
          err <- "Cannot match gene symbols to PIN."
          output$patherror <- renderText({err})
                 
        }
         
        incProgress(1/4, message = "Generating Top Gene Pathways")
        output$genepath <- renderPlot({
          pathplot()
        })
        
        output$pathplotsave <- downloadHandler(
          filename = function(){
            paste0("Top_", input$nterms, "_Gene_Pathways.png", sep="")
          },
          content = function(file){
            ggsave(file, plot=pathplot(), width = 9, height = 6, dpi = 300, device='png')
          }
        )
        
        incProgress(1/4, message = "Generating Top Gene Pathways") 
        output$genepathdata <- renderDataTable({
          pathtab()
        }, options=list(scrollY="600px", scrollX=TRUE))
        
        
        output$pathtabsave <- downloadHandler(
          filename = function() {
            paste("Top_Gene_Pathway_Data", ".csv", sep="") 
          }, 
          content = function(file) {
            fwrite(pathtab(), file, row.names = FALSE)
          }
        )          
        
        
           
      })
  })
  
  observeEvent(input$plsda_button,{
    withProgress(message="Generating PLS-DA Analysis", value=0,{
      
      incProgress(1/2, message = "Generating PLS-DA Analysis")
      plsda_data <- gene_plsda(GSEassay(), GSEpheno.clean(), input$meta)
      
      plsdaplot <- function(){
        gene_plsda_plot(plsda_data$sacPlsda)
      }
      
      output$geneplsda <- renderPlot({
      plsdaplot()
    })
      
      output$plsdasave <- downloadHandler(
        filename = function(){
          paste("PLS-DA_Plot", ".png", sep="")
        },
        content = function(file){
          png(file)
          plsdaplot()
          dev.off()
        }
      ) 
      
      
  })  
  })


  
}

# Run the application 
shinyApp(ui = ui, server = server)
