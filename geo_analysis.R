#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# To run: open this file in RStudio and click 'Run App', or run
#    shiny::runApp("path/to/geo_analysis.R")
# Any missing R packages are installed automatically the first time the app is run.
# The 'Top Gene Pathways' tab also needs Java 8 or newer installed (see the tab for details).
# Tests: testthat::test_dir("tests/testthat")
#

required_pkgs <- c("shiny", "plyr", "tidyverse", "readr", "cli", "GEOquery", "shinyWidgets",
                   "shinycssloaders", "limma", "dplyr", "stringr", "DescTools", "ggplot2",
                   "gridExtra", "umap", "pathfindR", "tibble", "ropls", "Biobase", "DT",
                   "shinyFiles", "data.table")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(missing_pkgs, update = FALSE, ask = FALSE)
}

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
# pathfindR is not attached here: attaching it fails without Java and would stop the whole app.
# Its functions are called as pathfindR::... in diff_cell.R only when the pathway tab is used.
library(tibble)
library(ropls)
library(Biobase)
library(DT)
library(shinyFiles)
library(data.table)

source("diff_cell.R")

app_version <- "1.1.0"
java_available <- nzchar(Sys.which("java"))  # pathfindR needs Java on the PATH that R sees

#### UI ####


# Define UI 
ui <- fluidPage(

    # Application title
    titlePanel("GEO Accession Data Analysis"),

    # SIDEBAR 
    sidebarLayout(
        sidebarPanel(
            helpText("Enter a GEO Series accession (e.g. GSE148426) and press 'Load Data'.",
                     "Then choose a classification variable and press 'Run Analysis'.",
                     "Loading a new accession clears all previous results."),
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
                     helpText("Sample information for the loaded series, with GEO administrative columns removed.",
                              "Variables with two or more groups can be chosen as the classification variable in the side bar."),
                     DTOutput("geodata") %>% withSpinner(),
                     br(),
                     downloadButton("phenotabsave", "Download Data Table")),
            tabPanel("Top Gene Expression Table", 
                     helpText("limma results for the pairwise comparison selected in the side bar, ordered by the B statistic",
                              "(log-odds that the gene is differentially expressed). logFC is the log2 fold change of the first",
                              "group relative to the second group in the comparison name (e.g. A-B means A minus B).",
                              "adj.P.Val is the Benjamini-Hochberg adjusted p-value, i.e. the false discovery rate (FDR)."),
                     DTOutput("genetables") %>% withSpinner(),
                     br(),
                     downloadButton("geneexpsave", "Download Data Table")),
            tabPanel("Histogram of P-Values", 
                     helpText("Distribution of Benjamini-Hochberg adjusted p-values (FDR) for all genes, as in GEO2R.",
                              "With two groups this is the single comparison; with three or more groups it is the moderated",
                              "F-test across all groups. A large peak near 0 indicates many differentially expressed genes."),
                     plotOutput("genehist") %>% withSpinner(),
                     downloadButton("histsave", "Download Graph")),
            tabPanel("Q-Q Plot", 
                     helpText("The moderated t-statistic of every gene with a valid test result is plotted against the",
                              "quantiles of a Student's t-distribution with the same degrees of freedom. If most genes are",
                              "not differentially expressed, their moderated t-statistics follow this t-distribution, so the",
                              "points lie on the diagonal; points that leave the line at either end are differentially expressed",
                              "genes. This checks the test statistics, not whether the expression data are normally distributed.",
                              "With three or more groups the first comparison is shown."),
                     plotOutput("geneqq") %>% withSpinner(), 
                     br(),
                     downloadButton("qqplotsave", "Download Graph")),
            tabPanel("Volcano Plots", 
                     helpText("Each point is a gene: log2 fold change against -log10(p-value). Genes are labelled UP or DOWN",
                              "when their Benjamini-Hochberg adjusted p-value is below the chosen false discovery rate (FDR)",
                              "threshold. 0.05 is a common choice, but any threshold is a convention, so choose one to suit your study."),
                     fluidRow(uiOutput("volcanopval")),
                     plotOutput("genevolcano") %>% withSpinner(),
                     br(),
                     downloadButton("volplotsave", "Download Graph"), 
                     br(),
                     DTOutput("volctop"),
                     br(),
                     downloadButton("voltabsave", "Download Data Table")),
            tabPanel("Mean Difference Plots", 
                     helpText("Each point is a gene: log2 fold change against average log2 expression. Genes are labelled UP",
                              "or DOWN using the chosen false discovery rate (FDR) threshold, as in the volcano plot."),
                     fluidRow(uiOutput("meandiffpval")), 
                     plotOutput("genemd") %>% withSpinner(),
                     br(),
                     downloadButton("mdplotsave", "Download Graph"),
                     br(),
                     DTOutput("mdtop"),
                     br(),
                     downloadButton("mdtabsave", "Download Data Table")),
            tabPanel("Density and Intensity Plot", 
                     helpText("Distribution of expression values in each sample, coloured by group. Similar curves suggest",
                              "the samples are comparably normalised; a clearly shifted curve can indicate a problem sample."),
                     plotOutput("genedense") %>% withSpinner(),
                     br(),
                     downloadButton("disave", "Download Graph")),
            tabPanel("UMAP",
                     helpText("UMAP is a non-linear dimension reduction that places samples with similar expression profiles",
                              "close together. The layout depends strongly on the number of neighbours (small values emphasise",
                              "local structure, large values global structure), and distances between clusters are not directly",
                              "meaningful, so clusters should not be read as evidence of biological relatedness on their own.",
                              "The default of 15 neighbours matches GEO2R.",
                              tags$a(href = "https://pair-code.github.io/understanding-umap/", target = "_blank",
                                     "Understanding UMAP")),
                     fluidRow(uiOutput("neighbours")),
                     br(),
                     textOutput("umaperr"),
                     plotOutput("geneumap"),
                     br(),
                     downloadButton("umapsave", "Download Graph")),
            tabPanel("Mean Variance Trend", 
                     helpText("Square root of the residual standard deviation of each gene against its average log2",
                              "expression (limma plotSA). The line is the mean-variance trend used by limma's empirical",
                              "Bayes moderation of the gene variances."),
                     plotOutput("genemeanvar"),
                     br(),
                     downloadButton("mvsave", "Download Graph")),
            tabPanel("Top Gene Pathways",
                     helpText("Pathway enrichment using active subnetworks (pathfindR), based on the gene identifiers and",
                              "adjusted p-values. Requires Java 8 or newer. Gene identifiers must be gene symbols to match the",
                              "protein interaction network, so datasets whose rows are probe IDs may find no match."),
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
                     helpText("Partial least squares discriminant analysis (ropls) of the expression data, using the",
                              "classification variable as the response. If a significant model is found the sample score plot",
                              "is shown; otherwise the model overview is shown."),
                     fluidRow(uiOutput("plsdabutt")),
                     plotOutput("geneplsda"),
                     br(),
                     downloadButton("plsdasave", "Download Graph")),
            tabPanel("About",
                     br(),
                     helpText("App version, dataset, settings and software versions for the current analysis.",
                              "Download this file and keep it with your results so the analysis can be reproduced."),
                     verbatimTextOutput("about"),
                     downloadButton("aboutsave", "Download Analysis Details"))
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
  
  geo_info <- reactiveVal(NULL)  # accession and download time, for the About tab
  
  observeEvent(input$button, { 
    geo_info(list(accession = input$geoacc, loaded = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
    
    # Clear all results from a previous analysis, so the app never shows results for an
    # earlier accession next to newly loaded data (no restart of the app needed)
    for (id in c("genects", "paths", "geneset", "neighbours", "volcanopval", "meandiffpval", "plsdabutt")) {
      output[[id]] <- renderUI(NULL)
    }
    for (id in c("genehist", "geneqq", "genedense", "genemeanvar", "geneumap", "genevolcano",
                 "genemd", "genepath", "geneplsda")) {
      output[[id]] <- renderPlot(NULL)
    }
    for (id in c("genetables", "volctop", "mdtop", "genepathdata")) {
      output[[id]] <- renderDataTable(NULL)
    }
    for (id in c("umaperr", "patherror")) {
      output[[id]] <- renderText(NULL)
    }
    for (id in c("geneexpsave", "histsave", "qqplotsave", "volplotsave", "voltabsave", "mdplotsave",
                 "mdtabsave", "disave", "umapsave", "mvsave", "pathplotsave", "pathtabsave", "plsdasave")) {
      output[[id]] <- downloadHandler(
        filename = "No_Results.txt",
        content = function(file) writeLines("No results yet: run this analysis for the loaded dataset first.", file))
    }
    
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
      if (!java_available) {
        return(p(strong("Pathway analysis is unavailable because Java was not found."),
                 "pathfindR needs Java 8 or newer, installed and on the PATH used by R",
                 "(check with Sys.which('java') in R). Install Java, restart R and reload the app."))
      }
      sets <- c("KEGG", "Reactome", "BioCarta", "GO-All", "GO-BP", "GO-CC", "GO-MF", "cell_markers", "mmu_KEGG")
      fluidPage(selectInput("set", 
                  label = "Select Database for Gene Search", 
                  choices = as.list(sets)),
      actionButton(inputId="term_button", label="Generate Top Gene Pathways"))
    })

    
    output$neighbours <- renderUI({
      fluidPage(numericInput("nneighs", label="Number of neighbours for UMAP analysis", 15, min=3, max=40, step=1),
                actionButton(inputId="neighs_button", label="Run UMAP Analysis"))
    })
    
    output$volcanopval <- renderUI({
      fluidPage(numericInput("volpval", label="Choose False Discovery Rate (FDR) Threshold", 0.05, min=0, max=1, step=0.005),
                actionButton(inputId="vol_button", label="Generate Volcano Plots"))
    })
    
    output$meandiffpval <- renderUI({
      fluidPage(numericInput("mdpval", label="Choose False Discovery Rate (FDR) Threshold", 0.05, min=0, max=1, step=0.005),
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

  about_text <- reactive({
    info <- geo_info()
    txt <- c(paste("GEO Accession Data Analysis app, version", app_version),
             paste("Report created:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")))
    if (!is.null(info)) {
      last_update <- tryCatch(
        max(as.character(pluck(GSE(), paste0(info$accession, "_series_matrix.txt.gz"))@phenoData@data$last_update_date)),
        error = function(e) "unknown")
      txt <- c(txt, "",
               paste("GEO accession:", info$accession),
               paste("Data downloaded from GEO:", info$loaded),
               paste("Most recent sample update date in GEO:", last_update),
               paste("Classification variable:", if (is.null(input$meta)) "not chosen" else input$meta),
               paste("Volcano plot FDR threshold:", if (is.null(input$volpval)) "not set" else input$volpval),
               paste("Mean difference plot FDR threshold:", if (is.null(input$mdpval)) "not set" else input$mdpval),
               paste("UMAP neighbours:", if (is.null(input$nneighs)) "not set" else input$nneighs),
               paste("Pathway database:", if (is.null(input$set)) "not set" else input$set))
    }
    c(txt, "",
      "Methods: limma (lmFit, contrasts.fit, eBayes with proportion = 0.01, topTable with Benjamini-Hochberg",
      "adjustment), following the GEO2R script; pathfindR for pathway enrichment; ropls for PLS-DA.",
      paste("Java available for pathfindR:", java_available),
      "", "Software versions:", capture.output(sessionInfo()))
  })
  
  output$about <- renderText({
    paste(about_text(), collapse = "\n")
  })
  
  output$aboutsave <- downloadHandler(
    filename = function() {
      paste0("Analysis_Details_", if (nzchar(input$geoacc)) input$geoacc else "no_dataset", ".txt")
    },
    content = function(file) {
      writeLines(about_text(), file)
    }
  )
  
}

# Run the application 
shinyApp(ui = ui, server = server)

