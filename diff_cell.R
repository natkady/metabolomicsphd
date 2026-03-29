#function for if statement for pheno clean with diff cell count



diff_cell <- function(GSEpheno){
  useless_cols <-  c("title","geo_accession","status","submission_date","last_update_date","type","channel_count", "source_name_ch1","organism_ch1","characteristics_ch1","characteristics_ch1.1","molecule_ch1","extract_protocol_ch1","label_ch1","label_protocol_ch1","taxid_ch1","hyb_protocol","scan_protocol","data_processing","platform_id","contact_name","contact_institute","contact_address","contact_city","contact_zip/postal_code","contact_country","supplementary_file","data_row_count")
  if("diff cell count:ch1" %in% colnames(GSEpheno)){
    GSEpheno.nodiff <- str_split(GSEpheno$`diff cell count:ch1`, ";")
    GSEpheno.nodiff.frame <- do.call(rbind.data.frame, GSEpheno.nodiff)
    colvec <- c()
      for (i in 1:length(GSEpheno.nodiff.frame)){
        colvec[i] <- unlist(str_split(GSEpheno.nodiff.frame[1,i], ":"))[1]
      }
    colnames(GSEpheno.nodiff.frame) <- colvec
    GSEpheno.nodiff.new <- GSEpheno.nodiff.frame
      for (i in 1:ncol(GSEpheno.nodiff.frame)){
        GSEpheno.nodiff.new[,i] <- separate(GSEpheno.nodiff.frame,i,into=c(NA,"number"),sep=": ") %>% pull(number)
      }
    GSEpheno.test <- GSEpheno %>% select(-c(all_of(useless_cols), contains("characteristics_ch"), contains("contact"), contains("_ch1") , contains("diff cell count:ch1")))
    GSEnodiffclean <- bind_cols(GSEpheno.nodiff.new, GSEpheno.test)
    colnames(GSEnodiffclean) <- sub(":.*", "", colnames(GSEnodiffclean)) 
    GSEpheno.clean <- GSEnodiffclean %>% mutate_if(~all(str_detect(.,"^[0-9]+[\\.]?[0-9]*$")==TRUE),~as.numeric(.))
  } else {
    GSEpheno.clean <- GSEpheno %>% select(-c(all_of(useless_cols), contains("characteristics_ch"), contains("contact"), contains("_ch1")))
    colnames(GSEpheno.clean) <- sub(":.*", "", colnames(GSEpheno.clean))
    GSEpheno.clean <- GSEpheno.clean %>% mutate_if(~all(str_detect(.,"^[0-9]+[\\.]?[0-9]*$")==TRUE),~as.numeric(.))
  }
  GSEpheno.clean <- GSEpheno.clean %>% mutate_if(~is.character(.), ~as.factor(.))
}


# meta <- pull(GSEpheno.clean, var=meta)
# function(GSEpheno.clean, GSE.t.assay, GSEassay, meta=input$meta)


top_gene_exp <- function(GSEpheno.clean, GSE.t.assay, GSEassay, meta){
  
  GSEdata <- bind_cols(GSEpheno.clean, GSE.t.assay)
  GSEassay <- log10(GSEassay + 1 - min(GSEassay))
  any(is.na(GSEassay))
  
  
  gs <- factor(pull(GSEpheno.clean, meta))
  groups <- make.names(unique(factor(pull(GSEpheno.clean, meta))))
  levels(gs) <- groups
  GSE.t.assay$group <- gs 
  design <- model.matrix(~group + 0, GSE.t.assay)
  colnames(design) <- levels(gs)
  
  fit <- lmFit(GSEassay, design)
  
  if (length(unique(pull(GSEpheno.clean, meta)))==2){
    cts <- paste(CombPairs(groups)[,1], CombPairs(groups)[,2], sep="-")
    cont.matrix <- makeContrasts(contrasts=cts, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT.logfc <- topTable(fit2, adjust.method="fdr", sort.by="B", number=Inf)
    #tT.subset.logfc <- subset(tT.logfc, select=c("adj.P.Val","P.Value","t","B","logFC"))
    output <- list(cts=cts, fit2=fit2, tT.logfc=tT.logfc, tT=NULL, gs=gs)
    #write.table(tT.subset.logfc, file = paste("Contrast_", cts, ".txt"), row.names = TRUE)
  } else {
    cts <- paste(CombPairs(groups)[,1], CombPairs(groups)[,2], sep="-")
    cont.matrix <- makeContrasts(contrasts=cts, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    tT <- topTable(fit2, adjust.method="fdr", sort.by="B", number=Inf)
    
    tT.subset <- subset(tT, select=c("adj.P.Val","P.Value","F")) 
    cont.matrix.logfc <- NULL
    fit2.logfc <- NULL
    tT.logfc <- NULL
    for(i in 1:length(cts)){
      cont.matrix.logfc[[i]] <- makeContrasts(contrasts=cts[i], levels=design)
      fit2.logfc[[i]] <- contrasts.fit(fit, cont.matrix.logfc[[i]])
      fit2.logfc[[i]] <- eBayes(fit2.logfc[[i]], 0.01)
      tT.logfc[[i]] <- topTable(fit2.logfc[[i]], adjust.method="fdr", sort.by="B", number=Inf)
      #tT.subset.logfc[[i]] <- subset(tT.logfc[[i]], select=c("adj.P.Val","P.Value","t","B","logFC"))
      #write.table(tT.subset.logfc[[i]], file = paste("Contrast_", cts[i], ".txt"), row.names = TRUE)
    }
    output <- list(cts=cts, fit2=fit2, tT.logfc=tT.logfc, tT=tT, gs=gs)
  }
  return(output)
}

# object <- top_gene_exp(inputs)
# object$cts

top_gene_table <- function(tT.logfc, tT){
  tT.subset <- NULL
  tT.subset.logfc <- NULL
  if (is.null(tT)){
    tT.subset.logfc <- tT.logfc %>% select("adj.P.Val","P.Value","t","B","logFC")}
  else{
    tT.subset <- tT %>% select("adj.P.Val","P.Value","F")
    for(i in 1:length(tT.logfc)){
    tT.subset.logfc[[i]] <- tT.logfc[[i]] %>% select("adj.P.Val","P.Value","t","B","logFC")
    }
  }
  return(list(tT.subset=tT.subset, tT.subset.logfc=tT.subset.logfc))
}


#### plot functions

gene_hist <- function(GSEpheno.clean, meta, tT.logfc, tT){
  if (is.null(tT)){
    hist(tT.logfc$adj.P.Val, col = "grey", border = "white", xlab = "Adjusted P-value",
         ylab = "Number of genes", main = "Adjusted P-value distribution")
  } else {
    hist(tT$adj.P.Val, col = "grey", border = "white", xlab = "Adjusted P-value", 
         ylab= "Number of genes", main = " Adjusted P-value distribution")
  }
}

gene_qqplot <- function(fit2){
  t.good <- which(!is.na(fit2$F)) 
  qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")
}

gene_volcano <- function(cts, tT.logfc, volpval){
  if(length(cts)==1){
    tT.logfc <- tT.logfc %>% 
      mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < volpval ~ "UP",
                                       logFC < 0 & adj.P.Val < volpval ~ "DOWN",
                                       TRUE ~ "NO"))
    p <- ggplot(data=tT.logfc, aes(x=logFC, y=-log10(P.Value), col=diffexpressed)) + 
      geom_point() + 
      theme_minimal() +
      ggtitle(paste(cts))  + 
      xlab("Log2 Fold Change") +
      ylab("-log10(P-value)") +
      theme(
        plot.title = element_text(hjust=0.5, size=20),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        plot.background = element_rect(fill = "white",
                                       colour = "white"))
    p
  } else {
    p <- NULL
    for(i in 1:length(cts)){
      tT.logfc[[i]] <- tT.logfc[[i]] %>% 
        mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < volpval ~ "UP",
                                         logFC < 0 & adj.P.Val < volpval ~ "DOWN",
                                         TRUE ~ "NO"))
      p[[i]] <- ggplot(data=tT.logfc[[i]], aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed)) + 
        geom_point() + 
        theme_minimal() +
        ggtitle(paste(cts[i])) + 
        xlab("Log2 Fold Change") +
        ylab("-log10(P-value)") +
        theme(
          plot.title = element_text(hjust=0.5, size=20),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18), 
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          plot.background = element_rect(fill = "white",
                                         colour = "white"))
    }
    p
  }
  
}


gene_md <- function(cts, tT.logfc, mdpval){
  if(length(cts)==1){
    tT.logfc <- tT.logfc %>% 
      mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < mdpval ~ "UP",
                                       logFC < 0 & adj.P.Val < mdpval ~ "DOWN",
                                       TRUE ~ "NO"))
    p <- ggplot(data=tT.logfc, aes(x=AveExpr, y=logFC, col=diffexpressed)) + 
      geom_point() + 
      theme_minimal() +
      ggtitle(paste(cts)) +
      xlab("Average log expression") + 
      ylab("Log2 Fold Change") +
      geom_hline(yintercept=0) +
      theme(
        plot.title = element_text(hjust=0.5, size=20),
        axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18), 
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        plot.background = element_rect(fill = "white",
                                       colour = "white"))
    p
  } else {
    p <- NULL
    for(i in 1:length(cts)){
      tT.logfc[[i]] <- tT.logfc[[i]] %>% 
        mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < mdpval ~ "UP",
                                         logFC < 0 & adj.P.Val < mdpval ~ "DOWN",
                                         TRUE ~ "NO"))
      p[[i]] <- ggplot(data=tT.logfc[[i]], aes(x=AveExpr, y=logFC, col=diffexpressed)) + 
        geom_point() + 
        theme_minimal() +
        ggtitle(paste(cts[i])) +
        xlab("Average log expression") + 
        ylab("Log2 Fold Change") +
        geom_hline(yintercept=0) +
        theme(
          plot.title = element_text(hjust=0.5, size=20),
          axis.title.x = element_text(size=16),
          axis.title.y = element_text(size=16),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18), 
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
          plot.background = element_rect(fill = "white",
                                         colour = "white"))
    }
    p
  }
  
}


gene_vol_table <- function(cts, tT.logfc, volpval){
  if(length(cts)==1){
    tT.logfc <- tT.logfc %>% 
      mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < volpval ~ "UP",
                                       logFC < 0 & adj.P.Val < volpval ~ "DOWN",
                                       TRUE ~ "NO"))
    tT.logfc <- tT.logfc %>% filter(diffexpressed != "NO")
    tT.logfc
  } else {
    for(i in 1:length(cts)){
      tT.logfc[[i]] <- tT.logfc[[i]] %>% 
        mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < volpval ~ "UP",
                                         logFC < 0 & adj.P.Val < volpval ~ "DOWN",
                                         TRUE ~ "NO"))
      tT.logfc[[i]] <- tT.logfc[[i]] %>% filter(diffexpressed != "NO")
    }
    tT.logfc
  }
}

gene_md_table <- function(cts, tT.logfc, mdpval){
  if(length(cts)==1){
    tT.logfc <- tT.logfc %>% 
      mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < mdpval ~ "UP",
                                       logFC < 0 & adj.P.Val < mdpval ~ "DOWN",
                                       TRUE ~ "NO"))
    tT.logfc <- tT.logfc %>% filter(diffexpressed != "NO")
    tT.logfc
  } else {
    for(i in 1:length(cts)){
      tT.logfc[[i]] <- tT.logfc[[i]] %>% 
        mutate(diffexpressed = case_when(logFC > 0 & adj.P.Val < mdpval ~ "UP",
                                         logFC < 0 & adj.P.Val < mdpval ~ "DOWN",
                                         TRUE ~ "NO"))
      tT.logfc[[i]] <- tT.logfc[[i]] %>% filter(diffexpressed != "NO")
    }
    tT.logfc
  }
}

gene_densities <- function(GSEassay, gs){
  plotDensities(GSEassay, group=gs, legend ="topleft")
}

gene_umap <- function(GSE.t.assay, gs, nneighs){
  ump <- umap(GSE.t.assay, n_neighbors = nneighs, random_state = 123)
  df <- data.frame(x=ump$layout[,1],
                      y=ump$layout[,2])
  ggplot(df, aes(x, y, colour=gs)) + 
    geom_point(size=1.5) +
    labs(colour="Group") +
    xlab("") +
    ylab("") +
    ggtitle(paste0("UMAP (neighbours: ", nneighs, ")")) +
    theme(
      plot.title = element_text(hjust=0.5, size=20),
      axis.title.x = element_text(size=16),
      axis.title.y = element_text(size=16),
      axis.text.x = element_text(size=18),
      axis.text.y = element_text(size=18), 
      legend.title = element_text(size=14),
      legend.text = element_text(size=12),
      plot.background = element_rect(fill = "white",
                                     colour = "white"))
}


gene_meanvar <- function(fit2, geoacc){
  plotSA(fit2, main=paste0(geoacc," Mean Variance Trend"), cex=0.75, cex.axis=1.5, cex.lab=1.5, cex.main=1.5)
}






gene_paths <- function(tT.logfc, tT, set){
  if (is.null(tT)){
    tT.row <- rownames_to_column(tT.logfc, "Gene.symbol")
    input_df <- tT.row[,c('Gene.symbol','adj.P.Val')]
    output_df <- tryCatch({run_pathfindR(input_df, p_val_threshold = 0.05, gene_sets=set, plot_enrichment_chart=FALSE)}, 
                          error= function(cond) {
        message("Here's the original error message:")
        #message(cond)
        return(NA)
      }, warning=function(w){
        message("Here's the original error message:")
        #message(w)
        return(NA)
        })
  } else {
    tT.row <- rownames_to_column(tT, "Gene.symbol")
    input_df <- tT.row[,c('Gene.symbol','adj.P.Val')]
    output_df <- tryCatch({run_pathfindR(input_df, p_val_threshold = 0.05, gene_sets=set, plot_enrichment_chart=FALSE)}, 
                          error= function(cond) {
      message("Here's the original error message:")
      #message(cond)
      return(NA)
    }, warning=function(w){
      message("Here's the original error message:")
      #message(w)
      return(NA)
    })
    }
  return(output_df)
}






library(ropls)

gene_paths_plot <- function(output_df, nterms){
    enrichment_chart(result_df = output_df, top_terms = nterms)
}


gene_plsda <- function(GSEassay, GSEpheno.clean, meta){
  sacSet <- Biobase::ExpressionSet(assayData = as.matrix(GSEassay), 
                                   phenoData = new("AnnotatedDataFrame", 
                                                   data = GSEpheno.clean))
  sacPlsda <- opls(sacSet, meta)
  sacSet <- getEset(sacPlsda)
  return(list(sacPlsda=sacPlsda, sacSet=sacSet))
}

gene_plsda_plot <- function(sacPlsda){
  if(nrow(sacPlsda@modelDF)>1){
  plot(sacPlsda, typeVc="x-score")    
  }else{
  plot(sacPlsda)    
  }
}









