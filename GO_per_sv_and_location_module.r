# Author: Anna Gelbe
# Date: 13.07.2023

# GO Term enrichment 
go_per_sv_and_location <- function(data_sets, allvar_df, svtypes, split_locations, location, dir_out, file_name) {

  for (c in 1:length(data_sets)) {
    for (sv in svtypes) {
      # filter for current data and sv type
      temp <- allvar_df %>%
        filter(DATA == data_sets[c]) %>%
        filter(SVTYPE %in% c(sv)) 

      # extract gene_names for enrichment
      genes_to_test <- temp$SYMBOL
      number_of_genes <- length(unique(genes_to_test))

      # Define output directory for the specific SV type
      if (sv == "DUP:TANDEM") {
        dir_plot <- paste(dir_out, data_sets[c], "GO_enrichment", "DUP_TANDEM", sep = "/")
      } else if (sv == 'DUP:INT') {
        dir_plot <- paste(dir_out, data_sets[c], "GO_enrichment", "DUP_INT", sep = "/")
      } else {
        dir_plot <- paste(dir_out, data_sets[c], "GO_enrichment", sv, sep = "/")
      }

      # Create the dir_plot if it doesn't exist
      if (!dir.exists(dir_plot)) {
        dir.create(dir_plot, recursive = TRUE)
      }

      # Perform enrichment analysis for Molecular Function (MF)
      if (number_of_genes != 0) {
        GO_results_MF <- enrichGO(gene = unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "MF")
      }

      # Generate MF plot
      if (number_of_genes == 0 || is.null(GO_results_MF) || nrow(GO_results_MF) == 0) {
        # No enriched genes found for MF
        if (sv == 'DUP:TANDEM') {
          pdf(file = paste(dir_plot, "/", "DUP_TANDEM_MF_none.pdf", sep = ""))
        } else if (sv == 'DUP:INT') {
          pdf(file = paste(dir_plot, "/", "DUP_INT_MF_none.pdf", sep = ""))
        } else {
          pdf(file = paste(dir_plot, "/", sv, "_MF_none.pdf", sep = ""))
        }
        par(mar = c(0, 0, 0, 0))
        plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
        text(x = 0.5, y = 0.5, paste("No enriched genes found for Molecular Functions in", sv, sep = ""), cex = 0.8)
        print(plot)
        dev.off()
      } else {
        # Enriched genes found for MF, generate barplot
        if (sv == 'DUP:TANDEM') {
          pdf(file = paste(dir_plot, "/", "DUP_TANDEM_MF.pdf", sep = ""))
        } else if (sv == 'DUP:INT') {
          pdf(file = paste(dir_plot, "/", "DUP_INT_MF.pdf", sep = ""))
        } else {
          pdf(file = paste(dir_plot, "/", sv, "_MF.pdf", sep = ""))
        }
        plot <- barplot(GO_results_MF,font.size = 8, showCategory = 20)+
            labs (title = paste(data_sets[c], "/", sv, ": GO Enrichment Molecular Function With ", number_of_genes, " Genes To Test", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"))
        print(plot)
        dev.off()
      }

  # BP
  # mit unique gene names
  if(number_of_genes != 0) #&& any(!is.na(temp$LOCEND - temp$LOCSTART)) && any(abs(temp$LOCEND - temp$LOCSTART) >= 10) && any(abs(temp$LOCEND - temp$LOCSTART) <= 500))
  {
    GO_results_BP <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "BP")
    print("all regions: 1. if case BP")
    }

  # dont plot if there are no results
  if (number_of_genes == 0 || is.null(GO_results_BP) || nrow(GO_results_BP) == 0) {
    print(paste( " 2. if case: BP", sep="" ))
    if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM",  "_BP_none.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT",  "_BP_none.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv,  "_BP_none.pdf", sep = ""))}
    par(mar = c(0, 0, 0, 0))  
    plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n",type = "n",xaxt = "n",yaxt = "n")
    text(x = 0.5, y = 0.5, paste("no enriched Genes found for \n Biological Processes for ", sv, sep = ""), cex = 0.8)
    print(plot)
    dev.off()
  }
  #otherwise create plot
  else 
  {
    if (sv == 'DUP:TANDEM') {pdf( file = paste(dir_plot,"/", "DUP_TANDEM",  "_BP.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT",  "_BP.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv,  "_BP.pdf", sep = ""))}
    plot <- barplot(GO_results_BP, font.size = 8, showCategory = 20)+ 
            labs(title = paste(data_sets[c], "/", sv, ": GO Enrichment Biological Process With ", number_of_genes, " Genes To Test", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"))
    print(plot)
    dev.off()
  }


  # CC
  if(number_of_genes != 0) #&& any(!is.na(temp$LOCEND - temp$LOCSTART)) && any(abs(temp$LOCEND - temp$LOCSTART) >= 10) && any(abs(temp$LOCEND - temp$LOCSTART) <= 500))
  {
    GO_results_CC <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "CC")
    print("all regions: 1. if case CC")
    }
  
  # dont plot if there are no results
  if (number_of_genes == 0 || is.null(GO_results_CC) || nrow(GO_results_CC) == 0) {
    print(paste( " 2. if case: CC", sep="" ))
    if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM",  "_CC_none.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT",  "_CC_none.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv,  "_CC_none.pdf", sep = ""))}
    par(mar = c(0, 0, 0, 0))  
    plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n",type = "n",xaxt = "n",yaxt = "n")
    text(x = 0.5, y = 0.5, paste("no enriched Genes found for \n Cellular Components for ", sv, sep = ""), cex = 0.8)
    print(plot)
    dev.off()
  }
  #otherwise create plot
  else {
  if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM",  "_CC.pdf", sep = ""))}
  else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT",  "_CC.pdf", sep = ""))}
  else{pdf(file = paste(dir_plot, "/", sv, "_CC.pdf", sep = ""))}
  plot <- barplot(GO_results_CC, font.size = 8, showCategory = 20) +
          labs(title = paste(data_sets[c], "/", sv,  ": GO Enrichment Cellular Components With ", number_of_genes, " Genes To Test", sep = "")) + 
          theme(plot.title = element_text(size = 10, face = "bold"))
  print(plot)
  dev.off()
  }

  # all in one
  if(number_of_genes != 0) #&& any(!is.na(temp$LOCEND - temp$LOCSTART)) && any(abs(temp$LOCEND - temp$LOCSTART) >= 10) && any(abs(temp$LOCEND - temp$LOCSTART) <= 500))
  {
    GO_results <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "ALL")
    print("all regions: 1. if case ALL")}
  
  if (number_of_genes == 0 || is.null(GO_results) || nrow(GO_results)== 0) {
    print(paste(" 2. if case: ALL ", sep="" ))
    if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM",  "_GO_none.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT",  "_GO_none.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv,  "_GO_none.pdf", sep = ""))}
    par(mar = c(0, 0, 0, 0))  
    plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n",type = "n",xaxt = "n",yaxt = "n")
    text(x = 0.5, y = 0.5, paste("no enriched Genes found for ", sv, sep = ""), cex = 0.8)
    print(plot)
    dev.off()
  }

  else{
  #plot
  if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM",  "_GO.pdf", sep = ""))}
  else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT",  "_GO.pdf", sep = ""))}
  else{pdf(file = paste(dir_plot, "/", sv, "_GO.pdf", sep = ""))}
  plot <- barplot(GO_results, font.size = 6, showCategory = 10, split = "ONTOLOGY") + 
          facet_grid(ONTOLOGY~., scale = "free") +
          labs(title = paste(data_sets[c], "/", sv, ": GO Enrichment With ", number_of_genes, " Genes To Test", sep = "")) + 
          theme(plot.title = element_text(size = 10, face = "bold"))
  print(plot)
  dev.off()

  }



  # ------------------------
  # GO Term per SV Type and Per Region
  if (split_locations == TRUE && length(location > 0))
  {
    for (loc in location)
  {
  
  # get only the right genes to test for GO Term depending on loc
  temp2 <- temp %>%
              filter(LOCATION == loc)

  # extract gene_names for enrichment
  genes_to_test <- temp2$SYMBOL
  number_of_genes <- length(unique(genes_to_test))
  print(number_of_genes)

  # ------
  # GO-TERM enrichment 
  # needs current to work 

  # directory
  if (sv == 'DUP:TANDEM') {dir_plot <- paste(dir_out, data_sets[c], "GO_enrichment", "DUP_TANDEM", loc , sep = "/")}
  else if (sv == 'DUP:INT') {dir_plot <- paste(dir_out, data_sets[c], "GO_enrichment", "DUP_INT", loc , sep = "/")}
  else {dir_plot <-paste(dir_out, data_sets[c], "GO_enrichment", sv , loc, sep = "/")}

  # Create the dir_plot if it doesn't exist
  if (!dir.exists(dir_plot)) {
      dir.create(dir_plot, recursive = TRUE)    
      }


  # GO
  # --
  # MF
  if((number_of_genes != 0) && any(!is.na(temp2$LOCEND - temp2$LOCSTART)) && any(abs(temp2$LOCEND - temp2$LOCSTART) >= 10) && any(abs(temp2$LOCEND - temp2$LOCSTART) <= 500)){ 
    GO_results_MF <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "MF")
    print(paste(loc, sv,  " 1. if case: MF", sep="" ))
  }

  # dont plot if there are no results
  if (number_of_genes == 0 || is.null(GO_results_MF) || nrow(GO_results_MF) == 0) {
    print(paste(loc, " 2. if case: MF", sep="" ))
    if (sv == 'DUP:TANDEM') {pdf(file =paste(dir_plot,"/", "DUP_TANDEM_", loc, "_MF_none.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file =paste(dir_plot,"/", "DUP_INT_", loc,  "_MF_none.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv, "_", loc,  "_MF_none.pdf", sep = ""))}
    par(mar = c(0, 0, 0, 0))  
    plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n",type = "n",xaxt = "n",yaxt = "n")
    text(x = 0.5, y = 0.5, paste("no enriched Genes found for \n Molecular Functions for ", sv, "/", loc,  sep = ""), cex = 0.8)
    print(plot)
    dev.off()
    }
  #otherwise create plot
  else {
    # barplot
    if (sv == 'DUP:TANDEM') {pdf(file =paste(dir_plot,"/", "DUP_TANDEM_",loc,  "_MF.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file =paste(dir_plot,"/", "DUP_INT_",loc, "_MF.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot, "/", sv, "_", loc, "_MF.pdf", sep = ""))}
    plot <- barplot(GO_results_MF,font.size = 8, showCategory = 20)+
            labs (title = paste(data_sets[c], "/", sv, "/", loc, ": GO Enrichment Molecular Function With ", number_of_genes, " Genes To Test", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"))
    print(plot)
    dev.off()
  }

  # BP
  # mit unique gene names
  if((number_of_genes != 0) && any(!is.na(temp2$LOCEND - temp2$LOCSTART)) && any(abs(temp2$LOCEND - temp2$LOCSTART) >= 10) && any(abs(temp2$LOCEND - temp2$LOCSTART) <= 500)){
    GO_results_BP <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "BP")
    print(paste(loc, sv,  " 1. if case: BP", sep="" ))
    }

  # dont plot if there are no results
  if (number_of_genes == 0 || is.null(GO_results_BP) || nrow(GO_results_BP) == 0) {
    print(paste(loc, " 2. if case: BP", sep="" ))
    if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM_",loc, "_BP_none.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT_",loc, "_BP_none.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv, "_", loc,  "_BP_none.pdf", sep = ""))}
    par(mar = c(0, 0, 0, 0))  
    plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n",type = "n",xaxt = "n",yaxt = "n")
    text(x = 0.5, y = 0.5, paste("no enriched Genes found for \n Biological Processes for ", sv, "/", loc, sep = ""), cex = 0.8)
    print(plot)
    dev.off()
  }
  #otherwise create plot
  else 
  {
    if (sv == 'DUP:TANDEM') {pdf( file = paste(dir_plot,"/", "DUP_TANDEM_",loc, "_BP.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT_",loc,  "_BP.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv, "_", loc,  "_BP.pdf", sep = ""))}
    plot <- barplot(GO_results_BP, font.size = 8, showCategory = 20)+ 
            labs(title = paste(data_sets[c], "/", sv, "/", loc, ": GO Enrichment Biological Process With ", number_of_genes, " Genes To Test", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"))
    print(plot)
    dev.off()
  }


  # CC
  if((number_of_genes != 0) && any(!is.na(temp2$LOCEND - temp2$LOCSTART)) && any(abs(temp2$LOCEND - temp2$LOCSTART) >= 10) && any(abs(temp2$LOCEND - temp2$LOCSTART) <= 500)){
    GO_results_CC <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "CC")
    print(paste(loc, sv, " 1. if case: CC", sep="" ))
    }
  
  # dont plot if there are no results
  if (number_of_genes == 0 || is.null(GO_results_CC) || nrow(GO_results_CC) == 0) {
    print(paste(loc, " 2. if case: CC", sep="" ))
    if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM_",loc, "_CC_none.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT_",loc,  "_CC_none.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv, "_", loc,  "_CC_none.pdf", sep = ""))}
    par(mar = c(0, 0, 0, 0))  
    plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n",type = "n",xaxt = "n",yaxt = "n")
    text(x = 0.5, y = 0.5, paste("no enriched Genes found for \n Cellular Components for ", sv, "/", loc, sep = ""), cex = 0.8)
    print(plot)
    dev.off()
  }
  #otherwise create plot
  else {
  if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM_",loc,  "_CC.pdf", sep = ""))}
  else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT_",loc,  "_CC.pdf", sep = ""))}
  else{pdf(file = paste(dir_plot, "/", sv, "_", loc, "_CC.pdf", sep = ""))}
  plot <- barplot(GO_results_CC, font.size = 8, showCategory = 20) +
          labs(title = paste(data_sets[c],"/", sv, "/", loc, ": GO Enrichment Cellular Components With ", number_of_genes, " Genes To Test", sep = "")) + 
          theme(plot.title = element_text(size = 10, face = "bold"))
  print(plot)
  dev.off()
  }

  # all in one
  if((number_of_genes != 0) && any(!is.na(temp2$LOCEND - temp2$LOCSTART)) && any(abs(temp2$LOCEND - temp2$LOCSTART) >= 10) && any(abs(temp2$LOCEND - temp2$LOCSTART) <= 500))
  {
    GO_results <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "ALL")
    print(paste(loc, sv,  " 1. if case: ALL", sep="" ))
  }
  
  if (number_of_genes == 0 || is.null(GO_results) || nrow(GO_results)== 0) {
    print(paste(loc, " 2. if case: ALL", sep="" ))
    if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM_",loc,  "_GO_none.pdf", sep = ""))}
    else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT_",loc, "_GO_none.pdf", sep = ""))}
    else{pdf(file = paste(dir_plot,"/", sv, "_", loc,  "_GO_none.pdf", sep = ""))}
    par(mar = c(0, 0, 0, 0))  
    plot <- plot(x = 0:1, y = 0:1, ann = F, bty = "n",type = "n",xaxt = "n",yaxt = "n")
    text(x = 0.5, y = 0.5, paste("no enriched Genes found for ", sv, "/", loc, sep = ""), cex = 0.8)
    print(plot)
    dev.off()
  }

  else{
  #plot
  if (sv == 'DUP:TANDEM') {pdf(file = paste(dir_plot,"/", "DUP_TANDEM_",loc,  "_GO.pdf", sep = ""))}
  else if (sv == 'DUP:INT') {pdf(file = paste(dir_plot,"/", "DUP_INT_",loc, "_GO.pdf", sep = ""))}
  else{pdf(file = paste(dir_plot, "/", sv,"_", loc,  "_GO.pdf", sep = ""))}
  plot <- barplot(GO_results, font.size = 6, showCategory = 10, split = "ONTOLOGY") + 
          facet_grid(ONTOLOGY~., scale = "free") +
          labs(title = paste(data_sets[c],"/", sv, "/", loc, ": GO Enrichment With ", number_of_genes, " Genes To Test", sep = "")) + 
          theme(plot.title = element_text(size = 10, face = "bold"))
  print(plot)
  dev.off()
  }

  }
  }
  else{print("no locations chosen")}



}

}
}