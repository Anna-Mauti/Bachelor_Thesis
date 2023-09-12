# Author: Anna Gelbe
# Date: 13.07.2023


# GO-TERM enrichment 
# needs allvar_df to work 
general_go_term <- function(data_sets, allvar_df, dir_out, file_name)
{
for (c in 1:length(data_sets))
{
    

    temp <- allvar_df %>%
                filter(DATA == data_sets[c])
    # extract gene_names for enrichment
    genes_to_test <- temp$SYMBOL
    genes_to_test <- as.data.frame(genes_to_test)
    genes_to_test <- genes_to_test$genes_to_test

    # number of genes
    number_of_genes <- length(unique(genes_to_test))

    #
    # mit unique gene names
    GO_results_BP <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "BP")
    GO_results_CC <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "CC")
    GO_results_MF <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "MF")
    GO_results <- enrichGO(gene=unique(genes_to_test), keyType = "SYMBOL", OrgDb = "org.Hs.eg.db", ont = "ALL")

    # BP
  
    # Create the dir_out[c] if it doesn't exist
    if (!dir.exists(dir_out[c])) {
    dir.create(dir_out[c], recursive = TRUE)
    }

    #BP
    # plot 
    pdf(file = paste(dir_out[c], file_name[1], sep = "/"))
    plot <- barplot(GO_results_BP, font.size = 8, showCategory = 20) + 
            labs(title = paste(data_sets[c], ": GO Enrichment Biological Process With ", number_of_genes, " Genes", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"))
    print(plot)
    dev.off()
    
    
    # CC
    # plot 
    pdf(file = paste(dir_out[c],file_name[2], sep = "/"))
    plot <- barplot(GO_results_CC, font.size = 8, showCategory = 20) +
            labs(title = paste(data_sets[c], ": GO Enrichment Cellular Components With ", number_of_genes, " Genes", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"))
    print(plot)
    dev.off()

    # MF
    # plot 
    pdf(file = paste(dir_out[c], file_name[3], sep = "/"))
    plot <- barplot(GO_results_MF,font.size = 8, showCategory = 20) +
            labs(title = paste(data_sets[c], ": GO Enrichment Molecular Function With ", number_of_genes, " Genes", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"))
    print(plot)
    dev.off()

    # all in one
    pdf(file = paste(dir_out[c], file_name[4], sep = "/"))
    plot <- barplot(GO_results,font.size = 6, showCategory = 10, split = "ONTOLOGY") + 
            facet_grid(ONTOLOGY~., scale = "free") +
            labs(title = paste(data_sets[c], ": GO Enrichment With ", number_of_genes, " Genes", sep = "")) + 
            theme(plot.title = element_text(size = 10, face = "bold"), axis.text.x=element_text(size=5))
    print(plot)
    dev.off()

}
}