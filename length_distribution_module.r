#-----------------------------------------------------------------------------
## LENGTH DISTRIBUTION

length_distribution <- function(data_sets, all_data, dir_out, file_name, color_palette)
{
    if (nrow(all_data) == 0) {
    stop("Input data is empty. Please provide a non-empty data frame.")
  }

else{

# check for INVs because after SVIM they have NAs at SVLEN even though they have actual lengths
for(n in 1:nrow(all_data))
{
    if(all_data$SVTYPE[n] == "INV"){all_data$SVLEN[n] = all_data$END[n] - all_data$POS[n]}
    else{}
}

# add column with only absolute values of SVLEN
all_data$abs_SVLEN <- abs(all_data$SVLEN)

# add column with logs
all_data$log_SVLEN <- log2(all_data$abs_SVLEN)

for (c in 1:length(data_sets))
{
    temp <- all_data %>%
        filter(data == data_sets[c])


    # Create the dir_out if it doesn't exist
    if (!dir.exists(dir_out[c])) {
    dir.create(dir_out[c], recursive = TRUE)
    }


    # scatterplot: zeigt SV Type und frequency, mit Info auf welchem chr 
    pdf(file = paste(dir_out[c], file_name[1], sep = "/"))
    plot <- ggplot(temp, aes(x = SVLEN, y = CHROM, color = SVTYPE)) +
    geom_point() +
    scale_size(range = c(5, 15)) +
    scale_color_manual(values = color_palette) + 
    labs(title = paste(data_sets[c], ": SV Length And Chromosome Information", sep = ""), x = "SV Length", y = "Chromosome", color = "Type") +
    theme_classic()
    print(plot)
    dev.off()

    # boxplot log
    pdf(file = paste(dir_out[c], file_name[2], sep = "/"))
    plot <- ggplot(temp, aes(x = SVTYPE, y = log_SVLEN,  fill = SVTYPE)) +
    geom_boxplot() +
    scale_fill_manual(values = color_palette) +
    labs(title = paste(data_sets[c],": Log Length Distribution SV", sep = ""), x = "SV Type", y = "Log Length", fill= "SV Type") +
    theme_classic()
    print(plot)
    dev.off()

}
}
}

