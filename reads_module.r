

supporting_reads <- function(data_sets, all_data, dir_out, file_name, color_palette) {
    # Create the fraction column
    all_data$fraction <- all_data$SUPPORT / all_data$gt_DP
    
    for (c in 1:length(data_sets)) {
        temp <- all_data %>%
            filter(data == data_sets[c])
        
        # Create the dir_out if it doesn't exist
        if (!dir.exists(dir_out[c])) {
            dir.create(dir_out[c], recursive = TRUE)
        }
        
        # Fraction - Boxplot
        pdf(file = paste(dir_out, file_name[1]))
        plot <- ggplot(temp, aes(x = SVTYPE, y = fraction, fill = SVTYPE)) +
            geom_boxplot() +
            scale_fill_manual(values = color_palette) +
            labs(title = paste(data_sets[c], ": Fraction of Reads Supporting Per Event - Boxplot", sep = ""), x = "SV Type", y = "Support Per Read Depth") +
            theme_classic()
        print(plot)
        dev.off()
        
        # Fraction with Chrom - Scatter Plot
        pdf(file = paste(dir_out, file_name[2]))
        plot <- ggplot(temp, aes(x = fraction, y = CHROM, color = SVTYPE)) +
            geom_point() +
            scale_color_manual(values = color_palette) +
            scale_size(range = c(5, 15)) +
            labs(title = paste(data_sets[c], ": Fraction of Reads Supporting Per Event with Chromosome - Scatter Plot", sep = ""), x = "Support Per Read Depth", y = "Chromosome", color = "SV Type") +
            theme_classic()
        print(plot)
        dev.off()
        
        # Support per Read depth - Scatter Plot
        pdf(file = paste(dir_out, file_name[3]))
        plot <- ggplot(temp, aes(x = gt_DP, y = SUPPORT, color = SVTYPE)) +
            geom_point() +
            geom_abline(intercept = 0, slope = 1, color = "black") +
            scale_color_manual(values = color_palette) +
            scale_size(range = c(5, 15)) +
            labs(title = paste(data_sets[c], ": Support Per Read Depth - Scatter Plot", sep = ""), x = "Read Depth", y = "Support", color = "SV Type") +
            theme_classic() +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0))  # Adjust the expand and limits arguments
        print(plot)
        dev.off()
    }
}
