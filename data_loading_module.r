# Author: Anna Gelbe
# Date: 13.07.2023

load_data <- function(data_sets, dir_variations, support) {

    # Create an empty list to store the data frames
    data_list <- list()
    
    for (n in 1:length(data_sets)) {
        # Read VCF file
        vcf <- read.vcfR(dir_variations[n])
        
        # Convert vcfr object to tibble
        frame <- vcfR2tidy(
            vcf,
            info_only = FALSE,
            single_frame = TRUE,
            toss_INFO_column = FALSE
        )
        
        # Filter variants based on minimum support
        data <- filter(frame$dat, SUPPORT >= support)
        data$data <- data_sets[n]
        
        # Store the data frame in a list
        data_list[[n]] <- data
    }
    
    # Combine all the data frames in the list into a single data frame
    all_data <- do.call(rbind, data_list)
    
    return(all_data)
}


