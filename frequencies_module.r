# Author: Anna Gelbe
# Date: 13.07.2023


frequencies <- function(data_sets, all_data, svtypes, dir_out, file_name, color_palette) {

    # initialize vector to store maximum frequency values
    temp_max <- c()  

    # initialize empty dataframe to store frequencies
    all_freq <- data.frame()  
    
    # Calculate frequencies for each dataset
    for (c in 1:length(data_sets)) {

        # initialize vector to store frequencies for each variant type
        frequency_vec <- c()  

        # create dataframe to store frequencies and related information
        freq_df <- data.frame(svtypes = svtypes, freq = numeric(length(svtypes)), 
                              prozente = numeric(length(svtypes)), data = character(length(svtypes)))
        
        # Iterate through each variant type
        for (sv in svtypes) {
            
            # calculate frequency by filtering data based on dataset and variant type conditions
            frequency_vec <- append(frequency_vec, nrow(filter(all_data, data == data_sets[c] & SVTYPE == sv))) 
        }

        # assign calculated frequencies to dataframe column
        freq_df$freq <- frequency_vec 

        # calculate percentages 
        freq_df$prozente <- prop.table(freq_df$freq) * 100  

        # store dataset name
        freq_df$data <- data_sets[c]  

        # find maximum frequency for scaling plots
        temp_max <- append(temp_max, max(frequency_vec))   
        
        # rbind frequency dataframe to overall dataframe
        all_freq <- rbind(all_freq, freq_df) 
    }

    # find the maximum frequency value
    max <- max(temp_max)  

    # round up to the nearest thousand
    max <- ceiling(max / 1000) * 1000  
    
    # Generate plots for each dataset
    for (c in 1:length(data_sets)) {

        # filter the frequency dataframe for the current dataset
        temp_freq <- all_freq %>% filter(data ==  data_sets[c])  

        # create output directory if it doesn't exist
        if (!dir.exists(dir_out[c])) {
            dir.create(dir_out[c], recursive = TRUE)  
        }

        # save the plot as a PDF file
        pdf(file = paste(dir_out[c], file_name, sep = "/"))  
        plot <- ggplot(temp_freq, aes(x = svtypes, y = freq, fill = svtypes)) +
                geom_bar(stat = "identity") +
                scale_fill_manual(values = color_palette) +
                labs(title = paste(data_sets[c], ": Frequencies Of Structural Variant Types", sep = ""),
                     x = "SV Types", y = "Frequency", fill = "SV Types") +
                coord_cartesian(ylim = c(0, max))  # adjust y-axis limits
        print(plot)
        dev.off() 
    }
}
