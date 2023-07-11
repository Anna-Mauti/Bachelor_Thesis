

# FREQUENCIES 

#create tibble with frequencies and plot for every data
frequencies <- function(data_sets, all_data, svtypes, dir_out, file_name, color_palette){
    
temp_max <- c()

all_freq <- data.frame()

for (c in 1:length(data_sets))
{
    frequency_vec <- c()

    freq_df <- data.frame(svtypes = svtypes, freq = numeric(length(svtypes)), prozente = numeric(length(svtypes)), data = (character(length(svtypes))))

    for (sv in svtypes)
    {
        frequency_vec <- append(frequency_vec, nrow(filter(all_data, data == data_sets[c] & SVTYPE == sv))) 
    }
        
    freq_df$freq <- frequency_vec
    freq_df$prozente <- prop.table(freq_df$freq)*100
    freq_df$data <- data_sets[c]

    # search for maximum to get coord_cartesian
    temp_max <- append(temp_max, max(frequency_vec))
    all_freq <- (rbind(all_freq, freq_df))

}

max <- max(temp_max)
# aufrunden für 1000er Zahl 
max <- ceiling(max / 1000) * 1000

for ( c in 1:length(data_sets))
{
    # filter correct data
    temp_freq <- all_freq %>% 
                filter(data ==  data_sets[c])
    
    # plot frequencies 

    # Create the dir_out if it doesn't exist
    if (!dir.exists(dir_out[c])) {
    dir.create(dir_out[c], recursive = TRUE) 
    }


    # plot frequencies
    pdf(file = paste(dir_out[c], file_name, sep = "/"))
    plot <- ggplot(temp_freq, aes(x = svtypes, y = freq, fill = svtypes))+
        geom_bar(stat = "identity") +
        scale_fill_manual(values = color_palette) + 
        labs(title = paste( data_sets[c], ": Frequencies Of Structural Variant Types", sep = ""), x = "SV Types", y = "Frequency", fill = "SV Types" ) +
        #here the y axis get scaled automatically
        coord_cartesian(ylim = c(0, max))
    print(plot)
    dev.off()

      
}

}





    