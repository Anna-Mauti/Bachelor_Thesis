##---------------------------------------------------------------------------

# Frequencies of locations
# for locating variants 
freq_of_locations <- function(data_sets, allvar_df, location, dir_out, file_name)
{

temp_max <- c()
all_loc <- data.frame()

for (c in 1:length(data_sets))
{
 
  # for counting frequencies
  freq_loc <- c()

  # create empty df
  loc_df <- data.frame(location = character(length(location)),  freq_loc = numeric(length(location)), data = (character(length(location))))
  
  # get frequencies of locations
  for (loc in location )
  {
    freq_loc <- append (freq_loc, nrow(filter(allvar_df, LOCATION == loc & DATA == data_sets[c]))) 
  }

  #  collect frequencies of sv in loc df
  loc_df$freq_loc <- freq_loc
  loc_df$location <- location
  loc_df$data <- data_sets[c]

  # search for maximum to get coord_cartesian
  temp_max <- append(temp_max, max(freq_loc))
  all_loc <- (rbind(all_loc, loc_df))

}

# for y axis scaling
max <- max(temp_max)
# aufrunden für 1000er Zahl
maxi <- ceiling(max / 1000) * 1000


for ( c in 1:length(data_sets))
{ 
  # filter correct data
  temp_loc <- all_loc %>% 
                filter(data == data_sets[c])

  # plot frequencies 


  # Create the dir_out[c] if it doesn't exist
  if (!dir.exists(dir_out[c])) {
    dir.create(dir_out[c], recursive = TRUE)
  }

  pdf(file = paste(dir_out[c], file_name, sep = "/"))
  plot <- ggplot(temp_loc, aes(x = location, y = freq_loc, fill = "#139a43")) + 
              geom_bar(stat = "identity") + 
              guides(fill = "none") +
              labs(title = paste(data_sets[c], ": Gene Locations For All SV Types", sep = ""), x = "Location", y = "Frequency") +
              #here I scale the y axis automatically
              coord_cartesian(ylim = c(0, maxi))
  print(plot)
  dev.off()
}

}