


# --------------
# stacked svtypes and locations plot 

freq_of_locations_per_sv <- function(data_sets, allvar_df, location, svtypes, dir_out, file_name)
{
temp_max <- c()
all_loc <- data.frame()


for (c in 1:length(data_sets)) 
{

for (sv in svtypes)  {

  
  # collect frequencies of sv in locations
  freq_loc <- c()

  # create empty df
  loc_df <- data.frame(location = character(length(location)),  freq_loc = numeric(length(location)), svtype =  (character(length(location))), data = (character(length(location))))
  
  # get frequencies of locations
  for (loc in location )
  {
    foo <- allvar_df %>% 
            filter(DATA == data_sets[c]) %>%
            filter(SVTYPE %in% c(sv)) %>%
            filter( LOCATION == loc)
    freq_loc <- append (freq_loc, nrow(foo)) 
  }

  loc_df$freq_loc <- freq_loc
  loc_df$location <- location
  loc_df$data <-data_sets[c] 
  loc_df$svtype <- sv 

  all_loc <- (rbind(all_loc, loc_df))

  # search for maximum to get coord_cartesian
  temp_max <- append(temp_max, max(freq_loc))
  

}
}


# for y axis scaling
max <- max(temp_max)
# aufrunden für 1000er Zahl 
maxi <- ceiling(max / 1000) * 1000

for ( c in 1:length(data_sets))
{

# filter correct cell_line
temp_loc <- all_loc %>% 
              filter(data == data_sets[c])


# Create the dir_out[c] if it doesn't exist
if (!dir.exists(dir_out[c])) {
    dir.create(dir_out[c], recursive = TRUE)    
    }


#einmal nicht normalisieren 
# Stacked
pdf(file = paste(dir_out[c], file_name[1], sep = "/"))
plot <- ggplot(temp_loc, aes(fill = location, y = freq_loc, x = svtype)) + 
    geom_bar(position="stack", stat="identity") +
    labs(title = paste(data_sets[c], ": SV Types And Their Locations", sep = ""), x = "SV Type", y = "Frequency", fill = "Location" )+
    #here I scale the y axis automatically
    coord_cartesian(ylim = c(0, maxi))
print(plot)
dev.off()

#einmal normalisieren auf 100 % 
# Stacked + percent
pdf(file = paste(dir_out[c], file_name[2], sep = "/"))
plot <- ggplot(temp_loc, aes(fill = location, y = freq_loc, x = svtype))+ 
    geom_bar(position="fill", stat="identity") +
    labs(title = paste(data_sets[c], ": SV Types And Their Locations (Normalized)", sep = ""), x = "SV Type", y = "Frequency", fill = "Location" )
print(plot)
dev.off()

}

}
  
