# Author: Anna Gelbe
# Date: 13.07.2023

heatmap <- function(data_sets, allvar_df, svtypes, dir_out, file_name) {
        
# Filter data for the specified data sets
filter_data <- allvar_df %>%
         filter(DATA %in% data_sets)

  # Iterate over each SV type
  for (sv in 1:length(svtypes)) {
        
    # Filter data for the current SV type
    filtered_heat <- filter_data %>% 
      filter(SVTYPE %in% c(svtypes[sv])) 

    # Filter genes with at least 2 different data sets
    filtered_heat <- filtered_heat %>%
      group_by(SYMBOL) %>%
      filter(n_distinct(DATA) >= 2)

    # Convert filtered data frame to matrix with gene count
    matrix <- xtabs(~DATA + SYMBOL, data = filtered_heat)

    # Create the output directory if it doesn't exist
    if (!dir.exists(dir_out[sv])) {
      dir.create(dir_out[sv], recursive = TRUE)    
    }

    # Generate and save the heatmap as a PDF file
    pdf(file = paste(dir_out[sv], file_name, sep = "/"))
    plot <- Heatmap(matrix, name = paste(svtypes[sv], " Heatmap", sep = ""),
                    col = colorRamp2(c(0, 5), c("yellow", "red")),
                    na_col = "black",
                    row_names_side = "left",
                    column_names_gp = gpar(fontsize = 8),
                    row_km = 0.5,
                    column_km = 1)
    print(plot)
    dev.off()
  }
}