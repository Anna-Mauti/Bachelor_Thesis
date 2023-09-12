# Author: Anna Gelbe
# Date: 13.07.2023

# CIRCOS PLOT
circos <- function(data_sets, all_data, dir_out, file_name, color_palette){

for (c in 1:length(data_sets))
{
    # filter for DEL, INS
    circos_data <- all_data %>% 
                   filter(SVTYPE %in% c("DEL", "INS")) %>%
                   filter(data == data_sets[c]) %>%
                   dplyr::select(CHROM, POS, END, SVLEN, SVTYPE)

    del_value <-data.frame(value = filter(circos_data, SVTYPE == 'DEL')$SVLEN)
    ins_value <-data.frame(value = filter(circos_data, SVTYPE == 'INS')$SVLEN)

    # filter for BND
    bnd_data <- all_data %>%
                   filter( data == data_sets[c]) %>%
                   filter(SVTYPE %in% c("BND"))

    # infos about start/end/chr
    chrA <- c()
    startA <- c()
    endA <- c()
    chrB <- c()
    startB <- c()
    endB <- c()

    # Creating an empty tibble
    translocations <- tibble(
    chrA = character(),
    startA = numeric(),
    endA = numeric(),
    chrB = character(),
    startB = numeric(),
    endB= numeric()
    )

    # save infos in vectors
    for (n in 1:(dim(bnd_data)[1]))
    {
    chrA <- append(chrA, bnd_data$CHROM[n]) 
    startA <- append(startA, bnd_data$POS[n])
    endA <- append(endA, bnd_data$POS[n])

    chrB <- append(chrB, gsub(".*?(chr[0-9XY]+).*", "\\1", bnd_data$ALT[n]))
    startB <- append(startB, gsub(".*:(\\d+).*", "\\1", bnd_data$ALT[n]))
    endB <- append(endB, gsub(".*:(\\d+).*", "\\1", bnd_data$ALT[n]))
    }

    # save in one tibble
    translocations <- tibble(chrA, startA, endA, chrB, startB, endB)
    print(paste(data_sets[c],nrow(translocations)))
    print(paste(data_sets[c],nrow(del_value)))
    print(paste(data_sets[c],nrow(ins_value)))

    # PLOT

    # Create the dir_out if it doesn't exist
    if (!dir.exists(dir_out[c])) {
    dir.create(dir_out[c], recursive = TRUE)
    }

    pdf(file = paste(dir_out[c], file_name, sep = "/"))

    # function to plot circos 
    circos_plot <- function () 
    {
        
    
    # set up the generic graphical parameters:
    par(mar = c(1, 1, 1, 1))
    circos.par("start.degree" = 90)
    circos.par("track.height" = 0.05)
    circos.par("canvas.xlim" = c(-1.3, 1.3), "canvas.ylim" = c(-1.3, 1.3))

    # draw hg19 reference ideograms:
    circos.initializeWithIdeogram(species = "hg19")


    # DEL
    circos.genomicTrackPlotRegion(filter(circos_data, SVTYPE %in% c("DEL")),stack=TRUE, panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value=del_value, cex = 0.001, pch = 9,col= color_palette[2] , ...)
    })

    # INS
    circos.genomicTrackPlotRegion(filter(circos_data, SVTYPE %in% c("INS")),stack=TRUE, panel.fun = function(region, value, ...) {
        circos.genomicPoints(region, value=ins_value, cex = 0.001, pch = 9,col= color_palette[5] , ...)
    })

    # add translocations
    # Iterate over the dataframe and add links
    for (i in 1:nrow(translocations) ) {
    circos.link(
        sector.index1 = translocations$chrA[i],
        point1 = as.numeric(translocations$startA[i]),
        sector.index2 = translocations$chrB[i],
        point2 = as.numeric(translocations$startB[i]),
        col = rand_color(nrow(translocations))
    )
    }
    # add title
    title(paste(data_sets[c], ": Circos For Deletions, Insertions and Translocations", sep = ""), cex = 0.8)

    dev.off()

    circos.clear()
    }

    # use function to plot
    circos_plot()


}

}