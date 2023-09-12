# Author: Anna Gelbe
# Date: 13.07.2023


# Visualization of results after variant calling with SVIM
# MAIN




# *************************************************
# Install needed packages in personal library
#
# This installs every needed package for the entire Framework in the specified r_lib directory.
# 
# Arguments:
# - r_lib: directory path for installing packages
#
# Output:
# - packages are installed

# Define the Directory for installing packages 
r_lib <- "/project/meissner_training/Anna/Bachelor/R_libs"

# Load install function
source(paste(dir_modules, "install_packages_module.r", sep = "/"))

# Call the install function
install_packages(r_lib)




# *************************************************
# Libraries

# Set the path to custom R libraries
r_lib <- "/project/meissner_training/Anna/Bachelor/R_libs"

# Load the vcfR package from the custom library
library(vcfR, lib.loc = r_lib)
library(memuse, lib.loc = r_lib)

# Load necessary packages for plots and data manipulation
library(ggplot2)
library(tidyr)
library(dplyr)

# Load the circlize package for circular plots
library(circlize, lib = r_lib)

# Load packages for gene annotation
library(GenomicRanges, lib.loc = r_lib)
library(VariantAnnotation, lib.loc = r_lib)
library(annotate, lib.loc = r_lib)
library(NLP, lib.loc = r_lib)

# Load packages for GO analysis
library(clusterProfiler, lib.loc = r_lib)
library(AnnotationDbi, lib.loc = r_lib)
library(TxDb.Hsapiens.UCSC.hg19.knownGene, lib.loc = r_lib)
library(org.Hs.eg.db, lib.loc = r_lib)
library(annotate, lib.loc = r_lib)

# Load packages for circos plot with BND
library(StructuralVariantAnnotation, lib = r_lib)

# Load packages for heatmap
library(ComplexHeatmap, lib = r_lib)

# *************************************************
# Definition of Global Variables

# Define Directory for sourcing modules
dir_modules <- "/project/meissner_training/Anna/Bachelor/analysis/Bachelor_Thesis"

# Define results Directory
dir_results <- "/project/meissner_training/Anna/Bachelor/analysis/results"

# *************************************************
# Load Functions


# Load the data loading function
source(paste(dir_modules, "data_loading_module.r", sep = "/"))

# Load frequencies function
source(paste(dir_modules, "frequencies_module.r", sep = "/"))

# Load length distribution function
source(paste(dir_modules, "length_distribution_module.r", sep = "/"))

# Load fraction of reads supporting per event function
source(paste(dir_modules, "normalized_cov_module.r", sep = "/"))

# Load circos plot function
source(paste(dir_modules, "circos_module.r", sep = "/"))

# Load gene annotation function
source(paste(dir_modules, "gene_annotation_module.r", sep = "/"))

# Load frequencies of locations function
source(paste(dir_modules, "frequencies_of_locations_module.r", sep = "/"))

# Load frequencies per location and sv type function
source(paste(dir_modules, "frequencies_of_locations_per_sv_module.r", sep = "/"))

# Load general GO Term function
source(paste(dir_modules, "general_GO_term_module.r", sep = "/"))

# Load the GO term enrichment analysis for SV types and locations function
source(paste(dir_modules, "GO_per_sv_and_location_module.r", sep = "/"))

# Load the heatmap function
source(paste(dir_modules, "heatmap_module.r", sep = "/"))




# *************************************************
# load_data: Load data from VCF files
# 
# This function reads the VCF file for each data set, filters the variants based on the minimum support value
# and combines the filtered variants into a single data frame called all_data.
# 
# Arguments:
# - data_sets: Vector of data set names
# - dir_variations: Directory path of VCF files for each data set
# - support: Minimum support value for filtering variants
#
# Output:
# - all_data: Combined data frame containing the filtered variants from all data sets

# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the directories of the input VCF files
dir_variations <- c()

# Save the different input directories in the variable 'dir_variations'
for (n in 1:length(data_sets))
{
    dir_variations <- append(dir_variations, paste("/project/meissner_training/Anna/Bachelor/variations/", data_sets[n], "/variants.vcf", sep = ""))
}

# Define a variable to filter the VCF file for minimum support
support = 10 

# Load the variant data using the load_data function and save the result to all_data
all_data <- load_data(data_sets, dir_variations, support)



# *************************************************
# frequencies: Calculate Frequencies and Generate Plots
#
# This function calculates frequencies of structural variant types (svtypes) for multiple datasets,
# and generates individual plots for each dataset displaying the frequencies. The function takes 
# the input data, svtypes, output directory, and file names as arguments.
#
# Arguments:
#   - all_data: The input data containing information about structural variants. It should be a data 
#               frame with columns 'data', 'SVTYPE', and other relevant columns generated with the function: load_data
#   - svtypes: A character vector specifying the structural variant types to consider for frequency calculation.
#   - dir_out: A character vector specifying the output directory path where the generated plots will be saved.
#   - file_name: A character vector specifying the file name of the generated plot.
#   - color_palette: A color vector to ensure consistency in colors for each SV type across all modules, enabling easy visual comparison.
#
# Output:
#   The function generate a individual PDF plot for each dataset, displaying the frequencies of structural 
#   variant types. The plots are saved in the specified output directory with the specified file names.
#
# Example Usage:
#    frequencies(data_sets, all_data, svtypes, dir_out, file_name)
#
# Notes:
#   - Make sure to provide the correct input data, svtypes, output directory, and file names.
#   - The 'all_data' input should contain the necessary columns 'data' and 'SVTYPE' for filtering and calculations.
#   - Ensure that the required libraries (ggplot2, dplyr) are loaded before using this function.


# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the SV types
svtypes <- c('DEL', 'INS', 'INV','DUP:TANDEM', 'DUP:INT', 'BND')

# Define the output directory
dir_out <- c()
for (n in 1:length(data_sets))
{
    dir_out[n] <- paste(dir_results, data_sets[n], "frequencies", sep = "/")
}

# Define the output file name
file_name <- "frequencies.pdf"

# Color palette for SV types
color_palette <- c('BND' = "#ff91ed", 'DEL' = "#81ebe4", 'DUP:INT' = "#ffc738", 'DUP:TANDEM' = "#b4fa82", 'INS' = "#fa7f7f", 'INV' = "#c18ff7")

# Call the frequencies function
frequencies(data_sets, all_data, svtypes, dir_out, file_name, color_palette)

     


# *************************************************
# length_distribution: Calculate Length Distribution and Generate Plots
#
# This function calculates the length distribution of structural variants (SVs) for multiple datasets, and
# generates two types of plots for each dataset. The function takes the data set names,  output
# directory, file names and a color palette as arguments.
#
# Arguments:
#   - all_data: The input data containing information about structural variants. It should be a data frame
#               or tibble with relevant columns such as 'data', 'SVTYPE', 'POS', 'END', and 'SVLEN'.
#   - data_sets: A character vector specifying the dataset names to consider for length distribution calculation
#                and plotting.
#   - dir_out: A character vector specifying the output directory path where the generated plots will be saved.
#   - file_name: A character vector specifying the file names of the generated plots.
#   - color_palette: A color vector to ensure consistency in colors for each SV type across all modules, enabling easy visual comparison.
#
# Output:
#   The function generates two types of plots for each dataset:
#   - Scatterplot: Displays the SV Length against Chromosome with different SV types represented by colors.
#   - Boxplot: Displays the log-transformed length distribution of SVs for different SV types.
#   The plots are saved in the specified output directory with the specified file names.
#
# Example Usage:
#   length_distribution(all_data, data_sets, dir_out, file_name, color_palette)
#
# Notes:
#   - Make sure to provide the correct input data, data sets, output directory, file names and color palette.
#   - The 'all_data' input should contain the necessary columns 'data', 'SVTYPE', 'POS', 'END', and 'SVLEN'.
#   - The length of 'data_sets', 'dir_out' should match each other.
#   - Ensure that the required libraries (ggplot2, dplyr) are loaded before using this function.


# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the output directories
dir_out <- c()
for (n in 1:length(data_sets))
{
    dir_out[n] <- paste(dir_results, data_sets[n], "length_distribution", sep = "/")
}

# Define the output file names
file_name <- c("length_distribution_scatter.pdf", "length_distribution_box_log.pdf")

# Color palette for SV types
color_palette <- c('BND' = "#ff91ed", 'DEL' = "#81ebe4", 'DUP:INT' = "#ffc738", 'DUP:TANDEM' = "#b4fa82", 'INS' = "#fa7f7f", 'INV' = "#c18ff7")

# Call the length distribution function
length_distribution(data_sets, all_data, dir_out, file_name, color_palette)


# *************************************************
# Fraction of reads supporting per event
# normalized_coverage
# Arguments:
#   - data_sets: Vector of data set names
#   - all_data: Data frame containing the variant data
#   - dir_out: vector of Output directory paths
#   - file_name: Vector of output file names
#   - color_palette: A color vector to ensure consistency in colors for each SV type across all modules, enabling easy visual comparison
#   
# Output:
#    PDF files:
#      1. Fraction - Boxplot: A boxplot showing the distribution of support per read depth for different SV types.
#         The color of the boxplots corresponds to the SV types defined in the color_palette.
#      2. Fraction with Chrom - Scatter Plot: A scatter plot showing the relationship between the fraction of reads supporting each event and the chromosome.
#         Each data point represents an event, with the color indicating the SV type.
#      3. Support per Read depth - Scatter Plot: A scatter plot showing the relationship between the support and read depth for each event.
#         The data points are colored based on the SV type, and a diagonal line represents a 1:1 ratio between support and read depth.
#    The PDF files are saved in the dir_out directory, which should be provided as an argument to the function.


# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the output directories
dir_out <- c()
for (n in 1:length(data_sets))
{
    dir_out[n] <- paste(dir_results, data_sets[n], "reads", sep = "/")
}

# Define the output file names
file_name <- c("fraction_boxplot.pdf","fraction_with_chr_scatter.pdf","support_per_read_depth_scatter.pdf")

# Color palette for SV types
color_palette <- c('BND' = "#ff91ed", 'DEL' = "#81ebe4", 'DUP:INT' = "#ffc738", 'DUP:TANDEM' = "#b4fa82", 'INS' = "#fa7f7f", 'INV' = "#c18ff7")

# Call the normalized coverage function
normalized_coverage(data_sets, all_data, dir_out, file_name, color_palette)



# *************************************************
# Circos plot for deletions, insertions, and translocations
#
# This function generate a Circos plot for deletions, insertions, and translocations. 
# 
# Arguments:
#   - data_sets: Vector of data set names
#   - all_data: Data frame containing the variant data
#   - dir_out: vector of Output directory paths
#   - file_name: Vector of output file names
#   - color_palette: A color vector to ensure consistency in colors for each SV type across all modules, enabling easy visual comparison
# 
# Output:
#   - set of Circos plots saved as PDF files in the specified directory
#
# Notes:
#   - the color of DEL and INS is not defined in a legend, but you know it from the color_palette 


# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the output directories
dir_out <- c()
for (n in 1:length(data_sets))
{
    dir_out[n] <- paste(dir_results, data_sets[n], "circos", sep = "/")
}

# Define output file name
file_name <- c("circos.pdf")

# Color palette for SV types
color_palette <- c('BND' = "#ff91ed", 'DEL' = "#81ebe4", 'DUP:INT' = "#ffc738", 'DUP:TANDEM' = "#b4fa82", 'INS' = "#fa7f7f", 'INV' = "#c18ff7")

# Call the circos function
circos(data_sets, all_data, dir_out, file_name, color_palette)



# *************************************************
# Gene Annotation
# This function performs gene annotation for structural variants (SVs) in the input data sets.
# Arguments:
#   - data_sets: Vector of data set names
#   - all_data: Data frame containing the variant data
#   - svtypes: A vector of SV types to consider for annotation.
#
# Output:
#   - allvar_df: a data frame containing the annotated SVs with gene IDs, gene symbols, SV types, and data set information.


# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the SV types
svtypes <- c('DEL', 'INS', 'INV','DUP:TANDEM', 'DUP:INT', 'BND')

# Call the gene annotation function and save result to allvar_df
allvar_df <- gene_annotation(data_sets, all_data, svtypes )



# *************************************************
# Frequencies of locations 
# This function calculates the frequencies of locations for the variants in the given data sets.
# It generates bar plots of location frequencies for each data set and saves them as PDF files.
#
# Arguments:
#  - data_sets: A vector specifying the data sets to analyze.
#  - allvar_df: The dataframe containing the variant data with location annotations.
#  - location: A vector specifying the locations to analyze.
#  - dir_out: A vector specifying the output directories for each data set.
#  - file_name: The output file name for the bar plots.
#
# Output:
#  - Bar plots of location frequencies for each data set, saved as PDF files in the specified output directories.
#
# Requirements:
#  - The `gene_annotation` function should be executed prior to this function to annotate the gene locations in the `allvar_df` dataframe.


# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the locations
location <- c( 'coding','fiveUTR', 'threeUTR', 'intron', 'intergenic', 'spliceSite', 'promoter')

# Define the output directory 
dir_out <- c()
for (n in 1:length(data_sets))
{
    dir_out[n] <- paste(dir_results, data_sets[n], "locations", sep = "/")
}
# Define output file name
file_name <- c("frequencies_of_locations.pdf")

# Call the freq_of_location function
freq_of_locations(data_sets, allvar_df, location, dir_out, file_name)



# *************************************************
# Frequencies of locations per SV Type
# This function calculates the frequencies of locations per SV Type for the variants in the given data sets.
# It generates stacked bar plots to visualize the frequencies of SV (structural variant) types in different locations.
#
# Arguments:
#  - data_sets: A vector specifying the data sets to analyze.
#  - allvar_df: The dataframe containing the variant data with location annotations.
#  - location: A vector specifying the locations to analyze.
#  - svtypes: A vector specifying the svtypes to analyze.
#  - dir_out: A vector specifying the output directories for each data set.
#  - file_name: The output file name for the bar plots.
#
# Output:
#  - Two plots are generated for each data set: one without normalization and one with normalization.
#  - The plots are saved as PDF files in the specified output directories (dir_out) with the provided file names (file_name).
#  
# Requirements:
#  - The `gene_annotation` function should be executed prior to this function to annotate the gene locations in the `allvar_df` dataframe.


# Define the names of the data sets
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the locations
location <- c( 'coding','fiveUTR', 'threeUTR', 'intron', 'intergenic', 'spliceSite', 'promoter')

# Define the svtypes 
svtypes <- c('DEL', 'INS', 'INV','DUP:TANDEM', 'DUP:INT', 'BND')

# Define the output directory 
dir_out <- c()
for (n in 1:length(data_sets))
{
    dir_out[n] <- paste(dir_results, data_sets[n], "locations", sep = "/")
}

# Define the output file name
file_name <- c("freq_locations_per_sv_stacked.pdf", "freq_locations_per_sv_stacked_norm.pdf")

# Call the frequencies of locations per sv type function
freq_of_locations_per_sv (data_sets, allvar_df, location, svtypes, dir_out, file_name)



# *************************************************
# GO-TERM enrichment 
# This function performs Gene Ontology (GO) term enrichment analysis for a given set of data sets using the allvar_df dataframe.
# The allvar_df dataframe must contain SV data with gene annotations.
#
# Arguments:
# - data_sets: A vector specifying the data sets to analyze.
# - allvar_df: The dataframe containing SV data with gene annotations.
# - dir_out: The output directory where the results will be saved.
# - file_name: The name of the output file.
#
# Output:
# - PDF files: The function generates PDF files for each data set, containing bar plots of GO enrichment results for Molecular Function (MF),
#   Biological Process (BP) and Cellular Component (CC) and one bar plot for all 3.
#
# Requirements:
# - Gene annotation must be performed on the SV data and stored in the allvar_df dataframe.


# Define the data sets to analyze
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the output directory for each data set
dir_out <- c()
for (n in 1:length(data_sets)) {
  dir_out[n] <- paste(dir_results, data_sets[n], "GO_enrichment/all_sv", sep = "/")
}

# Define the output file names
file_name <- c("GO_BP.pdf", "GO_CC.pdf", "GO_MF.pdf", "GO.pdf")

# Perform GO term enrichment analysis for each data set using the general_go_term function
general_go_term(data_sets, allvar_df, dir_out, file_name)





# *************************************************
# go_per_sv_and_location function performs Gene Ontology (GO) term enrichment analysis
# for specified SV types and locations within the provided data sets.
#
# Arguments:
#   - data_sets: A character vector specifying the data sets to analyze.
#   - allvar_df: A data frame containing the variant data.
#   - svtypes: A character vector specifying the SV types to analyze.
#   - split_locations: A logical value indicating whether to look at split locations.
#   - location: A character vector specifying the locations to analyze.
#   - dir_out: A character strting specifying the output directory.
#   - file_name: A character string specifying the file name for the output plots.
#
# Output:
#   - The files for the GO Term per Svtype are saved in the output directory with the conevntion:<dir_out>/<data_set>/GO_enrichment>/all_sv/
#   - The files for the split locations are saved in the output directory with the conevntion: <dir_out>/<data_set>/GO_enrichment>/<Svtype>/<location> 
#     and with the naming convention: <SVtype>_<location>_<ontology>.pdf.
#   - Each file represents the GO term enrichment bar plot for a specific data set, SV type, location, and ontology.
#   - Empty plots if no results are found
# 
# Requirements:
# - The gene annotation must be performed prior to running this function.
# - The 'allvar_df' dataframe containing SV data with gene annotations must be available.
#
# Note: The function requires the org.Hs.eg.db package for gene annotation and the enrichGO function
# from the clusterProfiler package for performing GO term enrichment analysis.
#
# If split_locations is set to TRUE and locations are provided, the function performs GO term enrichment
# analysis separately for each location within each SV type. Plots are generated and saved for each location.


# Define the data sets to analyze
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the SV types to analyze
svtypes <- c('DEL', 'INS', 'INV','DUP:TANDEM', 'DUP:INT', 'BND')

# Specify whether to split locations for analysis
split_locations <-  TRUE

# Define the locations to analyze
location <- c("coding", "fiveUTR", "threeUTR", "intron", "intergenic", "spliceSite", "promoter")

# Specify the output directory
dir_out <- dir_results

# Perform GO term enrichment analysis for SV types and locations
go_per_sv_and_location(data_sets, allvar_df, svtypes, split_locations, location, dir_out)


# *************************************************
# HEATMAP
# This function generates a heatmap for each SVTYPE across different data sets.
# It filters the input data frame (allvar_df) based on the provided data sets (data_sets) and SV types (svtypes).
# The filtered data is then processed to create a matrix of gene counts for each data set and gene symbol.
# The resulting matrix is visualized as a heatmap using the Heatmap function from the ComplexHeatmap package.
# Each heatmap is saved as a PDF file in the specified output directory (dir_out) with the given file name (file_name).
#
# Arguments:
# Arguments:
#   - data_sets: A vector specifying the data sets to analyze.
#   - allvar_df: The dataframe containing SV data with gene annotations.
#   - svtypes: A character vector specifying the SV types to analyze.
#   - dir_out: The output directory where the results will be saved.
#   - file_name: The name of the output file.
#
# Output:
#   - Heatmap PDF files: A separate PDF file is generated for each SV type and saved in the specified output directory.
#   The PDF files contain the heatmaps visualizing the gene counts for each data set and gene symbol

# DCefine the data sets analyze
data_sets <- c('BL2', 'SUDHL5', 'DG75')

# Define the SV types 
svtypes <- c('DEL', 'INS', 'INV', 'DUP:TANDEM', 'DUP:INT', 'BND')

# Define the output directories
dir_out <- c()
for (sv in svtypes) {
  dir_out <- append(dir_out, paste(dir_results, "heatmap", sv, sep = "/"))
}

# Specify the file name for the heatmaps
file_name <- "heatmap.pdf"

# Call the heatmap function
heatmap(data_sets, allvar_df, svtypes, dir_out, file_name)















