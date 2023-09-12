# Author: Anna Gelbe
# Date: 13.07.2023

install_packages <- function(r_lib)

{

install.packages("vcfR", lib=r_lib)
install.packages("memuse", lib="/projtect/meissner_training/Anna/Bachelor/R_libs")
## circos plot 
install.packages("circlize", lib=r_lib)

## Gene Annotation
install.packages("GenomicRanges", lib=r_lib)

#
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("VariantAnnotation", lib=r_lib)

#
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene", lib=r_lib, force = TRUE)

#
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("annotate", lib=r_lib, force = TRUE)

#
install.packages("NLP", lib=r_lib)

# GO analysis
 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("clusterProfiler", lib=r_lib)

#
 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("AnnotationDbi", lib=r_lib)

#
 if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", lib=r_lib, force = TRUE)

#
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("StructuralVariantAnnotation", lib=r_lib)

# for heatmap 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap", lib=r_lib, force = TRUE)

}