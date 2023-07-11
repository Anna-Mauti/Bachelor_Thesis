# Gene Annotation
gene_annotation <- function(data_sets, all_data, svtypes )
{

allvar_df <- data.frame ()

# load annotation from USCS database
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

for ( c in 1:length(data_sets))
{


for (sv in svtypes)  {

  # subset dataframe to only include one Sv Type and the current data
  data_sub <- all_data %>% 
  filter(data == data_sets[c]) %>%
  filter(SVTYPE %in% c(sv)) %>%
  dplyr::select(CHROM, POS, END, SVTYPE)


  if (sv == 'BND') {
   data_sub <- data_sub[complete.cases(data_sub[, c("POS")]), ]
   # here i changed end to be the pos(start) position because otherwise there would be only NAs
   sv_granges <- GRanges(seqnames = data_sub$CHROM,
                          ranges = IRanges(start = data_sub$POS, end = data_sub$POS),
                          metadata = data_sub$SVTYPE)
                      }
  

  else {
   # Subset the dataframe to only include rows without NA values in "POS" and "END"
    data_sub <- data_sub[complete.cases(data_sub[, c("POS", "END")]), ]
    # Create a GRanges object from the SV dataframe
    sv_granges <- GRanges(seqnames = data_sub$CHROM,
                          ranges = IRanges(start = data_sub$POS, end = data_sub$END),
                          metadata = data_sub$SVTYPE)
                         }


  # extract the Variants -> result: grange object
  allvar <- locateVariants(sv_granges, txdb, AllVariants())
        
  # match gene ids and gene names
  # Subset allvar to only include rows without NA values in "GENEID"
  keep <- !is.na(mcols(allvar)$GENEID)
  allvar <- allvar[keep,]

  # get Symbols for Gene Ids
  # Entrez Genes
  require('org.Hs.eg.db')
  gene_names <- getSYMBOL(mcols(allvar)$GENEID, 'org.Hs.eg.db')

  # Add the gene names as a new column to GRanges object: allvar
  mcols(allvar)$SYMBOL <- gene_names
  mcols(allvar)$SVTYPE <- sv
  mcols(allvar)$DATA <- data_sets[c]

  # bind all together to one big dataframe
  temp <- as.data.frame(allvar)
  allvar_df <- rbind(allvar_df, temp)
  
}
}


  return(allvar_df)
}
