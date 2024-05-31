# Author: Shannon Laub

# Dependencies:
# GenomicRanges, libraryDesignFunctions.R, ...

# Description:
# Defines regions around each SNP. Combines SNPs that are close together into the same region.
# Exports a "regions" file with the region coordinates.
# Exports fasta files for CRISPick input.
# Exports a refseqLookup to convert Refseq chromosome names to chromosome numbers.

preprocessCRISPick <- function(){
# -----------------------------------------------------------------------------
  
  # Load Data
  
  wd <- getwd()
  wd # verify directory location
  
  # Load all SNPs.
  snp_list <- unique(read.xlsx(paste0(wd,"//input//snp_loc_hg19.xlsx")))
  rownames(snp_list) <- NULL
  snp_list$chr <- as.character(snp_list$chr)
  
  # Get refseq conversion for chromosome numbers. Needed for compatibility with CRISPick
  refseqLookup <- read.delim(paste0(wd,"//input//seq_report_hg19.tsv"))
  refseqLookup <- subset(refseqLookup, Role=="assembled-molecule")
  
  # Rename columns for convenience
  names(refseqLookup)[names(refseqLookup)=="Chromosome.name"] <- "chr"
  names(refseqLookup)[names(refseqLookup)=="RefSeq.seq.accession"] <- "refseq"
  refseqLookup <- refseqLookup %>%
    select(chr, refseq)
  
# -----------------------------------------------------------------------------
  # Preprocess Input for CRISPick
  
  # Combine nearby SNP regions.
  # Use 200bp for first pass. Expand region size for any 200bp regions 
  # that don't have enough good quality guides.
  regions200 <- findSnpRegions(200, snp_list)
  regions300 <- findSnpRegions(300, snp_list)
  regions400 <- findSnpRegions(400, snp_list)
  regions500 <- findSnpRegions(500, snp_list)
  
  # Convert to CRISPick format
  regions200 <- formatCrispick(regions200, refseqLookup)
  regions300 <- formatCrispick(regions300, refseqLookup)
  regions400 <- formatCrispick(regions400, refseqLookup)
  regions500 <- formatCrispick(regions500, refseqLookup)
  
  # Only need the 'crispick' column to run CRISPick for library design.
  crispick200Input <- regions200 %>% select(crispick)
  crispick300Input <- regions300 %>% select(crispick)
  crispick400Input <- regions400 %>% select(crispick)
  crispick500Input <- regions500 %>% select(crispick)
  
# -----------------------------------------------------------------------------
  # Export SNP regions as regions file and in fasta file for CRISPick Input.
  
  # Save SNP regions
  write.xlsx(regions200, file = paste0(wd, "//preprocessed//regions200.xlsx"))
  write.xlsx(regions300, file = paste0(wd, "//preprocessed//regions300.xlsx"))
  write.xlsx(regions400, file = paste0(wd, "//preprocessed//regions400.xlsx"))
  write.xlsx(regions500, file = paste0(wd, "//preprocessed//regions500.xlsx"))
  
  # Save CRISPick Inputs
  write.table(crispick200Input, paste0(wd, "//preprocessed//loc200.fasta"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(crispick300Input, paste0(wd, "//preprocessed//loc300.fasta"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(crispick400Input, paste0(wd, "//preprocessed//loc400.fasta"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  write.table(crispick500Input, paste0(wd, "//preprocessed//loc500.fasta"), 
              row.names=FALSE, col.names=FALSE, quote=FALSE)
  
  # Save refseqLookup for convenience.
  write.xlsx(refseqLookup, file = paste0(wd, "//input//refseqLookup.xlsx"))

}