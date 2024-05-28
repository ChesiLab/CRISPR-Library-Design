# Author: Shannon Laub

# Description: Contains functions for designScript.Rmd
# findSnpRegions()

# pickTopGDO()
# put description here...

# Function wishlist
# Check for overlaps + output dataframe w/ indices
# Guide preprocessing + add guide coords

# source('LibraryDesignFunctions.R')

# -----------------------------------------------------------------------------
# Generate Regions

findSnpRegions <- function(regionLength, snpList){
  # set half length parameter
  halfLength = regionLength/2-1
  
  # Create SNP "regions". +/- 1/2 the total region width.
  snpRegion <- snpList %>%
    mutate(start = pos-halfLength,
           end = pos+halfLength)
  
  # Use GenomicRanges to find region overlaps. Create GRanges objects to use it.
  gRegion <- GRanges(seqnames = snpRegion$chr,
                     ranges = IRanges(start = snpRegion$start, end=snpRegion$end))
  
  gSnps <- GRanges(seqnames = snpRegion$chr,
                   ranges = IRanges(start=snpRegion$pos, end=snpRegion$pos),
                   SNP = snpRegion$snp)
  
  # Merge overlapping SNP regions
  gRegion <- reduce(gRegion)
  
  
  # Map SNPs back to merged regions.
  gOverlap <- findOverlaps(gRegion, gSnps)
  
  # Convert regions from GRanges object to dataframe
  regions <- as.data.frame(gRegion) 
  regions$SNP <- NA
  names(regions)[names(regions)=="seqnames"] <- "chr"
  
  # Convert snps from GRanges object to dataframe
  snps <- as.data.frame(gSnps)
  
  # Map SNPs to their regions.
  for(i in 1:nrow(regions)){
    # i = 5
    reg_idx <- which(gOverlap@from==i) #get index for region
    snp_idx <- gOverlap@to[reg_idx] # map to index for SNP
    
    regions[i, "SNP"] <- paste0(snps$SNP[snp_idx],collapse="|")
  }
  
  
  return(regions)
}

# -----------------------------------------------------------------------------
# Description: converts chromosome name to the relevant ENCODE refseq lookup.
# Add crispick column that formats the region location for CRISPick input.

formatCrispick <- function(regions, refseqLookup){
  lookup <- refseqLookup %>%
    select(chr, refseq)
  
  crispickInput <- select(regions, chr, start, end, SNP)
  
  crispickInput <- left_join(crispickInput, lookup, by="chr")
  
  # Add range for CRISPick
  crispickInput <- crispickInput %>%
    mutate(crispick = paste0(refseq,":+:",start,"-",end))
  
  return(crispickInput)
}

# -----------------------------------------------------------------------------
# filterGDO()
# Removes guides with bad (<0.25 or >0.75) GC content
# Removes guides that have TTTT+ sequence. (RNA poly termination seq)
# Removes guides that have poor predicted On-Target Efficiency (<0.2 Score)
# Filter parameters for GC and TTTT content are from FlashFry.

filterGDO <- function(GDO){
  # Add GC_percent and polyT column to evaluate each guide
  GDO <- GDO %>%
    mutate(GC_percent = (str_count(GDO$sgRNA.Sequence,"G")+str_count(GDO$sgRNA.Sequence,"C"))/20,
           polyT = str_count(GDO$sgRNA.Sequence,"TTTT")) # TTT shows up, so this seems to work.
  
  # Remove guides that don't meet filter criteria.
  GDO <- GDO %>%
    subset((GC_percent>=0.25) & (GC_percent<=0.75) & (polyT < 1) & (`On-Target.Efficacy.Score` > 0.2)) #filtered ~40 guides of 2000; for 200bp_i30 w/ efficacy filter, went from ~1.7k to 700 guides.
  
  # Renumber the rows after removal of bad guides.
  rownames(GDO) = NULL
  
  return(GDO)
}

# -----------------------------------------------------------------------------
# addGDOLoc()
# Adds the chromosome number and the start and end coordinates for each guide.

addGDOLoc <- function(GDO, refseqLookup){
  # Add GDO chromosome number.
  GDO <- left_join(GDO, refseqLookup, join_by("Reference.Sequence"=="refseq"))
  
  # Add GDO start and end coordinate.
  GDO <- GDO %>%
    mutate(
      start = case_when(
        Orientation == "sense" ~ `sgRNA.'Cut'.Position`-17,
        Orientation == "antisense" ~ `sgRNA.'Cut'.Position`-3),
      end = case_when(
        Orientation == "sense" ~ `sgRNA.'Cut'.Position`+2,
        Orientation == "antisense" ~ `sgRNA.'Cut'.Position`+16)
    )
  
  return(GDO)
}

# -----------------------------------------------------------------------------
# findGDOSelfOverlap()
# DEPENDS on GenomicRanges package.
# Finds which guides overlap. Outputs dataframe with indices of overlapping guides.

findGDOSelfOverlap <- function(GDO, chrCol, startCol, endCol, strandBool=TRUE, selfBool=FALSE){
  # Create GRanges object
  gGuides <- GRanges(seqnames = GDO[[chrCol]],
                      ranges = IRanges(start = GDO[[startCol]], end=GDO[[endCol]]))
  
  # Find guides that are overlapping.
  gOverlap <- findOverlaps(gGuides, ignore.strand=strandBool, drop.self = selfBool) 
  
  # Convert from GRanges object to dataframe.
  overlap <- as.data.frame(gOverlap)
  
  return(overlap)
}

# -----------------------------------------------------------------------------
# pickTopGDO()
# Description: picks the top 'n' guides from a pool of guides
pickTopGDO <- function(numTopGDO, guidePool, regions, overlapGDO) {
  topGDO <- guidePool[0,]
  
  for (i in 1:numTopGDO){
    
    for (j in 1:nrow(regions)) {
      # j = 1
      if (length(subset(guidePool, SNP==regions$SNP[j]))>0){
        
        top <- guidePool %>%
          subset(SNP==regions$SNP[j]) %>%
          arrange(Pick.Order) %>%
          dplyr::slice(1)
        
        # get the top pick (based on Pick Order), add this to top_guides.
        # then remove all guides that are overlapping with the top pick.
        topGDO <- bind_rows(topGDO, top)
        
        remove <- overlapGDO[overlapGDO$queryHits==top$n,2] # gets indices of guides that overlap the current top pick
        
        guidePool <- guidePool[!guidePool$n %in% remove,] # is this okay? are we getting the right index?
        
      }
    }
  }
  
  return(topGDO)
}

# -----------------------------------------------------------------------------
# getBadRegion()
# QC Checks: missing regions or not enough guides.
# Creates list of regions that need a bigger guide pool.

getBadRegion <- function(testPool, regions, numGDO){
  reg <- data.frame(table(testPool$SNP)) %>% # Create dataframe that counts number of guides/SNP
    setNames(c("SNP","Freq")) # Set column names.
  
  badRegion <- data.frame(SNP=character())
  
  badRegion <- bind_rows(badRegion, 
                        anti_join(regions, reg, by = "SNP") %>%
                          select(SNP)) # find missing regions
  
  badRegion <- bind_rows(badRegion, reg %>% 
                          subset(Freq < numGDO) %>%
                          select(SNP)) # find regions with not enough GDOs
  
  return(badRegion)
}

# -----------------------------------------------------------------------------
# replaceBadRegion()

replaceBadRegion <- function(badPool, biggerPool, allRegions, badRegions){
  # Removes any guides that are in the "bad regions" from the original pool
  combinedGDO <- subset(badPool, !(SNP %in% badRegions$SNP))
  
  # For bad regions, adds guides from the expanded guide pool.
  combinedGDO <- bind_rows(combinedGDO, 
                         subset(biggerPool, SNP %in% badRegions$SNP))
  
  return(combinedGDO)
}

