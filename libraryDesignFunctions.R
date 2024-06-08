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
  
  # crispickInput <- select(regions, chr, start, end, SNP)
  
  crispickInput <- left_join(regions, lookup, by="chr")
  
  # Add range for CRISPick
  crispickInput <- crispickInput %>%
    mutate(crispick = paste0(refseq,":+:",start,"-",end))
  
  return(crispickInput)
}

# -----------------------------------------------------------------------------
# filterGDO()

# Default options:
# Removes guides with bad (<0.25 or >0.75) GC content
# Removes guides that have TTTT+ sequence. (RNA poly termination seq)
# Does not remove guides with poor On-Target Efficacy. (Suggested: remove <0.2 score)
# Filter parameters for GC and TTTT content are from FlashFry.

filterGDO <- function(GDO, GCbool = TRUE, Tbool = TRUE, minOnScore = -1000){
  # Add GC_percent and polyT column to evaluate each guide
  GDO <- GDO %>%
    mutate(GC_percent = (str_count(GDO$sgRNA.Sequence,"G")+str_count(GDO$sgRNA.Sequence,"C"))/20,
           polyT = str_count(GDO$sgRNA.Sequence,"TTTT")) # TTT shows up, so this seems to work.
  
  # Remove guides that don't meet filter criteria.
  if (GCbool == TRUE){
    GDO <- GDO %>%
      subset((GC_percent>=0.25) & (GC_percent<=0.75))
  }
  
  if (Tbool == TRUE){
    GDO <- GDO %>%
      subset(polyT < 1)
  }
  
  # Note that minOnScore can be < 0
  GDO <- GDO %>%
    subset((`On-Target.Efficacy.Score` > minOnScore)) #filtered ~40 guides of 2000; for 200bp_i30 w/ efficacy filter, went from ~1.7k to 700 guides.
  
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
# findSelfOverlap()
# DEPENDS on GenomicRanges package.
# Finds which coords overlap. Outputs dataframe with indices of overlapping coords.

findSelfOverlap <- function(regions, chrCol, startCol, endCol, strandBool=TRUE, selfBool=FALSE){
  # Create GRanges object
  gGuides <- GRanges(seqnames = regions[[chrCol]],
                      ranges = IRanges(start = regions[[startCol]], end=regions[[endCol]]))
  
  # Find guides that are overlapping.
  gOverlap <- findOverlaps(gGuides, ignore.strand=strandBool, drop.self = selfBool) 
  
  # Convert from GRanges object to dataframe.
  overlap <- as.data.frame(gOverlap)
  
  return(overlap)
}

# -----------------------------------------------------------------------------
# findOverlap()
# DEPENDS on GenomicRanges package.
# Similar to findSelfOverlap(), but compares two different lists of coords.

findOverlap <- function(reg1, reg2, 
                           chrCol1, startCol1, endCol1, 
                           chrCol2, startCol2, endCol2, 
                           strandBool=TRUE){
  
  gReg1 <- GRanges(seqnames = reg1[[chrCol1]],
                     ranges = IRanges(start = reg1[[startCol1]], end=reg1[[endCol1]]))
  
  gReg2 <- GRanges(seqnames = reg2[[chrCol2]],
                   ranges = IRanges(start = reg2[[startCol2]], end=reg2[[endCol2]]))
  
  gOverlap <- findOverlaps(gReg1, gReg2, ignore.strand=strandBool)
  
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
# pickTopGDObyStrand()
# todo: refactor/functionalize to simplify logic.

# numTopGDO = 5
# guidePool = GDO200
# regions = regions200
# overlapGDO = GDO200Overlap

pickTopGDObyStrand <- function(numTopGDO, guidePool, regions, overlapGDO){
  topGDO <- guidePool[0,]
  
  for (i in 1:nrow(regions)){ #nrow(regions)
    # i = 1
    for (j in 1:numTopGDO) {
      # j = 5
      # Check if there is more than one guide left in the pool for the given region.
      if (length(subset(guidePool, SNP==regions$SNP[i]))>0){
        
        # If less than half of the guides have been picked, orientation doesn't matter:
        if(nrow(subset(topGDO, SNP==regions$SNP[i])) < (numTopGDO/2)){
          
          top <- guidePool %>%
            subset(SNP==regions$SNP[i]) %>%
            arrange(Pick.Order) %>%
            dplyr::slice(1)
          
          # get the top pick (based on Pick Order), add this to top_guides.
          # then remove all guides that are overlapping with the top pick.
          topGDO <- bind_rows(topGDO, top)
          
          # gets indices of guides that overlap the current top pick
          remove <- overlapGDO[overlapGDO$queryHits==top$n,2] 
          
          # Removes guides that overlap the top hit from the remaining guide pool.
          guidePool <- guidePool[!guidePool$n %in% remove,] 
          
          # If more than half of the guides have been picked, orientation could matter:
        } else {
          numOrientation <- table(topGDO$Orientation)
          
          # If more than half of guides are antisense:
          if (numOrientation["antisense"] > (numTopGDO/2)){
            # subset for remaining guides that are antisense
            # then pick top guide
            if (length(subset(guidePool, SNP==regions$SNP[i] & Orientation == "sense")>0)){
              top <- guidePool %>%
                subset(SNP==regions$SNP[i] & Orientation=="sense") %>%
                arrange(Pick.Order) %>%
                dplyr::slice(1)
              
              topGDO <- bind_rows(topGDO, top)
              remove <- overlapGDO[overlapGDO$queryHits==top$n,2] 
              guidePool <- guidePool[!guidePool$n %in% remove,]
            }
            
            # If more than half of guides are sense:
          } else if (numOrientation["sense"] > (numTopGDO/2)){
            # subset for remaining guides that are sense
            # then pick top guide
            if (length(subset(guidePool, SNP==regions$SNP[i] & Orientation == "antisense")>0)){
              top <- guidePool %>%
                subset(SNP==regions$SNP[i] & Orientation=="antisense") %>%
                arrange(Pick.Order) %>%
                dplyr::slice(1)
              
              topGDO <- bind_rows(topGDO, top)
              remove <- overlapGDO[overlapGDO$queryHits==top$n,2] 
              guidePool <- guidePool[!guidePool$n %in% remove,]
            }
            
            # If there aren't too many guides of either orientation:
          } else {
            # use full pool of remaining guides,
            # pick top guide like normal
            top <- guidePool %>%
              subset(SNP==regions$SNP[i]) %>%
              arrange(Pick.Order) %>%
              dplyr::slice(1)
            
            topGDO <- bind_rows(topGDO, top)
            remove <- overlapGDO[overlapGDO$queryHits==top$n,2]
            guidePool <- guidePool[!guidePool$n %in% remove,]
          } 
          
        } #end if/else for guidepool > numTopGDO/2
      } #end if guidePool$region[i] > 0
    } #end for (j in 1:numTopGDO)
  }# end for (i in 1:nrow(regions))
  
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

