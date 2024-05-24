# Author: Shannon Laub

# Description: Contains functions for designScript.Rmd
# pickTopGDO()
# put description here...

# Function wishlist
# Check for overlaps + output dataframe w/ indices

# Guide preprocessing + add guide coords

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

