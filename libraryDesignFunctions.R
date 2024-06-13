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
  # For DEBUG:
  # numTopGDO <- 5
  # guidePool <- read.xlsx(paste0(wd,"//preprocessed//200bp_i50.xlsx"))
  # regions <- read.xlsx(paste0(wd,"//preprocessed//regions200.xlsx"))
  # overlapGDO <- findSelfOverlap(guidePool, "chr", "start", "end")
  
  # saveRDS(topGDO, file = "topGDOj3.rds")
  # topGDO <- readRDS(file = "topGDOj3.rds")
  
  # Initialize df for topGDO.
  topGDO <- guidePool[0,]
  
  for (i in 1:nrow(regions)){ # 1:nrow(regions)
    # i = 1
    numOrientation <- c(sense=0, antisense=0)
    
    for (j in 1:numTopGDO) {
      # j = 2
      
      # Are there any guides left at all?
      if (length(subset(guidePool, SNP==regions$SNP[i]))>0){
        # -----------------------------------------------------------------------
        # If less than half of the guides have been picked, orientation doesn't matter:
        if(nrow(subset(topGDO, SNP==regions$SNP[i])) < (numTopGDO/2)){
          
          top <- guidePool %>%
            subset(SNP==regions$SNP[i]) %>%
            arrange(Pick.Order) %>%
            dplyr::slice(1)
          
          topGDO <- bind_rows(topGDO, top) # Get top pick, add to topGDO
          remove <- overlapGDO[overlapGDO$queryHits==top$n,2] # Get indices of guides overlapping top pick
          guidePool <- guidePool[!guidePool$n %in% remove,] # Remove overlapping guides from pool
          
          # If more than half of the guides have been picked, orientation could matter:
        } else {
          
          orientCt <- table(subset(topGDO, SNP==regions$SNP[i])$Orientation)
          numOrientation[names(orientCt)] <- orientCt
          
          # If more than half of guides are antisense:
          if (numOrientation["antisense"] >= (numTopGDO/2)){
            # subset for remaining guides that are antisense
            # then pick top guide
            if ( length(subset(guidePool, SNP==regions$SNP[i] & Orientation == "sense")) > 0){
              top <- guidePool %>%
                subset(SNP==regions$SNP[i] & Orientation=="sense") %>%
                arrange(Pick.Order) %>%
                dplyr::slice(1)
              
              topGDO <- bind_rows(topGDO, top) # Get top pick, add to topGDO
              remove <- overlapGDO[overlapGDO$queryHits==top$n,2] # Get indices of guides overlapping top pick
              guidePool <- guidePool[!guidePool$n %in% remove,] # Remove overlapping guides from pool
            }
            
            # If more than half of guides are sense:
          } else if (numOrientation["sense"] >= (numTopGDO/2)){
            # subset for remaining guides that are sense
            # then pick top guide
            if ( length(subset(guidePool, SNP==regions$SNP[i] & Orientation == "antisense")) > 0){
              top <- guidePool %>%
                subset(SNP==regions$SNP[i] & Orientation=="antisense") %>%
                arrange(Pick.Order) %>%
                dplyr::slice(1)
              
              topGDO <- bind_rows(topGDO, top) # Get top pick, add to topGDO
              remove <- overlapGDO[overlapGDO$queryHits==top$n,2] # Get indices of guides overlapping top pick
              guidePool <- guidePool[!guidePool$n %in% remove,] # Remove overlapping guides from pool
            }
            
            # If there aren't too many guides of either orientation:
          } else {
            # use full pool of remaining guides,
            # pick top guide like normal
            top <- guidePool %>%
              subset(SNP==regions$SNP[i]) %>%
              arrange(Pick.Order) %>%
              dplyr::slice(1)
            
            topGDO <- bind_rows(topGDO, top) # Get top pick, add to topGDO
            remove <- overlapGDO[overlapGDO$queryHits==top$n,2] # Get indices of guides overlapping top pick
            guidePool <- guidePool[!guidePool$n %in% remove,] # Remove overlapping guides from pool
          } #if/else for more than half of guides have been picked.
        } #if/else for how many guides have been picked already
        
        # -----------------------------------------------------------------------  
      } #end Check if any guides left for the region at all.
      
      
    } #end Loops till numTopGDO; for (j in 1:numTopGDO)
  }# end Loops through all regions; for (i in 1:nrow(regions))
  
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
# For regions that have been combined due to larger region width,
# 
replaceBadRegion <- function(badPool, biggerPool, allRegions, badRegions){
  if (nrow(badRegions>0)){
    # Removes any guides that are in the "bad regions" from the original pool
    combinedGDO <- subset(badPool, !(SNP %in% badRegions$SNP))
    
    # For bad regions, adds guides from the expanded guide pool.
    combinedGDO <- bind_rows(combinedGDO, 
                           subset(biggerPool, SNP %in% badRegions$SNP))
  } else {
    combinedGDO <- badPool
    print("No bad regions listed.")
  }
  
  if (!all(badRegions$SNP %in% allRegions$SNP)){
    changedReg <- subset(badRegions, !(SNP %in% allRegions$SNP))
    
    # print("Some regions have been combined. Please check: ")
    print(changedReg$SNP)
  }
  return(combinedGDO)
}

# -----------------------------------------------------------------------------
assembleLibrary <- function(numGDO, fileList, regionList, dirGDO, dirRegion,
                            onTargetFilter = -100, overlapLeniency = 0, pickByStrand = FALSE){
  initLib <- TRUE
  badRegion <- data.frame(SNP=character())
  
  for (regionName in regionList){
    # regionName <- regionList[2] # DEBUG
    
    if( (initLib==TRUE)|(nrow(badRegion)>0) ){
      
      # Get the region
      region <- read.xlsx(paste0(dirRegion, regionName, ".xlsx"))
      
      # Get coreesponding GDO
      width <- substring(regionName, nchar(regionName)-2, nchar(regionName)) # get the width
      file <- grep(width, fileList, value=TRUE) # get the corresponding GDO file
      GDO <- read.xlsx(paste0(dirGDO, file, ".xlsx"))
      
      if (onTargetFilter > -100){ # the values for this aren't straightforward.. # REFACTOR
        GDO <- filterGDO(GDO, GCbool = FALSE, Tbool = FALSE, minOnScore = onTargetFilter)
        if ("n" %in% names(GDO)){
          GDO <- subset(GDO, select = -c(n))
          GDO$n <- as.numeric(rownames(GDO))
        }
      }
      
      # Choose overlap leniency (or lack thereof)
      if (overlapLeniency == 0){
        overlapGDO <- findSelfOverlap(GDO, "chr","start","end")
      } else if (overlapLeniency > 0){
        GDO <- GDO %>%
          mutate(startL = start + overlapLeniency,
                 endL = end - overlapLeniency)
        overlapGDO <- findSelfOverlap(GDO, "chr","startL","endL")
      }
      
      # Choose whether picking is sensitive to even strand representation.
      if (pickByStrand == FALSE){
        topGDO <- pickTopGDO(numGDO, GDO, region, overlapGDO)
      } else if (pickByStrand == TRUE){
        topGDO <- pickTopGDObyStrand(numGDO, GDO, region, overlapGDO)
      }
      
      if (initLib == TRUE){
        initReg <- region
      } else if (initLib == FALSE){
        if (!all(badRegion$SNP %in% initReg$SNP)){ # Warning if new regions have been created.
          print(paste(regionName, "has merged bad regions together compared to initial regions."))
        }
        topGDO <- replaceBadRegion(prevTopGDO, topGDO, initReg, badRegion)
      }
      
      badRegion <- getBadRegion(topGDO, initReg, numGDO)
      
      if (nrow(badRegion)>0){
        print(paste(regionName, "has the following bad regions:"))
        print(badRegion$SNP)
      }
      
      prevTopGDO <- topGDO
      
      initLib <- FALSE
      
    }
  }
  return(topGDO)
}

# -----------------------------------------------------------------------------
# calcDistGDOtoSNP()
# Dependency: GenomicRanges, dplyr
# Calculates the distance from each GDO to its nearest SNP.
# Check that all GDOs are near something.
# For SNP list, expects snp, chr, and pos columns.
# For GDO list, expects chr, start, end, and SNP, sgRNA.Sequence columns.
calcDistGDOtoSNP <- function(SNP, GDO){
  gSNP <- GRanges(
    seqnames = SNP$chr,
    ranges = IRanges(start = SNP$pos, end = SNP$pos)
  )
  
  gGDO <- GRanges(
    seqnames = GDO$chr,
    ranges = IRanges(start = GDO$start, end = GDO$end)
  )
  
  nearest <- nearest(gGDO, gSNP)
  dist <- distance(gGDO, gSNP[nearest])
  
  nearestSNP <- SNP[nearest,]
  
  result <- cbind(
    select(GDO, SNP, sgRNA.Sequence, chr, start, end),
    nearestSNP,
    distance = dist
  )
  
  return(result)
}

# -----------------------------------------------------------------------------
# qc()
# Outputs some quick metrics on a given GDO library input.

qcSimple <- function(lib){
  print(paste("Mean Num GDO per Region:", mean(table(lib$SNP))))
  print(paste("Mean On-Target:", mean(lib$`On-Target.Efficacy.Score`)))
}

qcMeanOff <- function(lib){
  print(paste("Mean Off-Target, Tier I Bin I:", mean(lib$`#.Off-Target.Tier.I.Match.Bin.I.Matches`)))
  print(paste("Mean Off-Target, Tier II Bin I:", mean(lib$`#.Off-Target.Tier.II.Match.Bin.I.Matches`)))
  print(paste("Mean Off-Target, Tier III Bin I:", mean(lib$`#.Off-Target.Tier.III.Match.Bin.I.Matches`)))
  print(paste("Mean Off-Target, Tier I Bin II:", mean(lib$`#.Off-Target.Tier.I.Match.Bin.II.Matches`)))
  print(paste("Mean Off-Target, Tier II Bin II:", mean(lib$`#.Off-Target.Tier.II.Match.Bin.II.Matches`)))
  print(paste("Mean Off-Target, Tier III Bin II:", mean(lib$`#.Off-Target.Tier.III.Match.Bin.II.Matches`)))
}

qcSumOff <- function(lib){
  print(paste("Sum Off-Target, Tier I Bin I:", sum(lib$`#.Off-Target.Tier.I.Match.Bin.I.Matches`)))
  print(paste("Sum Off-Target, Tier II Bin I:", sum(lib$`#.Off-Target.Tier.II.Match.Bin.I.Matches`)))
  print(paste("Sum Off-Target, Tier III Bin I:", sum(lib$`#.Off-Target.Tier.III.Match.Bin.I.Matches`)))
  print(paste("Sum Off-Target, Tier I Bin II:", sum(lib$`#.Off-Target.Tier.I.Match.Bin.II.Matches`)))
  print(paste("Sum Off-Target, Tier II Bin II:", sum(lib$`#.Off-Target.Tier.II.Match.Bin.II.Matches`)))
  print(paste("Sum Off-Target, Tier III Bin II:", sum(lib$`#.Off-Target.Tier.III.Match.Bin.II.Matches`)))
}

qcSumRegOff <- function(lib){
  print(paste("Reg w OffTarget, Tier I Bin I:", nrow(subset(lib, `#.Off-Target.Tier.I.Match.Bin.I.Matches`>0)) ))
  print(paste("Reg w OffTarget, Tier II Bin I:", nrow(subset(lib, `#.Off-Target.Tier.II.Match.Bin.I.Matches`>0)) ))
  print(paste("Reg w OffTarget, Tier III Bin I:", nrow(subset(lib, `#.Off-Target.Tier.III.Match.Bin.I.Matches`>0)) ))
  
  print(paste("Reg w OffTarget, Tier I Bin II:", nrow(subset(lib, `#.Off-Target.Tier.I.Match.Bin.II.Matches`>0)) ))
  print(paste("Reg w OffTarget, Tier II Bin II:", nrow(subset(lib, `#.Off-Target.Tier.II.Match.Bin.II.Matches`>0)) ))
  print(paste("Reg w OffTarget, Tier III Bin II:", nrow(subset(lib, `#.Off-Target.Tier.III.Match.Bin.II.Matches`>0)) ))
}

qcPercRegOff <- function(lib){
  print(paste("%Reg OffTarget, Tier I Bin I:", nrow(subset(lib, `#.Off-Target.Tier.I.Match.Bin.I.Matches`>0))/nrow(lib) ))
  print(paste("%Reg OffTarget, Tier II Bin I:", nrow(subset(lib, `#.Off-Target.Tier.II.Match.Bin.I.Matches`>0))/nrow(lib) ))
  print(paste("%Reg OffTarget, Tier III Bin I:", nrow(subset(lib, `#.Off-Target.Tier.III.Match.Bin.I.Matches`>0))/nrow(lib) ))
  
  print(paste("%Reg OffTarget, Tier I Bin II:", nrow(subset(lib, `#.Off-Target.Tier.I.Match.Bin.II.Matches`>0))/nrow(lib) ))
  print(paste("%Reg OffTarget, Tier II Bin II:", nrow(subset(lib, `#.Off-Target.Tier.II.Match.Bin.II.Matches`>0))/nrow(lib) ))
  print(paste("%Reg OffTarget, Tier III Bin II:", nrow(subset(lib, `#.Off-Target.Tier.III.Match.Bin.II.Matches`>0))/nrow(lib) ))
}

qcStrand <- function(lib){
  orient <- lib %>%
    group_by(SNP, Orientation) %>%
    summarize(count = n(), .groups='drop') %>%
    pivot_wider(names_from = Orientation, values_from = count, values_fill = list(count = 0))
  
  return(orient)
}



