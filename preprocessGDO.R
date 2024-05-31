# Author: Shannon Laub

# Dependencies...

# Description:


preprocessGDO <- function(){
# -----------------------------------------------------------------------------
  wd <- getwd()
  wd
  
  refseqLookup <- read.xlsx(paste0(wd, "//input//refseqLookup.xlsx"))
  
  regions200 <- read.xlsx(paste0(wd,"//preprocessed//regions200.xlsx"))
  regions300 <- read.xlsx(paste0(wd,"//preprocessed//regions300.xlsx"))
  regions400 <- read.xlsx(paste0(wd,"//preprocessed//regions400.xlsx"))
  
# Import CRISPick Output
  GDO200 <- read.xlsx(paste0(wd,"//crispick_output//sgrna_designs_200bp_i30.xlsx"))
  GDO300 <- read.xlsx(paste0(wd,"//crispick_output//sgrna_designs_300bp_i30.xlsx"))
  GDO400 <- read.xlsx(paste0(wd,"//crispick_output//sgrna_designs_400bp_i30.xlsx"))
  
  keepCol <- c("Input", "Reference.Sequence", "Orientation", "sgRNA.'Cut'.Position", 
               "sgRNA.Sequence", "sgRNA.Context.Sequence", "On-Target.Efficacy.Score", "Pick.Order")
  
  GDO200 <- select(GDO200, all_of(keepCol))
  GDO300 <- select(GDO300, all_of(keepCol))
  GDO400 <- select(GDO400, all_of(keepCol))
  
  # Alternative column selection:
  # guides <- guides %>%
  #   subset(select = -c(Quota, Target.Taxon, Target.Gene.ID, Target.Gene.Symbol, Target.Alias, 

# -----------------------------------------------------------------------------
  
# GC, TTTT, On-Target Efficacy >0.2 filter
  # Good GC content, no poly T, On-Target Efficacy > 0.2
  GDO200 <- filterGDO(GDO200)
  GDO300 <- filterGDO(GDO300)
  GDO400 <- filterGDO(GDO400)
  
  # Add GDO locations: chromosome number + start/end coordinates.
  GDO200 <- addGDOLoc(GDO200, refseqLookup)
  GDO300 <- addGDOLoc(GDO300, refseqLookup)
  GDO400 <- addGDOLoc(GDO400, refseqLookup)
  
  # Add the SNP Region
  GDO200 <- left_join(GDO200, select(regions200, SNP, crispick), join_by("Input"=="crispick"))
  GDO300 <- left_join(GDO300, select(regions300, SNP, crispick), join_by("Input"=="crispick"))
  GDO400 <- left_join(GDO400, select(regions400, SNP, crispick), join_by("Input"=="crispick"))
  
  # Add rownumber for reference later
  GDO200$n <- as.numeric(rownames(GDO200))
  GDO300$n <- as.numeric(rownames(GDO300))
  GDO400$n <- as.numeric(rownames(GDO400))
  
  # Debug
  # all(GDO200$start < GDO200$end) # CHECK that coordinate order is correct. TRUE if okay.
  # all(guides$n == as.numeric(rownames(guides))) #CHECK that n == rownumber. TRUE if okay.
  
# -----------------------------------------------------------------------------
  
# Export preprocessed GDO pools.
  write.xlsx(GDO200, file = paste0(wd, "//preprocessed//200bp_i30.xlsx"))
  write.xlsx(GDO300, file = paste0(wd, "//preprocessed//300bp_i30.xlsx"))
  write.xlsx(GDO400, file = paste0(wd, "//preprocessed//400bp_i30.xlsx"))
  
}