
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

