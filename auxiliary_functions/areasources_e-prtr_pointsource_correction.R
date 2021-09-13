area_point_correction <- function(area,point,SNAP)
{
  #SNAP <- 2
  #area <- SNAP2
  #point <- psrc
  
  pollutants <- names(area)
  #skip CH4 (not in EPER, not needed for CityChem)
  pollutants <- pollutants[2:8]
  
  for(i in 1:length(pollutants))
  {
    #i <- 1
    if(sum(values(area[[pollutants[i]]]), na.rm = T)-sum(subset(point, snap == SNAP)[,pollutants[i]], na.rm = T) <= 0) 
    {
      area[[pollutants[i]]] <- 0
      
    } else
    {
      area[[pollutants[i]]] <- area[[pollutants[i]]]/sum(values(area[[pollutants[i]]]), na.rm = T)*
        (sum(values(area[[pollutants[i]]]), na.rm = T)-sum(subset(point, snap == SNAP)[,pollutants[i]], na.rm = T))
    }
  }
  return(area)
  plot(area)
}

