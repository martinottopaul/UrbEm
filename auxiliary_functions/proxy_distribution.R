proxy_distribution <- function(emissions,proxy,proxy_method)
{
  temp <- list()
  for (j in 1:nlayers(emissions))
  {
    if (proxy_method == "top_down_proxy")
    {
      temp[[j]] <- sum(getValues(emissions[[j]]), na.rm = T)*proxy
    } 
    if (proxy_method == "coarse_cells_proxy")
    {
      temp[[j]] <- emissions[[j]]*proxy
    }
    names(temp[[j]]) <- names(emissions[[j]])
  }
  
  ### create rasterbrick out of different pollutant rasterlayers and replace original rasterstack in GNFR_raster
  emissions <- brick(temp)
  plot(emissions)
  
  return(emissions)
}  
