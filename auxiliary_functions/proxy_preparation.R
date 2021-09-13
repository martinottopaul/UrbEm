proxy_cwd <- function(cams_domain,cams_origin,proxy)
{
  #### PROXY PREPROCESSING 
  # The 1kmx1km emission grid is defined so that each coarse grid cell of 7kmx7km 
  # is divided into 7x7 parts. From each of the proxy datasets a factor is then 
  # calculated indicating the proportion of each proxy data type in one high 
  # resolution grid cell within one coarse grid cell. These factors are used 
  # in order to downscale the respective emissions in the respective area.
  print(paste0("normalizing proxy: ", names(proxy)))
  #plot(proxy)
  
  ### create proxyulation proxy matrix
  ## gather unique concentration values to identify "bigger" cells from original CAMS raster
  conc_val <- unique(values(cams_domain))
  ## calculate amount of "smaller" cells in the resampled CAMS data, 
  ## that fit into the "bigger" cells of the original CAMS raster 
  ## based on resolutions of the different raster
  
  max_cell <- prod(res(cams_origin)/res(cams_domain))

  ## initialize lists
  ids <- list()
  proxy_norm_tiles <- list()
  
  ## for every unique concentration value (or "big" block from the original CAMS raster) ...
  for(c in 1:length(conc_val))
  {
    ## loop to collect cell ids of identical concentration values from concentration
    ## raster with the same extent & resolution (but with a higher original resolution)
    ## to identify blocks of same concentration values
    
    n <- 1
    i_id <- array()
    j_id <- array()
    
    for(i in 1:ncol(cams_domain))
    {
      for(j in 1:nrow(cams_domain))
      {
        if(cams_domain[i,j]==conc_val[c])
        {
          #print(paste(i,j, sep = ","))
          i_id <- c(i_id,i)
          j_id <- c(j_id,j)
          n=n+1
        } else
        {
          n=n+1
        }
      }
    }
    ## gather blocks of same concentration values in one raster
    ## by cropping it with the min and max cell ids from the original 
    ## concentration raster (same extent & resolution)
    ## and store it in a list of rasters
    ids[[c]] <- crop(cams_domain,extent(cams_domain,
                                        min(i_id, na.rm = T),
                                        max(i_id, na.rm = T),
                                        min(j_id, na.rm = T),
                                        max(j_id, na.rm = T)))
    #plot(proxy)
    #plot(cams_domain==conc_val[c])
    #plot(crop(proxy, ids[[c]]))
    
    ### use cropped raster with blocks of same concentration to
    ### crop the proxyulation density raster (same extent & resolution) and
    ### normalize it with the sum of the proxyulation density in each block AND
    ### additionally normalize it with the amount of "smaller cells" that fit into bigger "cells" to avoid
    ### too high values at the border (= incomplete "big cells"),
    ### (if there is no proxyulation value in a block divide it through number of cells for uniform distribution)
    ### and store it in a list of rasters
    
    #### avoid NA or zero cells
    if (sum(values(crop(proxy, ids[[c]])/sum(values(crop(proxy, ids[[c]])), na.rm = T)), na.rm= T) > 0)
    {
      proxy_norm_tiles[[c]] <- crop(proxy, ids[[c]])
      proxy_norm_tiles[[c]] <- proxy_norm_tiles[[c]]/sum(values(crop(proxy, ids[[c]])), na.rm = T)
    } else
    {
      proxy_norm_tiles[[c]] <- crop(proxy, ids[[c]])
      proxy_norm_tiles[[c]][] <- 1/ncell(crop(proxy, ids[[c]]))
    }
    #plot(proxy_norm_tiles[[c]])
    #print(sum(values(proxy_norm_tiles[[c]]), na.rm = T))
    
    ### correct grid cell values with less than 80% coverage of the original CAMS coarse resolution
    ### by applying the proportion of resolutions to all grid cells, that are affected
    ### otherwise we have values that are too high (especially at the boundaries)
    if (ncell(ids[[c]])/max_cell < 1)
    {
      proxy_norm_tiles[[c]] <- proxy_norm_tiles[[c]]*(ncell(ids[[c]])/max_cell)
    }
    #plot(proxy_norm_tiles[[c]])
    #print(sum(values(proxy_norm_tiles[[c]]), na.rm = T))
  }
  ### now the normalized proxyulation density blocks need to be merged to one raster
  proxy_norm <- proxy_norm_tiles[[1]]
  for(i in 2:length(proxy_norm_tiles))
  {
    proxy_norm <- merge(proxy_norm, proxy_norm_tiles[[i]])  
  }
  plot(proxy_norm)
  return(proxy_norm)
}