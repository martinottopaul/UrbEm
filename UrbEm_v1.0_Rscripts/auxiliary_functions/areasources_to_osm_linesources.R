area_to_osm_lines <- function(domain, area_emissions)
{  
  ### prepare spatial dataframe objects for processing of line sources
  # our domain
  domain_sf <- rasterToPolygons(domain)
  domain_sf <- st_as_sf(domain_sf)
  # our emissions as sf in right projection to read osm data
  emis <- rasterToPolygons(area_emissions)
  emis <- st_as_sf(emis)
  emis <- st_transform(emis, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
  
  ### get osm road network data with ALL street types based on the extent of our domain
  ### for railways use "railway" instead of "highway"
  print("downloading OSM data for specified domain")
  
  ### bounding box for osm data
  osm <- opq(bbox = st_bbox(emis))
  ### string with osm road types we are interested in!
  ### this could be modified with adding "residential", etc.
  street_types <- c("motorway", "motorway_link", "trunk", "trunk_link", 
                    "primary", "primary_link", "secondary", "secondary_link")
  ## create query for osm data download
  osm <- add_osm_feature(osm, "highway", street_types)
  ### download osm data
  osm <- osmdata_sf(osm)
  ## select only linesources
  osm <- osm$osm_lines[, c("highway")]
  
  ### show me what osm data is left
  #plot(osm, axes =T)
  
  ### transform back to UTM projection and crop it to the correct domain exten
  osm <- st_transform(osm, crs = crs(domain))
  
  ### distribution of all emissions to selected road types including 
  ### (a) road elements lengths
  ### (b) weighting by road type and for each cell separately
  options(warn=-1)
  remove(out)
  print("distributing area emissions: each grid cell value is attributed to road elements that fall into a grid cell")
  
  # initiate progress bar
  progress = txtProgressBar(min = 0, max = length(domain_sf$geometry), initial = 0, style = 3) 
  
  for(i in 1:length(domain_sf$geometry))
  {
    
    #algorithm based on VEIN package's "emis_dist" function by S. Ibara
    #modified by M. Ramacher
    
    ### gather all roads, which are in cell i from emissions raster
    osm_cell <- st_intersection(osm,domain_sf$geometry[i])
    
    ### what road types are in the gathered roads
    st <- unique(osm_cell$highway)
    
    ### create road-type-weighting for cell i based on road types in cell i
    if(length(st)>0)
    {
      ### how do we weight each road type?
      ## motorway & motorway link     = 5
      ## trunk & trunk linkl          = 3
      ## primary & primary link       = 2
      ## secondary & secondary link   = 1
      
      ### initalize weights with 0
      road_weights <- c(0,0,0,0)
      
      ### each road type gets a weighting if the road type exists in cell i
      if(any(grepl("motorway", st)))
      {
        road_weights[1] <- 10
      }
      if(any(grepl("trunk", st)))
      {
        road_weights[2] <- 5
      }
      if(any(grepl("primary", st)))
      {
        road_weights[3] <- 2
      }
      if(any(grepl("secondary", st)))
      {
        road_weights[4] <- 2
      }
      
      ### normalize weighting factors individually for cell i
      normalized_road_weight <- road_weights/sum(road_weights)
      
      ## length of streets (line elements) in cell i
      osm_cell$lkm1 <- as.numeric(sf::st_length(osm_cell))
      
      ### prepare table for distributing emissions of cell i:
      osm_cell <- osm_cell[, c("highway", "lkm1")]
      
      ### gather emissions for cell i from emissions raster
      NOx <- cams$NOx[i]
      NMVOC <- cams$NMVOC[i]
      CO <- cams$CO[i]
      SO2 <- cams$SO2[i]
      NH3 <- cams$NH3[i]
      PM25 <- cams$PM2.5[i]
      PM10 <- cams$PM10[i]
      
      ### select motorways and distribute emissions of cell i to motorways using normalized road weight AND street length
      ### lenght of motorway / sum of motorway length * weighting factor * total emissions
      osm_cell_m <- osm_cell[grep("motorway",osm_cell$highway),]
      osm_cell_m$NOx <- osm_cell_m$lkm1/sum(osm_cell_m$lkm1) * normalized_road_weight[1] * NOx
      osm_cell_m$NMVOC <- osm_cell_m$lkm1/sum(osm_cell_m$lkm1) * normalized_road_weight[1] * NMVOC
      osm_cell_m$CO <- osm_cell_m$lkm1/sum(osm_cell_m$lkm1) * normalized_road_weight[1] * CO
      osm_cell_m$SO2 <- osm_cell_m$lkm1/sum(osm_cell_m$lkm1) * normalized_road_weight[1] * SO2
      osm_cell_m$NH3 <- osm_cell_m$lkm1/sum(osm_cell_m$lkm1) * normalized_road_weight[1] * NH3
      osm_cell_m$PM25 <- osm_cell_m$lkm1/sum(osm_cell_m$lkm1) * normalized_road_weight[1] * PM25
      osm_cell_m$PM10 <- osm_cell_m$lkm1/sum(osm_cell_m$lkm1) * normalized_road_weight[1] * PM10
      
      ### select trunks (same as for motorways)
      osm_cell_t <- osm_cell[grep("trunk",osm_cell$highway), ]
      osm_cell_t$NOx <- osm_cell_t$lkm1/sum(osm_cell_t$lkm1) * normalized_road_weight[2] * NOx
      osm_cell_t$NMVOC <- osm_cell_t$lkm1/sum(osm_cell_t$lkm1) * normalized_road_weight[2] * NMVOC
      osm_cell_t$CO <- osm_cell_t$lkm1/sum(osm_cell_t$lkm1) * normalized_road_weight[2] * CO
      osm_cell_t$SO2 <- osm_cell_t$lkm1/sum(osm_cell_t$lkm1) * normalized_road_weight[2] * SO2
      osm_cell_t$NH3 <- osm_cell_t$lkm1/sum(osm_cell_t$lkm1) * normalized_road_weight[2] * NH3
      osm_cell_t$PM25 <- osm_cell_t$lkm1/sum(osm_cell_t$lkm1) * normalized_road_weight[2] * PM25
      osm_cell_t$PM10 <- osm_cell_t$lkm1/sum(osm_cell_t$lkm1) * normalized_road_weight[2] * PM10
      
      ### primarys (same as for motorways)
      osm_cell_p <- osm_cell[grep("primary",osm_cell$highway), ]
      osm_cell_p$NOx <- osm_cell_p$lkm1/sum(osm_cell_p$lkm1) * normalized_road_weight[3] * NOx
      osm_cell_p$NMVOC <- osm_cell_p$lkm1/sum(osm_cell_p$lkm1) * normalized_road_weight[3] * NMVOC
      osm_cell_p$CO <- osm_cell_p$lkm1/sum(osm_cell_p$lkm1) * normalized_road_weight[3] * CO
      osm_cell_p$SO2 <- osm_cell_p$lkm1/sum(osm_cell_p$lkm1) * normalized_road_weight[3] * SO2
      osm_cell_p$NH3 <- osm_cell_p$lkm1/sum(osm_cell_p$lkm1) * normalized_road_weight[3] * NH3
      osm_cell_p$PM25 <- osm_cell_p$lkm1/sum(osm_cell_p$lkm1) * normalized_road_weight[3] * PM25
      osm_cell_p$PM10 <- osm_cell_p$lkm1/sum(osm_cell_p$lkm1) * normalized_road_weight[3] * PM10
      
      ### secondaries (same as for motorways)
      osm_cell_s <- osm_cell[grep("secondary",osm_cell$highway), ]
      osm_cell_s$NOx <- osm_cell_s$lkm1/sum(osm_cell_s$lkm1) * normalized_road_weight[4] * NOx
      osm_cell_s$NMVOC <- osm_cell_s$lkm1/sum(osm_cell_s$lkm1) * normalized_road_weight[4] * NMVOC
      osm_cell_s$CO <- osm_cell_s$lkm1/sum(osm_cell_s$lkm1) * normalized_road_weight[4] * CO
      osm_cell_s$SO2 <- osm_cell_s$lkm1/sum(osm_cell_s$lkm1) * normalized_road_weight[4] * SO2
      osm_cell_s$NH3 <- osm_cell_s$lkm1/sum(osm_cell_s$lkm1) * normalized_road_weight[4] * NH3
      osm_cell_s$PM25 <- osm_cell_s$lkm1/sum(osm_cell_s$lkm1) * normalized_road_weight[4] * PM25
      osm_cell_s$PM10 <- osm_cell_s$lkm1/sum(osm_cell_s$lkm1) * normalized_road_weight[4] * PM10
      
      ### unite all weigthed road type emissions
      osm_cell_all <- rbind(osm_cell_m,osm_cell_t, osm_cell_p, osm_cell_s)
      head(osm_cell_all)
      
      ### dismiss length of road types column
      osm_cell_all <- osm_cell_all[, c("NOx", "NMVOC", "CO", "SO2", "NH3", "PM25", "PM10", "highway")]
      
      ### give correct names
      names(osm_cell_all) <- c("NOx", "NMVOC", "CO", "SO2", "NH3", "PM25", "PM10", "roadtype", "geometry")
      
      #plot(osm_cell_all["NOx"])
      #print(sum(osm_cell_all$NOx))
      
      if(!exists("out"))
      {
        out <- osm_cell_all
      } else {
        out <- rbind(out, osm_cell_all)
      }
    }
    setTxtProgressBar(progress,i)
  }
  options(warn=0)
  return(out)
}