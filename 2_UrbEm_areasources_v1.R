# UrbEm v1.0 areasource module (2/3)
# by M.O.P. Ramacher & A. Kakouri, 2021

#####################
### INPUT section ###
#####################

#libraries
library(raster)

# location of auxiliary functions
setwd("/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/Urbem_v1/auxiliary_functions/")
source("./proxy_preparation.R")
source("./proxy_distribution.R")
source("./areasources_e-prtr_pointsource_correction.R")

# Set folder locations
# where is the CAMS-REG-AP (nc) data located? (all files in one folder)?
emissions <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/CAMS-REG/CAMS-REG-AP_TNO_0.05x0.1_anthro/"

### choose an emissions downscaling option:
## 1 "top_down_proxy" distribution of total grid emission value using a normalized proxy
#proxy_method <- "top_down_proxy"
## 2 "coarse_cells_proxy" distribution applies a normalized proxy grid that was created based on the coarse grid-cells of the original CAMS resolution
proxy_method <- "coarse_cells_proxy"

# where are the proxy files (tif) located?
population_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif"
industry_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/LU_Industry.tif"
wastefacility_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/LU_Waste.tif"
agriculture_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/LU_Agriculture.tif"
shipping_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/LU_Shipping.tif"
aviation_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/LU_Airpots.tif"
offroad_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/Non_Road_Mob_Sources.tif"
ghsl_urbancentre <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/ghs_europe_iso3.tif"

# to which folder should the output be written to?
output <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_outputs/"

# set text string of output file
outputstring <- "Antweerp_UrbEm_CAMS-REG-AP_v31_2016_asrc"

# name of pointsource file created with Rscript 1_UrbEm_pointsources
psrc <- "Antweerp_e-prtr_2016_psrc.csv"

### define target raster definition: "our domain"
### Antweerp - Belgium
domain <- raster(nrow = 30, ncol = 30,
                 ymn = 5668000, xmn = 581000,
                 ymx = 5698000, xmx = 611000,
                 crs = "+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

### which year are we looking at
year <- 2016

####################
### END OF INPUT ###
####################

### fill the domain with values of "1"
domain <- setValues(domain,rep(1,ncell(domain)))

### read CAMS nc files
setwd(emissions)
files <- list.files(pattern = paste0(year,".nc"))

##which sectors
## choose from 
GNFR <- c("A_PublicPower", "B_Industry", "C_OtherStationaryComb", 
          "D_Fugitives", "E_Solvents", "F_RoadTransport", "G_Shipping", 
          "H_Aviation", "I_OffRoad", "J_Waste", "K_AgriLivestock", 
          "L_AgriOther", "SumAllSectors")

# create empty list for CAMS-REG emissions in our domain
GNFR_raster <- list()
# read and proxess CAMS *.nc files
for (i in 1:length(GNFR))
{
  ### read all pollutants for selected sector
  cams <- stack(files, varname = GNFR[i])
  
  ### rename each layer to is pollutant name (only)
  for (j in 1:nlayers(cams))
  {
    if(grepl(pattern = "ch4", names(cams[[j]]))) {names(cams[[j]]) <- "CH4"}
    else if (grepl(pattern = "co", names(cams[[j]]))) {names(cams[[j]]) <- "CO"}
    else if (grepl(pattern = "nh3", names(cams[[j]]))) {names(cams[[j]]) <- "NH3"}
    else if (grepl(pattern = "nmvoc", names(cams[[j]]))) {names(cams[[j]]) <- "NMVOC"}
    else if (grepl(pattern = "nox", names(cams[[j]]))) {names(cams[[j]]) <- "NOx"}
    else if (grepl(pattern = "pm10", names(cams[[j]]))) {names(cams[[j]]) <- "PM10"}
    else if (grepl(pattern = "pm2_5", names(cams[[j]]))) {names(cams[[j]]) <- "PM2.5"}
    else if (grepl(pattern = "so2", names(cams[[j]]))) {names(cams[[j]]) <- "SO2"}
  }
  
  ### project it to domain projection
  cams <- projectRaster(cams, crs = crs(domain))
  ### crop it to domain extent (with overlap at the boundaries)
  cams <- crop(cams, domain, snap = "out")
  ### multiply with 1E9 to get from Tg/year to kg/year
  cams <- cams*1E9
  
  ## save origin & resample it to target domain with same extent and resolution for proxy processing
  cams_origin <- cams[[1]]
  cams_domain <- resample(cams[[1]], domain, method = "ngb")
  
  ### resampling, cropping and mass-consistent disaggregation to our domains resolution and extent
  if (proxy_method == "top_down_proxy")
  {
    GNFR_raster[[i]] <- resample(cams, domain, method = "ngb")  
    GNFR_raster[[i]] <- GNFR_raster[[i]]/prod(res(cams)/res(domain))  
  } 
  
  if (proxy_method == "coarse_cells_proxy")
  {
    GNFR_raster[[i]] <- resample(cams, domain, method = "ngb")  
  }  
}

### assign names of GNFR sectors to  rasterstacks in list of rasters for (easier) further processing
names(GNFR_raster) <- GNFR
#plot(GNFR_raster[["A_PublicPower"]])
#plot(GNFR_raster[["B_Industry"]])
#plot(GNFR_raster[["C_OtherStationaryComb"]])
#plot(GNFR_raster[["D_Fugitives"]])
#plot(GNFR_raster[["E_Solvents"]])
#plot(GNFR_raster[["F_RoadTransport"]])
#plot(GNFR_raster[["G_Shipping"]])
#plot(GNFR_raster[["H_Aviation"]])
#plot(GNFR_raster[["I_OffRoad"]])
#plot(GNFR_raster[["J_Waste"]])
#plot(GNFR_raster[["K_AgriLivestock"]])
#plot(GNFR_raster[["L_AgriOther"]])

### PROXY Preparation
### prepare all proxy inputs as defined in the input section
## gather all proxies (paths) from environment variables, containing "_proxy"
proxies <- mget(ls()[grep("_proxy", ls())])
### create proxy raster with same extent, projection and resolution as domain raster
for (i in 1:length(proxies))
{
  ### read proxy tif in any resolution and projection
  proxy <- raster(proxies[[i]])
  
  ## two possibilities: 
  ## 1. project our domain to proxy projection, process and project back to our domain
  ## --> this leads to blurry proxies but makes sure no spatial is "lost" due to projecting data
  # project our domain to projection of the proxy
  domain_etrs <- projectRaster(domain,crs = crs(proxy))
  # crop proxy to extent of our domain
  proxy <- crop(proxy,domain_etrs,snap = "out")
  # bring it to our domains resolution with bilinear sampling
  proxy <- resample(proxy,domain_etrs,method = "bilinear")
  # project proxy to our domains projection and extent
  proxy <- projectRaster(proxy,domain,method = "ngb")
  #plot(proxy)
  
  ## 2. project our proxy directly to our domain's projection, process it
  ## --> this leads to blurry proxies but makes sure no spatial is "lost" due to projecting data
  #proxy <- projectRaster(proxy,domain, method ="bilinear")
  #plot(proxy)
  
  if (proxy_method == "top_down_proxy")
  {
    ### create uniformly normalized proxy grid for simple top-down distribution of grid total values
    assign(substring(names(proxies)[i],1,nchar(names(proxies)[i])-6), 
           proxy/sum(getValues(proxy),na.rm = T)) 
  } 
  
  if (proxy_method == "coarse_cells_proxy")
  {
    ### create normalized proxies for top-down distribution of values in each coarse grid cell:
    # The 1kmx1km emission grid is defined so that each coarse grid cell of 7kmx7km 
    # is divided into 7x7 parts. From each of the proxy datasets a factor is then 
    # calculated indicating the proportion of each proxy data type in one high 
    # resolution grid cell within one coarse grid cell. These factors are used 
    # in order to downscale the respective emissions in the respective area.
    # this concept is realized in the "proxy_cwd" function!
    assign(substring(names(proxies)[i],1,nchar(names(proxies)[i])-6), 
           proxy_cwd(cams_domain,cams_origin,proxy))
  } 
}

setwd(output)
### process GNFR sectors with sector-specific proxies
for (i in 1:(length(GNFR)-1))
{
  #### PUBLIC POWER ETC. WITH CLC Landuse category Industry
  if (names(GNFR_raster[i]) == "A_PublicPower" || 
      names(GNFR_raster[i]) == "B_Industry" ||
      names(GNFR_raster[i]) == "D_Fugitives")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rasterstack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = industry, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot(GNFR_raster[[i]])
  }
  if (names(GNFR_raster[i]) == "J_Waste")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rastersack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = wastefacility, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot( GNFR_raster[[i]])     
  }
  #### Residential Heating and Solvents with Population Density
  if (names(GNFR_raster[i]) == "C_OtherStationaryComb" || 
      names(GNFR_raster[i]) == "E_Solvents")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rastersack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = population, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot( GNFR_raster[[i]])      
  }
  #### Residential Heating and Solvents with Population Density
  if (names(GNFR_raster[i]) == "F_RoadTransport")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rastersack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = population, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot(GNFR_raster[[i]])      
    
    ### GHS map urban core
    #print(paste0("Total Nox before urban-centre weighting: ",
    #             sum(getValues(GNFR_raster[[i]]$NOx), na.rm = T) / 1000 / 1000," kt NOx"))
    
    ghs <- raster(ghsl_urbancentre)
    domain_wgs <- projectRaster(domain, crs = crs(ghs))
    ghs_domain <- crop(ghs, domain_wgs)
    ghs_domain <- projectRaster(ghs_domain, domain)
    ghs_domain[ghs_domain > 0] <- 3
    ghs_domain[ghs_domain == 0] <- 1
    #plot(ghs_domain)
    
    ### increase all pollutatnts  in urban core
    layernames <- names(GNFR_raster[[i]])
    GNFR_raster[[i]] <- ghs_domain * GNFR_raster[[i]]
    names(GNFR_raster[[i]]) <- layernames
    
    #plot(GNFR_raster[[i]]$NOx)
    #print(paste0("Total Nox after urban-centre weighting with a factor of ", 3, ": ",
    #sum(getValues(GNFR_raster[[i]]$NOx), na.rm = T) / 1000 / 1000," kt NOx"))
  }
  #### Shipping with Port Areas and Water Bodies at the coast
  if (names(GNFR_raster[i]) == "G_Shipping")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rastersack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = shipping, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot(GNFR_raster[[i]])     
  }
  #### Aviation with Airports
  if (names(GNFR_raster[i]) == "H_Aviation")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rastersack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = aviation, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot(GNFR_raster[[i]]) 
  }
  #### Non Road Mobility with Offroad Layer
  if (names(GNFR_raster[i]) == "I_OffRoad")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rastersack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = offroad, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot( GNFR_raster[[i]])
  }
  ### Agricultural sectors with Agriculture Landuse
  if (names(GNFR_raster[i]) == "K_AgriLivestock" || 
      names(GNFR_raster[i]) == "L_AgriOther")
  {
    ### plot sector emissions before applying proxy
    #plot(GNFR_raster[[i]])
    ### distribute sum of each pollutant rasterlayer in sector's rastersack with normalized proxy
    GNFR_raster[[i]] <- proxy_distribution(emissions = GNFR_raster[[i]], proxy = agriculture, proxy_method = proxy_method)
    ### plot sector emissions after applying proxy
    #plot( GNFR_raster[[i]])     
  }
}

### translate GNFR to SNAP (a sector specific function is necessary to split SNAP3 and SNAP4 based on E-PRTR?)
SNAP1 <- GNFR_raster[["A_PublicPower"]]
SNAP2 <- GNFR_raster[["C_OtherStationaryComb"]]
SNAP3 <- GNFR_raster[["B_Industry"]]*0.8
SNAP4 <- GNFR_raster[["B_Industry"]]*0.2
SNAP5 <- GNFR_raster[["D_Fugitives"]]
SNAP6 <- GNFR_raster[["E_Solvents"]]
SNAP7 <- GNFR_raster[["F_RoadTransport"]]
SNAP8 <- GNFR_raster[["G_Shipping"]]
SNAP9 <- GNFR_raster[["J_Waste"]]
SNAP10 <- GNFR_raster[["K_AgriLivestock"]]+GNFR_raster[["L_AgriOther"]]
SNAP11 <- GNFR_raster[["H_Aviation"]]+GNFR_raster[["I_OffRoad"]]

# plot(SNAP1$NOx)
# plot(SNAP2$NOx)
# plot(SNAP3$NOx)
# plot(SNAP4$NOx)
# plot(SNAP5$NOx)
# plot(SNAP6$NOx)
# plot(SNAP7$NOx)
# plot(SNAP8$NOx)
# plot(SNAP9$NOx)
# plot(SNAP10$NOx)
# plot(SNAP11$NOx)

# Process each SNAP sector specific based on information from E-PRTR
### Set area emissions to zero OR change them (difference of E-PRTR points and CAMS), 
### when the point source information is higher OR CAMS is higher
### this needs it's own function
setwd(output)
### read point source information and determine the split for SNAP 3 and 4 from pointsources based on BIMSCHVG
psrc <- read.csv(psrc, sep = ",", dec = ".", na.strings = "-999")
snap_raster <- list(SNAP1,SNAP2,SNAP3,SNAP4,SNAP5,SNAP6,SNAP7,SNAP8,SNAP9,SNAP10,SNAP11)

for(i in 1:length(snap_raster))
{
  #plot(snap_raster[[i]])
  snap_raster[[i]] <- area_point_correction(snap_raster[[i]], psrc, SNAP = i)
  #plot(snap_raster[[i]])
}

### create vector with snap sectors

## all sectors
snap_sectors <- c(1,2,3,4,5,6,7,8,9,10,11)

## without SNAP7
snap_sectors <- c(1,2,3,4,5,6,8,9,10,11)
snap_raster <- snap_raster[-7]

### create UECT area output table for EPISODE-CityChem aq simulations
### for each snap sector you can change the emission height
uect_area <- data.frame()

for (i in 1:length(snap_sectors))
{
  
  if(snap_sectors[[i]]== 1) {
    snaps <- 1  
    z <- 10
  } else if (snap_sectors[[i]]== 2) {
    snaps <- 2  
    z <- 10
  } else if (snap_sectors[[i]]== 3) {
    snaps <- 3  
    z <- 10
  } else if (snap_sectors[[i]]== 4) {
    snaps <- 4  
    z <- 10
  } else if (snap_sectors[[i]]== 5) {
    snaps <- 5  
    z <- 10
  } else if (snap_sectors[[i]]== 6) {
    snaps <- 6  
    z <- 0
  } else if (snap_sectors[[i]]== 7) {
    snaps <- 7  
    z <- 0
  } else if (snap_sectors[[i]]== 8) {
    snaps <- 8  
    z <- 10
  } else if (snap_sectors[[i]]== 9) {
    snaps <- 1
    z <- 10
  } else if (snap_sectors[[i]]== 10) {
    snaps <- 10  
    z <- 0
  } else if (snap_sectors[[i]]== 11) {
    snaps <- 10  
    z <- 10
  }
  
  ## coordinates
  sw <- coordinates(snap_raster[[i]])-res(domain)/2
  ne <- coordinates(snap_raster[[i]])+res(domain)/2
  
  ## pollutants
  nox <- values(snap_raster[[i]]$NOx)
  nmvoc <- values(snap_raster[[i]]$NMVOC)
  co <- values(snap_raster[[i]]$CO)
  so2 <- values(snap_raster[[i]]$SO2)
  nh3 <- values(snap_raster[[i]]$NH3)
  pm25 <- values(snap_raster[[i]]$PM2.5)
  pm10 <- values(snap_raster[[i]]$PM10)
  pollutants <- cbind(nox,nmvoc,co,so2,nh3,pm25,pm10)
  
  data <- cbind(snaps,sw,z,ne,z,pollutants)
  data[,8:14][is.na(data[,8:14])] <- 0
  colnames(data) <- c("snap","xcor_sw","ycor_sw","zcor_sw","xcor_ne","ycor_ne","zcor_ne","NOx","NMVOC","CO","SO2","NH3","PM2.5","PM10")
  uect_area <- rbind(uect_area,data)
}
head(uect_area)

### remove rows which contain only 0
uect_input_csv <- c()
for (i in 1:nrow(uect_area)) 
{
  if (sum(uect_area[i,c("NOx", "NMVOC","CO","SO2","NH3","PM2.5","PM10")], na.rm = F)>0)
    uect_input_csv <- rbind(uect_input_csv,uect_area[i,])
}

## write it to area csv
setwd(output)
write.csv(uect_input_csv, paste0(outputstring,".csv"), row.names = F, quote = F, na = "-999")

#rm(list = ls())
gc()

#####################
### END OF SCRIPT ###
#####################

### OPTIONAL section to check emission totals in domain:
### table with NOx sums of all snap sectors after processing(in kt/year)
#uect_input_csv[uect_input_csv == -999] <- 0
#sums <- vector()
#for (i in 1:length(unique(uect_input_csv$snap)))
#{
#  sums[i] <- sum(uect_input_csv$NOx[uect_input_csv$snap == unique(uect_input_csv$snap)[i]])/1000/1000 
#}
#sums_table <- (data.frame(unique(uect_input_csv$snap),sums))
#sums_table[NROW(sums_table)+1,] <- c("total", sum(sums_table$sums[1:NROW(sums_table)]))
#sums_table
#write.csv(sums_table, paste0(outputstring,"_SNAP_kt_year.csv"), quote = F, row.names = F)