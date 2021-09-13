# UrbEm v1.0 linserource module (3/3)
# by M.O.P. Ramacher & A. Kakouri, 2021

#####################
### INPUT section ###
#####################

#libraries
library(raster)
library(sf)
library(osmdata)
library(lwgeom)

# auxiliary functions
setwd("/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/Urbem_v1/auxiliary_functions/")
source("./areasources_to_osm_linesources.R")
source("./proxy_preparation.R")

### Set folder locations
# where is the CAMS-REF-AP (nc) data located? (all files in one folder)?
emissions <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/CAMS-REG/CAMS-REG-AP_TNO_0.05x0.1_anthro/"

### do you want to apply population density as proxy before distributing to line sources?
### if no, emissions are just disaggregated to the resolution of our domain and distributed uniformly
pop_proxy <- "yes"
#pop_proxy <- "no"

### choose an emissions downscaling option:
##1 "top_down_proxy" distribution of total grid emission value using a normalized proxy
##2 "coarse_cells_proxy" distribution applies a normalized proxy grid that was created based on the coarse grid-cells of the original CAMS resolution
#proxy_method <- "top_down_proxy"
proxy_method <- "coarse_cells_proxy"

### do you want to use the GHSL urban-centre area to increase traffic emissions in the urban centre, following Kuik et al. 2019
### for "all" pollutants or for "nox" only?
centre <- "yes"
#centre <- "no"
centre_factor <- 3
pollutants <- "all"
#pollutants <- "nox"

# where are the proxy files (tif) located?
population_proxy <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif"
ghsl_urbancentre <- ("/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/proxies/ghs_europe_iso3.tif")

# to which folder should the output be written to?
output <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_outputs/"

# set text string of output file
outputstring <- "Antweerp_CAMS-REG-AP_v31_2016_lsrc"

### define target raster definition: "our domain"
### Antweerp - Belgium
domain <- raster(nrow = 30, ncol = 30,
                 ymn = 5668000, xmn = 581000,
                 ymx = 5698000, xmx = 611000,
                 crs = "+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# set year of CAMS-REG files
year <- 2016

####################
### END OF INPUT ###
####################

### fill the domain with values of "1"
domain <- setValues(domain,rep(1,ncell(domain)))

### read CAMS nc files
setwd(emissions)
files <- list.files(pattern = paste0(year,".nc"))

##which sector: F_RoadTransport for SNAP7
GNFR <- "F_RoadTransport"

### read all pollutants for selected sector
cams <- stack(files, varname = GNFR)

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
cams <- cams * 1E9
### resampling, cropping and mass-consistent disaggregation to our domains resolution and extent
cams_origin <- cams
cams_domain <- resample(cams, domain, method = "ngb")

if (pop_proxy == "no" || proxy_method == "top_down_proxy")
{
  cams <- resample(cams, domain, method = "ngb")
  cams <- cams / prod(res(cams_origin) / res(cams_domain))
  cams_sum <- sum(getValues(cams$NOx)) / 1000 / 1000
} else if (proxy_method == "coarse_cells_proxy")
{
  cams <- resample(cams, domain, method = "ngb")
  cams_sum <- sum(getValues(cams$NOx)) / 1000 / 1000
  cams_sum <- sum(getValues(cams$NOx)) / 1000 / 1000 / prod(res(cams_origin) / res(cams_domain))
}

### create a table for GNFR sector NOx sums in our domain [kt]
GNFR_sums_table <- data.frame(GNFR, cams_sum)
head(GNFR_sums_table)
#write.csv(GNFR_sums_table,paste0(output, "GNFR_TRAFFIC_kt_year.csv"),quote = F,row.names = F)

### process GNFR sectors with sector-specific proxies
### normalize the population density grid with the sum of population
if (pop_proxy == "yes")
{
  ### preparation of population proxy (reading from GHS population density)
  pop <- raster(population_proxy)
  ### project population density to the projection of our domain
  options(warn=-1)
  #pop <- projectRaster(pop, domain, method = "bilinear")
  ### crop population to extent of our domain
  #pop <- crop(pop, domain, snap = "out")
  
  ## 1. project our domain to proxy projection, process and project back to our domain
  ## --> this leads to blurry proxies but makes sure no spatial is "lost" due to projecting data
  # project our domain to projection of the proxy
  domain_etrs <- projectRaster(domain,crs = crs(pop))
  # crop proxy to extent of our domain
  pop <- crop(pop,domain_etrs,snap = "out")
  # bring it to our domains resolution with bilinear sampling
  pop <- resample(pop,domain_etrs,method = "bilinear")
  # project proxy to our domains projection and extent
  pop <- projectRaster(pop,domain,method = "ngb")
  
  plot(pop)
  
  if (proxy_method == "top_down_proxy")
  {
    ### create proxy (not normalized) for simple top-down distribution of grid total values
    pop_norm <- pop / sum(getValues(pop), na.rm = T)
    #plot(pop_norm * sum(getValues(cams$NOx)))
    sum(getValues(pop_norm * sum(getValues(cams$NOx))))
    
    temp <- list()
    for (j in 1:nlayers(cams))
    {
      temp[[j]] <- sum(getValues(cams[[j]]), na.rm = T) * pop_norm
      names(temp[[j]]) <- names(cams[[j]])
    }
  } else if (proxy_method == "coarse_cells_proxy")
  {
    ### create normalized proxies for top-down distribution of values in each coarse grid cell:
    # The 1kmx1km emission grid is defined so that each coarse grid cell of 7kmx7km
    # is divided into 7x7 parts. From each of the proxy datasets a factor is then
    # calculated indicating the proportion of each proxy data type in one high
    # resolution grid cell within one coarse grid cell. These factors are used
    # in order to downscale the respective emissions in the respective area.
    # this concept is realized in the "proxy_cwd" function!
    pop_norm <- proxy_cwd(cams_domain, cams_origin, pop)
    #plot(pop_norm * cams$NOx)
    sum(getValues(pop_norm * cams$NOx)) / 1000 / 1000
    
    temp <- list()
    for (j in 1:nlayers(cams))
    {
      temp[[j]] <- cams[[j]] * pop_norm
      names(temp[[j]]) <- names(cams[[j]])
    }
  }
  cams <- brick(temp)
}
#plot(cams)

### GHS map urban core
if (centre == "yes")
{
  print(paste0(
    "Total Nox before urban-centre weighting: ",
    sum(getValues(cams$NOx),na.rm = T) / 1000 / 1000,
    " kt NOx"
  ))
  ghs <- raster(ghsl_urbancentre)
  domain_wgs <- projectRaster(domain, crs = crs(ghs))
  ghs_domain <- crop(ghs, domain_wgs)
  ghs_domain <- projectRaster(ghs_domain, domain)
  ghs_domain[ghs_domain > 0] <- centre_factor
  ghs_domain[ghs_domain == 0] <- 1
  #plot(ghs_domain)
  
  if (pollutants == "all")
  {
    ### increase all pollutatnts  in urban core
    layernames <- names(cams)
    cams <- ghs_domain * cams
    names(cams) <- layernames
  } else if (pollutants == "nox")
  {
    ### only NOx increase in urban core
    cams$NOx <- cams$NOx * ghs_domain
  }
  #plot(cams$NOx)
  print(
    paste0(
      "Total Nox after urban-centre weighting with a factor of ",
      centre_factor,
      ": ",
      sum(getValues(cams$NOx),na.rm = T) / 1000 / 1000,
      " kt NOx"
    )
  )
}

### distribute area sources to line emissions with "area_to_osm_lines" function
### based on VEIN mechanism by Ibarra et al., modified by M.O.P.Ramacher
out <- area_to_osm_lines(domain = domain,area_emissions = cams)

### correct for "lost" emissions" that were in grid cells without roads
### due to resampling/downscaling procedure
out$NOx <- out$NOx*sum(getValues(cams$NOx), na.rm = T)/sum(out$NOx)
#plot(out["NOx"])
out$NMVOC <- out$NMVOC*sum(getValues(cams$NMVOC), na.rm = T)/sum(out$NMVOC)
out$CO <- out$CO*sum(getValues(cams$CO), na.rm = T)/sum(out$CO)
out$SO2 <- out$SO2*sum(getValues(cams$SO2), na.rm = T)/sum(out$SO2)
out$NH3 <- out$NH3*sum(getValues(cams$NH3), na.rm = T)/sum(out$NH3)
out$PM25 <- out$PM25*sum(getValues(cams$PM2.5), na.rm = T)/sum(out$PM25)
out$PM10 <- out$PM10*sum(getValues(cams$PM10), na.rm = T)/sum(out$PM10)
sum(out$NOx)/1000/1000
sum(getValues(cams$NOx), na.rm = T)/1000/1000
head(out)
plot(out["NOx"], axes = T, breaks = c(0,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000,10000))

### add road witdths based on road types to out table
unique(out$roadtype)
out$roadtype[out$roadtype == "motorway"] <- 20
out$roadtype[out$roadtype == "motorway_link"] <- 20
out$roadtype[out$roadtype == "trunk"] <- 16
out$roadtype[out$roadtype == "trunk_link"] <- 16
out$roadtype[out$roadtype == "primary"] <- 12
out$roadtype[out$roadtype == "primary_link"] <- 12
out$roadtype[out$roadtype == "secondary"] <- 12
out$roadtype[out$roadtype == "secondary_link"] <- 12
unique(out$roadtype)
names(out) <- c("NOx", "NMVOC", "CO", "SO2", "NH3", "PM25", "PM10", "width", "geometry")

### calcualte start and start and end points from line elements
start <- st_startpoint(out$geometry)
start <- unlist(st_geometry(start)) %>% 
  matrix(ncol=2,byrow=TRUE) %>% 
  data.frame() %>% 
  setNames(c("x_start","y_start"))

end <- st_endpoint(out$geometry)
end <- unlist(st_geometry(end)) %>% 
  matrix(ncol=2,byrow=TRUE) %>% 
  data.frame() %>% 
  setNames(c("x_end","y_end"))

### prepare table for output to UECT
snap <- 7
xcor_start <- as.integer(start$x_start)
ycor_start <- as.integer(start$y_start)
xcor_end <- as.integer(end$x_end)
ycor_end <- as.integer(end$y_end)
elevation <- 0
width <- out$width
### bring emissions from kg/year to g/s as required by UECT
NOx <- out$NOx*1E3/(365*24*3600)
NMVOC <- out$NMVOC*1E3/(365*24*3600)
CO <- out$CO*1E3/(365*24*3600)
SO2 <- out$SO2*1E3/(365*24*3600)
NH3 <- out$NH3*1E3/(365*24*3600)
PM25 <- out$PM25*1E3/(365*24*3600)
PM10 <- out$PM10*1E3/(365*24*3600)

### create uect format table
uect_lines <- cbind(snap,xcor_start,ycor_start,xcor_end,ycor_end,elevation,width,
                    NOx,NMVOC,CO,SO2,NH3,PM25,PM10)
colnames(uect_lines) <- c("snap","xcor_start","ycor_star","xcor_end","ycor_end","elevation",
                          "width","NOx","NMVOC","CO","SO2","NH3","PM2.5","PM10")
head(uect_lines)
setwd(output)
write.csv(uect_lines, paste0(outputstring,".csv"), row.names = F, quote = F, na = "-999")

#rm(list = ls())
gc()

#####################
### END OF SCRIPT ###
#####################

### OPTIONAL section to check emission totals in domain:
### table with NOx sums of all snap sectors after processing(in kt/year)
sums_table <- (data.frame(GNFR,sum(out$NOx)/1000/1000))
sums_table
#write.csv(sums_table, "UECT_input_TRAFFIC_kt_year_SNAP.csv", quote = F, row.names = F)
