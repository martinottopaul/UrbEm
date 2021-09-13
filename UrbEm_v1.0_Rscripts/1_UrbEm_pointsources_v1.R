# UrbEm v1.0 pointsource module (1/3)
# by M.O.P. Ramacher & A. Kakouri, 2021

#####################
### INPUT section ###
#####################

# load libraries
require(raster)
require(sf)
require(rgdal)

# set folder with e-prtr files
setwd("/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_inputs/E-PRTR/")

# set year of interest (must be included in e-prtr files)
year <- 2016

# define output string of csv file
outputstring <- "Antweerp_e-prtr_2016_psrc"

# to which folder should the output be written to?
output <- "/Users/martinramacher/Dropbox/martin/hzg/projects/UrbEm/_outputs/"

# pollutants of interest (must be included in e-prtr files)
pollutants <- c("CH4", "CO", "NH3", "NMVOC", "NOX", "PM10", "PM25", "SOX")

# which country?
country <- "Belgium"

# define domain extent for Antwerp - BE
domain <- raster(nrow = 1, ncol = 1,
                 ymn = 5668000, xmn = 581000,
                 ymx = 5698000, xmx = 611000,
                 crs = "+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

####################
### END OF INPUT ###
####################

# set domain value to 1
domain <- setValues(domain,rep(1,ncell(domain)))

# project UTM domain to WGS1984 lat/lon for cropping e-prtr data
domain_wgs <- projectRaster(domain, crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# read e-prtr facility reports with meta information
e_prtr_meta <- read.csv("dbo.PUBLISH_FACILITYREPORT.csv", sep = ",", dec = ".", stringsAsFactors = F)

# select only meta information for urban domain definition
e_prtr_meta_domain <- subset(e_prtr_meta, 
                             Long > xmin(domain_wgs) & 
                               Long < xmax(domain_wgs) &
                               Lat > ymin(domain_wgs) &
                               Lat < ymax(domain_wgs))

# read e-prtr pollutant report
e_prtr_pollutant_report <- read.csv("dbo.PUBLISH_POLLUTANTRELEASEANDTRANSFERREPORT.csv", sep = ",", 
                                    dec = ".", stringsAsFactors = F)
# get ID for pollutant reports by year and country
PRTR_ID <- e_prtr_pollutant_report$PollutantReleaseAndTransferReportID[e_prtr_pollutant_report$ReportingYear == year & 
                                                                         e_prtr_pollutant_report$CountryName == country]
# subset meta information by year and country ID
e_prtr_meta_domain_year <- subset(e_prtr_meta_domain, PollutantReleaseAndTransferReportID == PRTR_ID)

# read e-prtr pollutant release
e_prtr_pollutants <- read.csv("dbo.PUBLISH_POLLUTANTRELEASE.csv", sep = ",", dec = ".", stringsAsFactors = F)

# subset only air pollutants
e_prtr_pollutants <- subset(e_prtr_pollutants, ReleaseMediumCode == "AIR")

# create empty output vector and list
sums <- vector()
e_prtr_point_uect <- list()

# for each pollutant
for (i in 1:length(pollutants))
{
  # subset pollutant release information
  myData <- subset(e_prtr_pollutants, PollutantCode == pollutants[i])
  # merge with meta information by facility report ID
  myData <- merge(e_prtr_meta_domain_year, myData, by = "FacilityReportID")
  # convert to kg
  sums[i] <- sum(myData$TotalQuantity)/1000/1000
  
  # check for entries = zero
  if(sums[i]>0) {
    # get coordinates in latlon
    pts <- as.matrix(data.frame(myData$Long, myData$Lat))
    # project to domain projection (crs)
    pts <- project(pts, proj = as.character(crs(domain)))
    
    # fill list with dataframe for SNAP, coordinates, n.a. stack parameters, and emission total for each stack in domain
    e_prtr_point_uect[[i]] <- data.frame(
      as.numeric(myData$MainIASectorCode),
      pts[,1],
      pts[,2],
      rep(-999,length(pts[,1])),
      rep(-999,length(pts[,1])),
      rep(-999,length(pts[,1])),
      rep(-999,length(pts[,1])),
      myData$TotalQuantity)
 
    # if entry = 0
  } else {
    #list entry with 0
    e_prtr_point_uect[[i]] <- data.frame(0,0,0,0,0,0,0,0)
  }
  # set names
  names(e_prtr_point_uect[[i]]) <- c("snap","xcor","ycor","Hi","Vi","Ti","radi",pollutants[i])
}

# temporal array
temp <- e_prtr_point_uect[[1]]
# merge all tables in list to one table
for(i in 2:length(e_prtr_point_uect))
{
  temp <- merge(temp,e_prtr_point_uect[[i]], all = T)
}

## set sector classification E-PRTR to SNAP
temp$snap[temp$snap == 1] <- 1
temp$snap[temp$snap == 2] <- 3
temp$snap[temp$snap == 3] <- 3
temp$snap[temp$snap == 4] <- 4
temp$snap[temp$snap == 5] <- 9
temp$snap[temp$snap == 6] <- 3
temp$snap[temp$snap == 7] <- 10
temp$snap[temp$snap == 8] <- 3
temp$snap[temp$snap == 9] <- 6

# rearrange data
temp <- data.frame(temp[,1:7], temp$NOX, temp$NMVOC, temp$CO, temp$SOX, temp$NH3, temp$PM25, temp$PM10)
names(temp) <- c("snap","xcor","ycor","Hi","Vi","Ti","radi","NOx","NMVOC","CO","SO2","NH3","PM2.5","PM10")

## replace NA with -999
temp[is.na(temp)] <- -999

## remove empty rows
e_prtr_point_uect <-temp[!(temp$xcor==0),]

head(e_prtr_point_uect)

## write UECT_points.csv file
setwd(output)
write.csv(e_prtr_point_uect, paste0(outputstring,".csv"), row.names = F, quote = F)

#rm(list = ls())
gc()

#####################
### END OF SCRIPT ###
#####################

### OPTIONAL section to check emission totals in domain:
### table with NOx sums of all snap sectors after processing(in kt/year)
#e_prtr_point_uect[e_prtr_point_uect == -999] <- 0
#sums <- vector()
#for (i in 1:length(unique(e_prtr_point_uect$snap)))
#{
#  sums[i] <- sum(e_prtr_point_uect$NOx[e_prtr_point_uect$snap == unique(e_prtr_point_uect$snap)[i]])/1000/1000 
#}
#sums_table <- (data.frame(unique(e_prtr_point_uect$snap),sums))
#sums_table[NROW(sums_table)+1,] <- c("total", sum(sums_table$sums[1:NROW(sums_table)]))
#sums_table