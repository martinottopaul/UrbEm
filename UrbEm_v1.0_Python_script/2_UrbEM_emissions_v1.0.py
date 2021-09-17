import pandas as pd
import numpy as np
import gdal, osr
import geopandas as gpd
import os, sys, fnmatch, inspect
import rasterio as rio
import rasterio.mask
import earthpy.spatial as es
from shapely.geometry import Polygon, Point
from geopandas import GeoDataFrame
from earthpy import clip as cl
import fiona
from osgeo import ogr
import gc
import geotable
from pyproj import Proj, transform
from geotable.projections import LONGITUDE_LATITUDE_PROJ4
from pandas import DataFrame
from shapely import wkt
import time
from os import listdir
from os.path import isfile, join
import xarray as xr
from rasterio import features


start_time = time.time()


### Set working dir & create folders

theWD = os.path.dirname(inspect.getfile(inspect.currentframe())) + "\\"
theWD = str(theWD.replace("\\", "/"))
print(theWD)

##################################################
Year = 2016
Region = 'Athens'
uc_increase_factor = 1
InFolder = theWD + "Athens/Results_Folder"
Proxy_Folder = theWD + "Athens/Results_Folder/Proxies"
OutFolder = InFolder + "/Emissions/" + str(Year) + "/V3/Increase_Factor_" + str(uc_increase_factor) + "/Results/"
Domain_with_sea = 'YES'
xmin = 716397
xmax = 761397
ymin = 4191261
ymax = 4236261
cell_size = 1000
epsg_code = 32634
crs_utm = 'EPSG:32634'
crs_wgs = 'EPSG:4326'
crs_corine = 'EPSG:3035'
country_code = "GRC"
source_type = "A"
country = "Greece"
Eurostat_countries_shp = theWD + "Input_Data/ref-countries-2020-01m.shp/CNTR_RG_01M_2020_4326.shp/CNTR_RG_01M_2020_4326.shp"
##################################################







if not os.path.exists(OutFolder):
    os.makedirs(OutFolder)

##pngFolder = OutFolder + "/png/"
##if not os.path.exists(pngFolder):
##    os.makedirs(pngFolder)

csvFolder = OutFolder + "/Results_CSVs/"
if not os.path.exists(csvFolder):
    os.makedirs(csvFolder)



### Change scientific format from default
pd.set_option('display.float_format', lambda x: '%.16f' % x)





#**************************************************************************************************************
#*************************************** Create small grid for domain ***************************************#
#**************************************************************************************************************



### Set cell size and resolution - Coordinates and Grid Parameters
# For proxies


lat_point_list = [ymin, ymin, ymax, ymax]
lon_point_list = [xmin, xmax, xmax, xmin]

rows = int(np.ceil((ymax-ymin) /  cell_size))
cols = int(np.ceil((xmax-xmin) / cell_size))

proxy_width = rows
proxy_height = cols

XleftOrigin = xmin
XrightOrigin = xmin + cell_size
YtopOrigin = ymax
YbottomOrigin = ymax- cell_size
polygons = []
for i in range(cols):
	Ytop = YtopOrigin
	Ybottom =YbottomOrigin
	for j in range(rows):
		polygons.append(Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])) 
		Ytop = Ytop - cell_size
		Ybottom = Ybottom - cell_size
	XleftOrigin = XleftOrigin + cell_size
	XrightOrigin = XrightOrigin + cell_size

### Set the coordinate system
# EPSG Code for WGS84: 4326
# EPSG Code for WGS84_UTM_Zone_34N: 32634
# EPSG Code for ETRS_1989_LAEA: 3035 ~ *********** Corine Coordinate System *********** ~

# Create grid from coordinates and write shapefile
grid = gpd.GeoDataFrame({'geometry':polygons}, crs=crs_utm)
grid['grid_index'] = grid.index
grid.to_file(OutFolder + "grid.shp", driver="ESRI Shapefile")
grid_corine = grid.to_crs(crs_corine)
grid_wgs = grid.to_crs(crs_wgs)


# Create point form polygon - grid centroid coordinates
grid_points = grid.copy()

grid_points['geometry'] = grid_points['geometry'].centroid
grid_points['xcor'] = grid_points['geometry'].x
grid_points['ycor'] = grid_points['geometry'].y
print(grid_points)
grid_points.index = range(len(grid_points))
grid_points_coords = [(x,y) for x, y in zip(grid_points.ycor, grid_points.xcor)]



grid_wgs_bounds = grid_wgs.bounds
xmin_wgs = grid_wgs_bounds['minx'].min()
xmax_wgs = grid_wgs_bounds['maxx'].max()
ymin_wgs = grid_wgs_bounds['miny'].min()
ymax_wgs = grid_wgs_bounds['maxy'].max()

print(xmin_wgs)
print(xmax_wgs)
print(ymin_wgs)
print(ymax_wgs)
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








#**************************************************************************************************************
#************************************* Create domain for sea/land clip ***************************************#
#**************************************************************************************************************



xmin_sea=(xmin - 1000)
xmax_sea=(xmax + 1000)
ymin_sea=(ymin - 1000)
ymax_sea=(ymax + 1000)

lat_point_list = [ymin_sea, ymin_sea, ymax_sea, ymax_sea]
lon_point_list = [xmin_sea, xmax_sea, xmax_sea, xmin_sea]

rows = int(np.ceil((ymax_sea-ymin_sea) /  cell_size))
cols = int(np.ceil((xmax_sea-xmin_sea) / cell_size))

proxy_width = rows
proxy_height = cols

XleftOrigin = xmin_sea
XrightOrigin = xmin_sea + cell_size
YtopOrigin = ymax_sea
YbottomOrigin = ymax_sea - cell_size
polygons_sea = []
for i in range(cols):
	Ytop = YtopOrigin
	Ybottom =YbottomOrigin
	for j in range(rows):
		polygons_sea.append(Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])) 
		Ytop = Ytop - cell_size
		Ybottom = Ybottom - cell_size
	XleftOrigin = XleftOrigin + cell_size
	XrightOrigin = XrightOrigin + cell_size

### Set the coordinate system
# EPSG Code for WGS84: 4326
# EPSG Code for WGS84_UTM_Zone_34N: 32634
# EPSG Code for ETRS_1989_LAEA: 3035 ~ *********** Corine Coordinate System *********** ~

# Create grid from coordinates and write shapefile
grid_sea = gpd.GeoDataFrame({'geometry':polygons_sea}, crs=crs_utm)
grid_sea['grid_index'] = grid_sea.index
grid_sea.to_file(OutFolder + "grid_sea.shp", driver="ESRI Shapefile")
grid_sea_wgs = grid_sea.to_crs(crs_wgs)

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








#**************************************************************************************************************
#*************************************** Create CAMS grid for domain ***************************************#
#**************************************************************************************************************



cams_width = 0.1
cams_height = 0.05

xmin_cams_wgs=xmin_wgs
xmax_cams_wgs=xmax_wgs
ymin_cams_wgs=ymin_wgs
ymax_cams_wgs=ymax_wgs

lat_point_list = [ymin_cams_wgs, ymin_cams_wgs, ymax_cams_wgs, ymax_cams_wgs]
lon_point_list = [xmin_cams_wgs, xmax_cams_wgs, xmax_cams_wgs, xmin_cams_wgs]

rows = int(np.ceil((ymax_cams_wgs-ymin_cams_wgs) /  cams_height))
cols = int(np.ceil((xmax_cams_wgs-xmin_cams_wgs) / cams_width))


XleftOrigin = xmin_cams_wgs
XrightOrigin = xmin_cams_wgs + cams_width
YtopOrigin = ymax_cams_wgs
YbottomOrigin = ymax_cams_wgs- cams_height
cams_polygons = []
for i in range(cols):
	Ytop = YtopOrigin
	Ybottom =YbottomOrigin
	for j in range(rows):
		cams_polygons.append(Polygon([(XleftOrigin, Ytop), (XrightOrigin, Ytop), (XrightOrigin, Ybottom), (XleftOrigin, Ybottom)])) 
		Ytop = Ytop - cams_height
		Ybottom = Ybottom - cams_height
	XleftOrigin = XleftOrigin + cams_width
	XrightOrigin = XrightOrigin + cams_width

### Set the coordinate system
# EPSG Code for WGS84: 4326
# EPSG Code for WGS84_UTM_Zone_34N: 32634
# EPSG Code for ETRS_1989_LAEA: 3035 ~ *********** Corine Coordinate System *********** ~

### ~ CAMS ~ ###
### ~ Create grid from coordinates and write shapefile ~ ###
cams_grid = gpd.GeoDataFrame({'geometry':cams_polygons}, crs=crs_wgs)
cams_grid_utm = cams_grid.to_crs(crs_utm)
cams_grid.to_file(OutFolder + "cams_grid.shp", driver="ESRI Shapefile")

cams_grid_utm_cliped_temp = cl.clip_shp(cams_grid_utm, grid)
cams_grid_utm_cliped_temp.to_file(OutFolder + "cams_grid_utm_cliped.shp", driver="ESRI Shapefile")
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
















#**************************************************************************************************************
#**************************************************************************************************************
#*************************************** Calculate Point Source Emissions ***************************************#
#**************************************************************************************************************
#**************************************************************************************************************
















e_prtr_meta = pd.read_csv(InFolder + "/E-PRTR/dbo.PUBLISH_FACILITYREPORT.csv", sep=',', encoding = "ISO-8859-1")

e_prtr_meta_croped = e_prtr_meta.loc[(e_prtr_meta['Long'] > xmin_wgs) & (e_prtr_meta['Long'] < xmax_wgs) & (e_prtr_meta['Lat'] > ymin_wgs) & (e_prtr_meta['Lat'] < ymax_wgs)]

e_prtr_pollutant_report = pd.read_csv(InFolder + "/E-PRTR/dbo.PUBLISH_POLLUTANTRELEASEANDTRANSFERREPORT.csv", sep=',', encoding = "ISO-8859-1")




if Year == 2016:
    PRTR_ID = e_prtr_pollutant_report['PollutantReleaseAndTransferReportID'].loc[(e_prtr_pollutant_report['ReportingYear'] == Year) & (e_prtr_pollutant_report['CountryName'] == country)]
elif Year == 2017:
    PRTR_ID = e_prtr_pollutant_report['PollutantReleaseAndTransferReportID'].loc[(e_prtr_pollutant_report['ReportingYear'] == (Year - 1)) & (e_prtr_pollutant_report['CountryName'] == country)]
elif Year == 2018:
    PRTR_ID = e_prtr_pollutant_report['PollutantReleaseAndTransferReportID'].loc[(e_prtr_pollutant_report['ReportingYear'] == (Year - 2)) & (e_prtr_pollutant_report['CountryName'] == country)]
elif Year == 2019:
    PRTR_ID = e_prtr_pollutant_report['PollutantReleaseAndTransferReportID'].loc[(e_prtr_pollutant_report['ReportingYear'] == (Year - 3)) & (e_prtr_pollutant_report['CountryName'] == country)]




e_prtr_meta_domain_year = e_prtr_meta_croped.loc[e_prtr_meta_croped['PollutantReleaseAndTransferReportID'] == PRTR_ID.iloc[0]]

e_prtr_pollutants_full = pd.read_csv(InFolder + "/E-PRTR/dbo.PUBLISH_POLLUTANTRELEASE.csv", sep=',', encoding = "ISO-8859-1")

e_prtr_pollutants = e_prtr_pollutants_full.loc[e_prtr_pollutants_full['ReleaseMediumCode'] == "AIR"]
print(e_prtr_pollutants)
print(len(e_prtr_pollutants))

sums = []
e_prtr_point_gdf = pd.DataFrame(columns = ["snap","xcor","ycor","Hi","Vi","Ti","radi","CH4", "CO", "NH3", "NMVOC", "NOX", "PM10", "PM25", "SOX"])
e_prtr_point_uect = []





pollutants = ["CH4", "CO", "NH3", "NMVOC", "NOX", "PM10", "PM25", "SOX"]

for i in range(0, len(pollutants)):
    myData_temp = e_prtr_pollutants.loc[e_prtr_pollutants['PollutantCode'] == pollutants[i]]

    myData = pd.merge(e_prtr_meta_domain_year, myData_temp, on = ["FacilityReportID"], how='inner', left_index=True, right_index=False)

    sums.append(myData['TotalQuantity'].sum()/1000/1000)



    if sums[i] > 0:
        pts_df = myData[['Long', 'Lat', 'TotalQuantity', 'MainIASectorCode']]

        pts_gdf_wgs = gpd.GeoDataFrame(pts_df, geometry=gpd.points_from_xy(pts_df.Long, pts_df.Lat))

        pts_gdf_wgs.crs = crs_wgs

        pts_gdf_utm = pts_gdf_wgs.to_crs(crs_utm)

        pts_gdf_copy = pts_gdf_utm.copy()
        pts_gdf_copy['MainIASectorCode'] = pts_gdf_copy['MainIASectorCode'].astype(int)
        pts_gdf_copy['snap'] = 0
        pts_gdf_copy['xcor'] = pts_gdf_utm.geometry.x
        pts_gdf_copy['ycor'] = pts_gdf_utm.geometry.y
        pts_gdf_copy['Hi'] = int(-999)
        pts_gdf_copy['Vi'] = int(-999)
        pts_gdf_copy['Ti'] = int(-999)
        pts_gdf_copy['radi'] = int(-999)
        pts_gdf_copy[pollutants[i]] = pts_gdf_copy['TotalQuantity']
##        pts = pts_gdf_copy.drop(['Long', 'Lat', 'TotalQuantity', 'MainIASectorCode', 'geometry'], axis=1)
        pts = pts_gdf_copy.drop(['Long', 'Lat', 'TotalQuantity', 'geometry'], axis=1)

        e_prtr_point_uect.append(pts)

    else:

        pts_df =  pd.DataFrame(0, index=np.arange(1), columns=["snap","xcor","ycor","Hi","Vi","Ti","radi", pollutants[i]])
        e_prtr_point_uect.append(pts_df)


        
temp = pd.concat(e_prtr_point_uect)		

cams_data_list = []

cams_temp = temp.fillna(-999)

if len(cams_temp.loc[(cams_temp['xcor'] != 0) & (cams_temp['ycor'] != 0)]) != 0:
    print("E-PRTR point information EXIST for this domain")



    if (Year == 2017) & (country == "Greece"):

        cams_temp['snap'] = 0
        
        #*** E-PRTR == Energy sector ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 1] = 1
        cams_test = cams_temp.loc[cams_temp['snap'] == 1]
        cams_test['NOX'] = cams_test['NOX'] * 1.09
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.12
        cams_test['SOX'] = cams_test['SOX'] * 1.24
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1.21
        cams_test['PM10'] = cams_test['PM10'] * 1.24
        cams_test['CO'] = cams_test['CO'] * 1.10
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Production and processing of metals ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 2] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1.17
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.17
        cams_test['SOX'] = cams_test['SOX'] * 1.17
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1.17
        cams_test['PM10'] = cams_test['PM10'] * 1.17
        cams_test['CO'] = cams_test['CO'] * 1.17
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Mineral industry ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 3] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 0.63
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.08
        cams_test['SOX'] = cams_test['SOX'] * 1.17
        cams_test['NH3'] = cams_test['NH3'] * 1.27
        cams_test['PM25'] = cams_test['PM25'] * 0.99
        cams_test['PM10'] = cams_test['PM10'] * 0.99
        cams_test['CO'] = cams_test['CO'] * 1.18
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Chemical industry ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 4] = 4
        cams_test = cams_temp.loc[cams_temp['snap'] == 4]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 0.92
        cams_test['NH3'] = cams_test['NH3'] * 1.02
        cams_test['PM25'] = cams_test['PM25'] * 1.36
        cams_test['PM10'] = cams_test['PM10'] * 1.36
        cams_test['CO'] = cams_test['CO'] * 1.36
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Waste and waste water management ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 5] = 9
        cams_test = cams_temp.loc[cams_temp['snap'] == 9]
        cams_test['NOX'] = cams_test['NOX'] * 1.78
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.39
        cams_test['SOX'] = cams_test['SOX'] * 1.78
        cams_test['NH3'] = cams_test['NH3'] * 1.23
        cams_test['PM25'] = cams_test['PM25'] * 1.40
        cams_test['PM10'] = cams_test['PM10'] * 1.40
        cams_test['CO'] = cams_test['CO'] * 1.78
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Paper & wood production processing ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 6] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 1
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1
        cams_test['PM10'] = cams_test['PM10'] * 1
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Intensive livestock production & agriculture ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 7] = 10
        cams_test = cams_temp.loc[cams_temp['snap'] == 10]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 1
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1
        cams_test['PM10'] = cams_test['PM10'] * 1
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Animal and vegetable products from the food and beverage sector ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 8] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1.03
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.04
        cams_test['SOX'] = cams_test['SOX'] * 1.03
        cams_test['NH3'] = cams_test['NH3'] * 1.05
        cams_test['PM25'] = cams_test['PM25'] * 1.05
        cams_test['PM10'] = cams_test['PM10'] * 1.05
        cams_test['CO'] = cams_test['CO'] * 1.05
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Other activities ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 9] = 6
        cams_test = cams_temp.loc[cams_temp['snap'] == 6]
        cams_test['NOX'] = cams_test['NOX'] * 1.04
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.04
        cams_test['SOX'] = cams_test['SOX'] * 1.04
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1.04
        cams_test['PM10'] = cams_test['PM10'] * 1.04
        cams_test['CO'] = cams_test['CO'] * 1.04
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)





    elif (Year == 2018) & (country == "Greece"):

        cams_temp['snap'] = 0
        
        #*** E-PRTR == Energy sector ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 1] = 1
        cams_test = cams_temp.loc[cams_temp['snap'] == 1]
        cams_test['NOX'] = cams_test['NOX'] * 0.99
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.04
        cams_test['SOX'] = cams_test['SOX'] * 0.82
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1.06
        cams_test['PM10'] = cams_test['PM10'] * 1.09
        cams_test['CO'] = cams_test['CO'] * 1.09
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Production and processing of metals ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 2] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1.27
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.27
        cams_test['SOX'] = cams_test['SOX'] * 1.27
        cams_test['NH3'] = cams_test['NH3'] * 1.27
        cams_test['PM25'] = cams_test['PM25'] * 1.27
        cams_test['PM10'] = cams_test['PM10'] * 1.27
        cams_test['CO'] = cams_test['CO'] * 1.27
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Mineral industry ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 3] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 0.41
        cams_test['NMVOC'] = cams_test['NMVOC'] * 0.86
        cams_test['SOX'] = cams_test['SOX'] * 1.55
        cams_test['NH3'] = cams_test['NH3'] * 0.83
        cams_test['PM25'] = cams_test['PM25'] * 0.44
        cams_test['PM10'] = cams_test['PM10'] * 0.44
        cams_test['CO'] = cams_test['CO'] * 1.67
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Chemical industry ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 4] = 4
        cams_test = cams_temp.loc[cams_temp['snap'] == 4]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 0.73
        cams_test['NH3'] = cams_test['NH3'] * 0.83
        cams_test['PM25'] = cams_test['PM25'] * 2.16
        cams_test['PM10'] = cams_test['PM10'] * 2.16
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Waste and waste water management ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 5] = 9
        cams_test = cams_temp.loc[cams_temp['snap'] == 9]
        cams_test['NOX'] = cams_test['NOX'] * 25.95
        cams_test['NMVOC'] = cams_test['NMVOC'] * 13.26
        cams_test['SOX'] = cams_test['SOX'] * 25.86
        cams_test['NH3'] = cams_test['NH3'] * 1.55
        cams_test['PM25'] = cams_test['PM25'] * 13.18
        cams_test['PM10'] = cams_test['PM10'] * 13.18
        cams_test['CO'] = cams_test['CO'] * 25.82
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Paper & wood production processing ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 6] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 1
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1
        cams_test['PM10'] = cams_test['PM10'] * 1
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Intensive livestock production & agriculture ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 7] = 10
        cams_test = cams_temp.loc[cams_temp['snap'] == 10]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 1
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1
        cams_test['PM10'] = cams_test['PM10'] * 1
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Animal and vegetable products from the food and beverage sector ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 8] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1.03
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.08
        cams_test['SOX'] = cams_test['SOX'] * 0.96
        cams_test['NH3'] = cams_test['NH3'] * 1.10
        cams_test['PM25'] = cams_test['PM25'] * 1.08
        cams_test['PM10'] = cams_test['PM10'] * 1.08
        cams_test['CO'] = cams_test['CO'] * 1.08
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Other activities ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 9] = 6
        cams_test = cams_temp.loc[cams_temp['snap'] == 6]
        cams_test['NOX'] = cams_test['NOX'] * 1.25
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.25
        cams_test['SOX'] = cams_test['SOX'] * 1.25
        cams_test['NH3'] = cams_test['NH3'] * 1.25
        cams_test['PM25'] = cams_test['PM25'] * 1.25
        cams_test['PM10'] = cams_test['PM10'] * 1.25
        cams_test['CO'] = cams_test['CO'] * 1.25
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)





    elif (Year == 2019) & (country == "Greece"):

        cams_temp['snap'] = 0
        
        #*** E-PRTR == Energy sector ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 1] = 1
        cams_test = cams_temp.loc[cams_temp['snap'] == 1]
        cams_test['NOX'] = cams_test['NOX'] * 0.9
        cams_test['NMVOC'] = cams_test['NMVOC'] * 0.99
        cams_test['SOX'] = cams_test['SOX'] * 0.93
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 0.73
        cams_test['PM10'] = cams_test['PM10'] * 0.65
        cams_test['CO'] = cams_test['CO'] * 0.95
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Production and processing of metals ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 2] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1.17
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.17
        cams_test['SOX'] = cams_test['SOX'] * 1.17
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1.17
        cams_test['PM10'] = cams_test['PM10'] * 1.17
        cams_test['CO'] = cams_test['CO'] * 1.16
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Mineral industry ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 3] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 0.48
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.39
        cams_test['SOX'] = cams_test['SOX'] * 1.6
        cams_test['NH3'] = cams_test['NH3'] * 2.8
        cams_test['PM25'] = cams_test['PM25'] * 0.96
        cams_test['PM10'] = cams_test['PM10'] * 0.96
        cams_test['CO'] = cams_test['CO'] * 1.81
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Chemical industry ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 4] = 4
        cams_test = cams_temp.loc[cams_temp['snap'] == 4]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 1.18
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1.66
        cams_test['PM10'] = cams_test['PM10'] * 1.66
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Waste and waste water management ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 5] = 9
        cams_test = cams_temp.loc[cams_temp['snap'] == 9]
        cams_test['NOX'] = cams_test['NOX'] * 26.84
        cams_test['NMVOC'] = cams_test['NMVOC'] * 13.72
        cams_test['SOX'] = cams_test['SOX'] * 26.84
        cams_test['NH3'] = cams_test['NH3'] * 1.57
        cams_test['PM25'] = cams_test['PM25'] * 13.79
        cams_test['PM10'] = cams_test['PM10'] * 13.79
        cams_test['CO'] = cams_test['CO'] * 26.84
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Paper & wood production processing ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 6] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 1
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1
        cams_test['PM10'] = cams_test['PM10'] * 1
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Intensive livestock production & agriculture ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 7] = 10
        cams_test = cams_temp.loc[cams_temp['snap'] == 10]
        cams_test['NOX'] = cams_test['NOX'] * 1
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 1
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1
        cams_test['PM10'] = cams_test['PM10'] * 1
        cams_test['CO'] = cams_test['CO'] * 1
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Animal and vegetable products from the food and beverage sector ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 8] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_test['NOX'] = cams_test['NOX'] * 0.99
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1
        cams_test['SOX'] = cams_test['SOX'] * 0.96
        cams_test['NH3'] = cams_test['NH3'] * 1.02
        cams_test['PM25'] = cams_test['PM25'] * 1.01
        cams_test['PM10'] = cams_test['PM10'] * 1.01
        cams_test['CO'] = cams_test['CO'] * 1.01
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)
        
        #*** E-PRTR == Other activities ***#

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 9] = 6
        cams_test = cams_temp.loc[cams_temp['snap'] == 6]
        cams_test['NOX'] = cams_test['NOX'] * 1.23
        cams_test['NMVOC'] = cams_test['NMVOC'] * 1.23
        cams_test['SOX'] = cams_test['SOX'] * 1.23
        cams_test['NH3'] = cams_test['NH3'] * 1
        cams_test['PM25'] = cams_test['PM25'] * 1.23
        cams_test['PM10'] = cams_test['PM10'] * 1.23
        cams_test['CO'] = cams_test['CO'] * 1.23
        cams_test['CH4'] = cams_test['CH4'] * 1
        cams_data_list.append(cams_test)





    else:

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 1] = 1
        cams_test = cams_temp.loc[cams_temp['snap'] == 1]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 2] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 3] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 4] = 4
        cams_test = cams_temp.loc[cams_temp['snap'] == 4]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 5] = 9
        cams_test = cams_temp.loc[cams_temp['snap'] == 9]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 6] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 7] = 10
        cams_test = cams_temp.loc[cams_temp['snap'] == 10]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 8] = 3
        cams_test = cams_temp.loc[cams_temp['snap'] == 3]
        cams_data_list.append(cams_test)

        cams_temp['snap'].loc[cams_temp['MainIASectorCode'] == 9] = 6
        cams_test = cams_temp.loc[cams_temp['snap'] == 6]
        cams_data_list.append(cams_test)

        print(temp.snap.unique())



    cams_data = pd.concat(cams_data_list)
    print(cams_data)
    pts_no_nan = cams_data.fillna(-999)

    pts_final = pts_no_nan.loc[pts_no_nan['xcor'] != 0]

    final_grouped = pts_final.groupby(['snap','xcor','ycor'], as_index=False).sum()
    final_grouped[final_grouped < 0] = -999



    pts_new_gdf = gpd.GeoDataFrame(final_grouped, geometry=gpd.points_from_xy(x=final_grouped.xcor, y=final_grouped.ycor), crs = crs_utm)
    pts_new_gdf.to_file(OutFolder + "cams_points.shp", driver="ESRI Shapefile")

    pts_new = final_grouped[['snap', 'xcor', 'ycor', 'Hi', 'Vi', 'Ti', 'radi', 'NOX', 'NMVOC', 'CO', 'SOX', 'NH3', 'PM25', 'PM10']]
    pts_new.columns = ['snap','xcor','ycor','Hi','Vi','Ti','radi','NOx','NMVOC','CO','SO2','NH3','PM2.5','PM10']

    pts_new[pts_new < 0] = -999
    print(pts_new)
    print(pts_new.snap.unique())

else:
    print("There are NOOOOOOOOO E-PRTR point information for this domain")


    cams_data = e_prtr_point_gdf.copy()

    pts_no_nan = cams_data.fillna(-999)

    pts_final = pts_no_nan.loc[pts_no_nan['xcor'] != 0]

    final_grouped = pts_final.drop(['CH4'], axis=1)
    final_grouped[final_grouped < 0] = -999
    pts_new_gdf = gpd.GeoDataFrame(final_grouped, crs = crs_utm)

    pts_new = final_grouped.copy()
    pts_new.columns = ['snap','xcor','ycor','Hi','Vi','Ti','radi','NOx','NMVOC','CO','SO2','NH3','PM2.5','PM10']

    pts_new[pts_new < 0] = -999
    print(pts_new)
    print(pts_new.snap.unique())


pts_new.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_point_sources_" + str(Year) + ".csv", sep = ',', mode='w', header=True, index=False)
















#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
















#**************************************************************************************************************
#**************************************************************************************************************
#*************************************** Calculate Area Source Emissions ***************************************#
#**************************************************************************************************************
#**************************************************************************************************************
















###################### EMISSIONS:

##emissions = InFolder + "/cams_v3_1/"
emissions = "X:/NOA/SMURBS/Regions/Hamburg/cams_v3_1/"

###################### choose an emissions downscaling option:
###### ~ 1 "top_down_proxy" distribution of total grid emission value using a normalized proxy
###### ~ 2 "coarse_proxy" distribution applies a normalized proxy grid that was created based on the coarse grid-cells of the original CAMS resolution

proxy_method = "coarse_cells_proxy"







if os.path.isfile(OutFolder + "cams_data_area.shp"):
    cams_temp = gpd.read_file(OutFolder + "cams_data_area.shp")
else:
    GNFR = ["A_PublicPower", "B_Industry", "C_OtherStationaryComb", "D_Fugitives", "E_Solvents", "F_RoadTransport", "G_Shipping", "H_Aviation", "I_OffRoad", "J_Waste", "K_AgriLivestock", "L_AgriOther", "SumAllSectors"]

    GNFR_raster_list = []
    for file in os.listdir(emissions):
        if file.endswith(".nc"):
            path = os.path.join(emissions, file)
            ds = xr.open_dataset(path)
            df = ds.to_dataframe().reset_index(level=['lat', 'lon', 'time'])

            for i in range(0, len(GNFR)):
                cams_temp = df[['lat', 'lon', GNFR[i]]]
                cams_temp_gdf = gpd.GeoDataFrame(cams_temp, geometry=gpd.points_from_xy(cams_temp.lon, cams_temp.lat))
                cams_crop = cl.clip_shp(cams_temp_gdf, cams_grid)
                
                if 'ch4' in file:
                    cams_crop['CH4'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)
                
                elif 'co' in file:
                    cams_crop['CO'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)

                elif 'nh3' in file:
                    cams_crop['NH3'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)

                elif 'nmvoc' in file:
                    cams_crop['NMVOC'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)

                elif 'nox' in file:
                    cams_crop['NOX'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)

                elif 'pm2_5' in file:
                    cams_crop['PM2_5'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)

                elif 'pm10' in file:
                    cams_crop['PM10'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)

                elif 'so2' in file:
                    cams_crop['SO2'] = cams_crop[GNFR[i]]
                    cams_crop['GNFR'] = GNFR[i]
                    cams = cams_crop.drop([GNFR[i]], axis=1)
                    print(cams)
                    GNFR_raster_list.append(cams)

                else:
                    print("nope")
                
    cams_temp = pd.concat(GNFR_raster_list)
    print(cams_temp)
    cams_temp.to_file(OutFolder + "cams_data_area.shp", driver="ESRI Shapefile")

print("**********")

cams_data_list = []





### GNFR to SNAP by YEAR ( + scaling factors)





if (Year == 2017) & (Region == 'Athens'):
    
    cams_temp['SNAP'] = 0
    
    #*** GNFR == 'A_PublicPower' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'A_PublicPower'] = 1
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 1]
    cams_test['NOX'] = cams_test['NOX'] * 1.09
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.12
    cams_test['SO2'] = cams_test['SO2'] * 1.24
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.21
    cams_test['PM10'] = cams_test['PM10'] * 1.24
    cams_test['CO'] = cams_test['CO'] * 1.10
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'B_Industry'  - SNAP == 3 ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'B_Industry'] = 34
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 34]

    cams_snap3 =cams_test.copy()
    cams_snap3.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.8
    cams_snap3['SNAP'] = 3
    cams_snap3['NOX'] = cams_snap3['NOX'] * 1.14
    cams_snap3['NMVOC'] = cams_snap3['NMVOC'] * 1.27
    cams_snap3['SO2'] = cams_snap3['SO2'] * 1.4
    cams_snap3['NH3'] = cams_snap3['NH3'] * 1.23
    cams_snap3['PM2_5'] = cams_snap3['PM2_5'] * 1.04
    cams_snap3['PM10'] = cams_snap3['PM10'] * 1.05
    cams_snap3['CO'] = cams_snap3['CO'] * 1.12
    cams_snap3['CH4'] = cams_snap3['CH4'] * 1
    cams_data_list.append(cams_snap3)
    
    #*** GNFR == 'B_Industry'  - SNAP == 4 ***#

    cams_snap4 =cams_test.copy()
    cams_snap4.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.2
    cams_snap4['SNAP'] = 4
    cams_snap4['NOX'] = cams_snap4['NOX'] * 1.14
    cams_snap4['NMVOC'] = cams_snap4['NMVOC'] * 1.27
    cams_snap4['SO2'] = cams_snap4['SO2'] * 1.4
    cams_snap4['NH3'] = cams_snap4['NH3'] * 1.23
    cams_snap4['PM2_5'] = cams_snap4['PM2_5'] * 1.04
    cams_snap4['PM10'] = cams_snap4['PM10'] * 1.05
    cams_snap4['CO'] = cams_snap4['CO'] * 1.12
    cams_snap4['CH4'] = cams_snap4['CH4'] * 1
    cams_data_list.append(cams_snap4)
    
    #*** GNFR == 'C_OtherStationaryComb' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'C_OtherStationaryComb'] = 2
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 2]
    cams_test['NOX'] = cams_test['NOX'] * 1.03
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.12
    cams_test['SO2'] = cams_test['SO2'] * 0.99
    cams_test['NH3'] = cams_test['NH3'] * 1.12
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.12
    cams_test['PM10'] = cams_test['PM10'] * 1.12
    cams_test['CO'] = cams_test['CO'] * 1.11
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
        
    #*** GNFR == 'D_Fugitives' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'D_Fugitives'] = 5
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 5]
    cams_test['NOX'] = cams_test['NOX'] * 1
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.95
    cams_test['SO2'] = cams_test['SO2'] * 1
    cams_test['NH3'] = cams_test['NH3'] * 1.04
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.16
    cams_test['PM10'] = cams_test['PM10'] * 1.16
    cams_test['CO'] = cams_test['CO'] * 1.04
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'E_Solvents' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'E_Solvents'] = 6
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 6]
    cams_test['NOX'] = cams_test['NOX'] * 0.58
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.94
    cams_test['SO2'] = cams_test['SO2'] * 0.99
    cams_test['NH3'] = cams_test['NH3'] * 0.58
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.79
    cams_test['PM10'] = cams_test['PM10'] * 0.79
    cams_test['CO'] = cams_test['CO'] * 0.58
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'F_RoadTransport' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'F_RoadTransport'] = 7
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 7]
    cams_test['NOX'] = cams_test['NOX'] * 1.06
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.92
    cams_test['SO2'] = cams_test['SO2'] * 0.96
    cams_test['NH3'] = cams_test['NH3'] * 1.01
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.99
    cams_test['PM10'] = cams_test['PM10'] * 0.99
    cams_test['CO'] = cams_test['CO'] * 1
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'G_Shipping' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'G_Shipping'] = 8
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 8]
    cams_test['NOX'] = cams_test['NOX'] * 1.02
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.02
    cams_test['SO2'] = cams_test['SO2'] * 1.06
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.02
    cams_test['PM10'] = cams_test['PM10'] * 1.02
    cams_test['CO'] = cams_test['CO'] * 1.02
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'J_Waste' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'J_Waste'] = 9
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 9]
    cams_test['NOX'] = cams_test['NOX'] * 1.78
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.39
    cams_test['SO2'] = cams_test['SO2'] * 1.78
    cams_test['NH3'] = cams_test['NH3'] * 1.23
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.4
    cams_test['PM10'] = cams_test['PM10'] * 1.4
    cams_test['CO'] = cams_test['CO'] * 1.78
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'K_AgriLivestock' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'K_AgriLivestock'] = 10
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 10]
    cams_test['NOX'] = cams_test['NOX'] * 1
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1
    cams_test['SO2'] = cams_test['SO2'] * 1
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1
    cams_test['PM10'] = cams_test['PM10'] * 1
    cams_test['CO'] = cams_test['CO'] * 1
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'L_AgriOther' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'L_AgriOther'] = 10
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 10]
    cams_test['NOX'] = cams_test['NOX'] * 1.01
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.99
    cams_test['SO2'] = cams_test['SO2'] * 1.01
    cams_test['NH3'] = cams_test['NH3'] * 0.99
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.99
    cams_test['PM10'] = cams_test['PM10'] * 0.99
    cams_test['CO'] = cams_test['CO'] * 1.01
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'H_Aviation' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'H_Aviation'] = 11
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 11]
    cams_test['NOX'] = cams_test['NOX'] * 1.03
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.99
    cams_test['SO2'] = cams_test['SO2'] * 1.04
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.04
    cams_test['PM10'] = cams_test['PM10'] * 1.04
    cams_test['CO'] = cams_test['CO'] * 1.05
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'I_OffRoad' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'I_OffRoad'] = 12
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 12]
    cams_test['NOX'] = cams_test['NOX'] * 0.96
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.18
    cams_test['SO2'] = cams_test['SO2'] * 0.98
    cams_test['NH3'] = cams_test['NH3'] * 0.99
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1
    cams_test['PM10'] = cams_test['PM10'] * 1
    cams_test['CO'] = cams_test['CO'] * 1.21
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'SumAllSectors' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'SumAllSectors'] = 0
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 0]
    cams_data_list.append(cams_test)

    ##print(cams_test.loc[cams_test['GNFR'] == 'A_PublicPower'])





elif (Year == 2018) & (Region == 'Athens'):
    
    cams_temp['SNAP'] = 0
    
    #*** GNFR == 'A_PublicPower' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'A_PublicPower'] = 1
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 1]
    cams_test['NOX'] = cams_test['NOX'] * 0.99
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.04
    cams_test['SO2'] = cams_test['SO2'] * 0.82
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.06
    cams_test['PM10'] = cams_test['PM10'] * 1.09
    cams_test['CO'] = cams_test['CO'] * 1.09
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'B_Industry'  - SNAP == 3 ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'B_Industry'] = 34
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 34]

    cams_snap3 =cams_test.copy()
    cams_snap3.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.8
    cams_snap3['SNAP'] = 3
    cams_snap3['NOX'] = cams_snap3['NOX'] * 1.09
    cams_snap3['NMVOC'] = cams_snap3['NMVOC'] * 1.34
    cams_snap3['SO2'] = cams_snap3['SO2'] * 1.14
    cams_snap3['NH3'] = cams_snap3['NH3'] * 1.07
    cams_snap3['PM2_5'] = cams_snap3['PM2_5'] * 1
    cams_snap3['PM10'] = cams_snap3['PM10'] * 1
    cams_snap3['CO'] = cams_snap3['CO'] * 1.16
    cams_snap3['CH4'] = cams_snap3['CH4'] * 1
    cams_data_list.append(cams_snap3)
    
    #*** GNFR == 'B_Industry'  - SNAP == 4 ***#

    cams_snap4 =cams_test.copy()
    cams_snap4.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.2
    cams_snap4['SNAP'] = 4
    cams_snap4['NOX'] = cams_snap4['NOX'] * 1.09
    cams_snap4['NMVOC'] = cams_snap4['NMVOC'] * 1.34
    cams_snap4['SO2'] = cams_snap4['SO2'] * 1.14
    cams_snap4['NH3'] = cams_snap4['NH3'] * 1.07
    cams_snap4['PM2_5'] = cams_snap4['PM2_5'] * 1
    cams_snap4['PM10'] = cams_snap4['PM10'] * 1
    cams_snap4['CO'] = cams_snap4['CO'] * 1.16
    cams_snap4['CH4'] = cams_snap4['CH4'] * 1
    cams_data_list.append(cams_snap4)
    
    #*** GNFR == 'C_OtherStationaryComb' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'C_OtherStationaryComb'] = 2
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 2]
    cams_test['NOX'] = cams_test['NOX'] * 1.05
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.67
    cams_test['SO2'] = cams_test['SO2'] * 0.89
    cams_test['NH3'] = cams_test['NH3'] * 1.78
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.73
    cams_test['PM10'] = cams_test['PM10'] * 1.73
    cams_test['CO'] = cams_test['CO'] * 1.91
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
        
    #*** GNFR == 'D_Fugitives' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'D_Fugitives'] = 5
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 5]
    cams_test['NOX'] = cams_test['NOX'] * 1
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.05
    cams_test['SO2'] = cams_test['SO2'] * 1
    cams_test['NH3'] = cams_test['NH3'] * 1.05
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.12
    cams_test['PM10'] = cams_test['PM10'] * 1.12
    cams_test['CO'] = cams_test['CO'] * 1.05
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'E_Solvents' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'E_Solvents'] = 6
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 6]
    cams_test['NOX'] = cams_test['NOX'] * 0.71
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.71
    cams_test['SO2'] = cams_test['SO2'] * 1.34
    cams_test['NH3'] = cams_test['NH3'] * 0.71
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.78
    cams_test['PM10'] = cams_test['PM10'] * 0.78
    cams_test['CO'] = cams_test['CO'] * 0.71
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'F_RoadTransport' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'F_RoadTransport'] = 7
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 7]
    cams_test['NOX'] = cams_test['NOX'] * 1.20
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.86
    cams_test['SO2'] = cams_test['SO2'] * 0.64
    cams_test['NH3'] = cams_test['NH3'] * 1.03
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.98
    cams_test['PM10'] = cams_test['PM10'] * 0.98
    cams_test['CO'] = cams_test['CO'] * 0.97
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'G_Shipping' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'G_Shipping'] = 8
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 8]
    cams_test['NOX'] = cams_test['NOX'] * 0.89
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.10
    cams_test['SO2'] = cams_test['SO2'] * 131.28
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.08
    cams_test['PM10'] = cams_test['PM10'] * 1.10
    cams_test['CO'] = cams_test['CO'] * 1.10
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'J_Waste' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'J_Waste'] = 9
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 9]
    cams_test['NOX'] = cams_test['NOX'] * 25.95
    cams_test['NMVOC'] = cams_test['NMVOC'] * 13.26
    cams_test['SO2'] = cams_test['SO2'] * 25.86
    cams_test['NH3'] = cams_test['NH3'] * 1.55
    cams_test['PM2_5'] = cams_test['PM2_5'] * 13.18
    cams_test['PM10'] = cams_test['PM10'] * 13.18
    cams_test['CO'] = cams_test['CO'] * 25.82
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'K_AgriLivestock' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'K_AgriLivestock'] = 10
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 10]
    cams_test['NOX'] = cams_test['NOX'] * 0.89
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.88
    cams_test['SO2'] = cams_test['SO2'] * 1
    cams_test['NH3'] = cams_test['NH3'] * 0.96
    cams_test['PM2_5'] = cams_test['PM2_5'] * 3
    cams_test['PM10'] = cams_test['PM10'] * 2.98
    cams_test['CO'] = cams_test['CO'] * 1
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'L_AgriOther' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'L_AgriOther'] = 10
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 10]
    cams_test['NOX'] = cams_test['NOX'] * 0.97
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.92
    cams_test['SO2'] = cams_test['SO2'] * 0.88
    cams_test['NH3'] = cams_test['NH3'] * 1.08
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.96
    cams_test['PM10'] = cams_test['PM10'] * 0.92
    cams_test['CO'] = cams_test['CO'] * 0.88
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'H_Aviation' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'H_Aviation'] = 11
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 11]
    cams_test['NOX'] = cams_test['NOX'] * 1.10
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.02
    cams_test['SO2'] = cams_test['SO2'] * 1.16
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.18
    cams_test['PM10'] = cams_test['PM10'] * 1.18
    cams_test['CO'] = cams_test['CO'] * 1.21
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'I_OffRoad' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'I_OffRoad'] = 12
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 12]
    cams_test['NOX'] = cams_test['NOX'] * 0.91
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.95
    cams_test['SO2'] = cams_test['SO2'] * 0.9
    cams_test['NH3'] = cams_test['NH3'] * 1.01
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1
    cams_test['PM10'] = cams_test['PM10'] * 1
    cams_test['CO'] = cams_test['CO'] * 0.86
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'SumAllSectors' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'SumAllSectors'] = 0
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 0]
    cams_data_list.append(cams_test)






elif (Year == 2019) & (Region == 'Athens'):
    
    cams_temp['SNAP'] = 0
    
    #*** GNFR == 'A_PublicPower' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'A_PublicPower'] = 1
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 1]
    cams_test['NOX'] = cams_test['NOX'] * 0.9
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.99
    cams_test['SO2'] = cams_test['SO2'] * 0.93
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.73
    cams_test['PM10'] = cams_test['PM10'] * 0.65
    cams_test['CO'] = cams_test['CO'] * 0.95
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'B_Industry'  - SNAP == 3 ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'B_Industry'] = 34
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 34]

    cams_snap3 =cams_test.copy()
    cams_snap3.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.8
    cams_snap3['SNAP'] = 3
    cams_snap3['NOX'] = cams_snap3['NOX'] * 1.01
    cams_snap3['NMVOC'] = cams_snap3['NMVOC'] * 1.41
    cams_snap3['SO2'] = cams_snap3['SO2'] * 1.04
    cams_snap3['NH3'] = cams_snap3['NH3'] * 1.66
    cams_snap3['PM2_5'] = cams_snap3['PM2_5'] * 1.07
    cams_snap3['PM10'] = cams_snap3['PM10'] * 1.07
    cams_snap3['CO'] = cams_snap3['CO'] * 1.19
    cams_snap3['CH4'] = cams_snap3['CH4'] * 1
    cams_data_list.append(cams_snap3)
    
    #*** GNFR == 'B_Industry'  - SNAP == 4 ***#

    cams_snap4 =cams_test.copy()
    cams_snap4.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.2
    cams_snap4['SNAP'] = 4
    cams_snap4['NOX'] = cams_snap4['NOX'] * 1.01
    cams_snap4['NMVOC'] = cams_snap4['NMVOC'] * 1.41
    cams_snap4['SO2'] = cams_snap4['SO2'] * 1.04
    cams_snap4['NH3'] = cams_snap4['NH3'] * 1.66
    cams_snap4['PM2_5'] = cams_snap4['PM2_5'] * 1.07
    cams_snap4['PM10'] = cams_snap4['PM10'] * 1.07
    cams_snap4['CO'] = cams_snap4['CO'] * 1.19
    cams_snap4['CH4'] = cams_snap4['CH4'] * 1
    cams_data_list.append(cams_snap4)
    
    #*** GNFR == 'C_OtherStationaryComb' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'C_OtherStationaryComb'] = 2
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 2]
    cams_test['NOX'] = cams_test['NOX'] * 1.13
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.62
    cams_test['SO2'] = cams_test['SO2'] * 0.97
    cams_test['NH3'] = cams_test['NH3'] * 1.72
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.67
    cams_test['PM10'] = cams_test['PM10'] * 1.67
    cams_test['CO'] = cams_test['CO'] * 1.9
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
        
    #*** GNFR == 'D_Fugitives' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'D_Fugitives'] = 5
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 5]
    cams_test['NOX'] = cams_test['NOX'] * 1
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.88
    cams_test['SO2'] = cams_test['SO2'] * 1
    cams_test['NH3'] = cams_test['NH3'] * 0.99
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.84
    cams_test['PM10'] = cams_test['PM10'] * 0.84
    cams_test['CO'] = cams_test['CO'] * 0.99
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'E_Solvents' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'E_Solvents'] = 6
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 6]
    cams_test['NOX'] = cams_test['NOX'] * 1.06
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.93
    cams_test['SO2'] = cams_test['SO2'] * 1.67
    cams_test['NH3'] = cams_test['NH3'] * 1.05
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.93
    cams_test['PM10'] = cams_test['PM10'] * 0.93
    cams_test['CO'] = cams_test['CO'] * 1.05
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'F_RoadTransport' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'F_RoadTransport'] = 7
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 7]
    cams_test['NOX'] = cams_test['NOX'] * 1.05
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.86
    cams_test['SO2'] = cams_test['SO2'] * 0.66
    cams_test['NH3'] = cams_test['NH3'] * 1.09
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.03
    cams_test['PM10'] = cams_test['PM10'] * 1.03
    cams_test['CO'] = cams_test['CO'] * 0.98
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'G_Shipping' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'G_Shipping'] = 8
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 8]
    cams_test['NOX'] = cams_test['NOX'] * 0.9
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.15
    cams_test['SO2'] = cams_test['SO2'] * 139.86
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.07
    cams_test['PM10'] = cams_test['PM10'] * 0.97
    cams_test['CO'] = cams_test['CO'] * 1.15
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'J_Waste' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'J_Waste'] = 9
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 9]
    cams_test['NOX'] = cams_test['NOX'] * 26.84
    cams_test['NMVOC'] = cams_test['NMVOC'] * 13.72
    cams_test['SO2'] = cams_test['SO2'] * 26.84
    cams_test['NH3'] = cams_test['NH3'] * 1.57
    cams_test['PM2_5'] = cams_test['PM2_5'] * 13.79
    cams_test['PM10'] = cams_test['PM10'] * 13.79
    cams_test['CO'] = cams_test['CO'] * 26.84
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'K_AgriLivestock' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'K_AgriLivestock'] = 10
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 10]
    cams_test['NOX'] = cams_test['NOX'] * 0.9
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.92
    cams_test['SO2'] = cams_test['SO2'] * 1
    cams_test['NH3'] = cams_test['NH3'] * 0.91
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.92
    cams_test['PM10'] = cams_test['PM10'] * 0.92
    cams_test['CO'] = cams_test['CO'] * 1
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'L_AgriOther' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'L_AgriOther'] = 10
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 10]
    cams_test['NOX'] = cams_test['NOX'] * 0.99
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.93
    cams_test['SO2'] = cams_test['SO2'] * 0.93
    cams_test['NH3'] = cams_test['NH3'] * 0.98
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.96
    cams_test['PM10'] = cams_test['PM10'] * 0.93
    cams_test['CO'] = cams_test['CO'] * 0.93
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'H_Aviation' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'H_Aviation'] = 11
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 11]
    cams_test['NOX'] = cams_test['NOX'] * 1.10
    cams_test['NMVOC'] = cams_test['NMVOC'] * 1.03
    cams_test['SO2'] = cams_test['SO2'] * 1.17
    cams_test['NH3'] = cams_test['NH3'] * 1
    cams_test['PM2_5'] = cams_test['PM2_5'] * 1.16
    cams_test['PM10'] = cams_test['PM10'] * 1.16
    cams_test['CO'] = cams_test['CO'] * 1.24
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'I_OffRoad' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'I_OffRoad'] = 12
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 12]
    cams_test['NOX'] = cams_test['NOX'] * 0.61
    cams_test['NMVOC'] = cams_test['NMVOC'] * 0.82
    cams_test['SO2'] = cams_test['SO2'] * 0.6
    cams_test['NH3'] = cams_test['NH3'] * 0.77
    cams_test['PM2_5'] = cams_test['PM2_5'] * 0.74
    cams_test['PM10'] = cams_test['PM10'] * 0.74
    cams_test['CO'] = cams_test['CO'] * 0.67
    cams_test['CH4'] = cams_test['CH4'] * 1
    cams_data_list.append(cams_test)
    
    #*** GNFR == 'SumAllSectors' ***#

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'SumAllSectors'] = 0
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 0]
    cams_data_list.append(cams_test)






else:
    
    cams_temp['SNAP'] = 0

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'A_PublicPower'] = 1
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 1]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'C_OtherStationaryComb'] = 2
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 2]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'B_Industry'] = 34
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 34]

    cams_snap3 =cams_test.copy()
    cams_snap3.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.8
    cams_snap3['SNAP'] = 3
    cams_data_list.append(cams_snap3)

    cams_snap4 =cams_test.copy()
    cams_snap4.loc[:, ['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']] *=0.2
    cams_snap4['SNAP'] = 4
    cams_data_list.append(cams_snap4)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'D_Fugitives'] = 5
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 5]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'E_Solvents'] = 6
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 6]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'F_RoadTransport'] = 7
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 7]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'G_Shipping'] = 8
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 8]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'J_Waste'] = 9
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 9]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[(cams_temp['GNFR'] == 'K_AgriLivestock') | (cams_temp['GNFR'] == 'L_AgriOther')] = 10
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 10]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'H_Aviation'] = 11
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 11]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'I_OffRoad'] = 12
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 12]
    cams_data_list.append(cams_test)

    cams_temp['SNAP'].loc[cams_temp['GNFR'] == 'SumAllSectors'] = 0
    cams_test = cams_temp.loc[cams_temp['SNAP'] == 0]
    cams_data_list.append(cams_test)



cams_data = pd.concat(cams_data_list)
print(cams_data['CO'].loc[cams_data['GNFR'] == 'J_Waste'].unique())

##filtered_data = cams_data.groupby(['lat', 'lon', 'GNFR'], as_index=False).sum()
filtered_data = cams_data.groupby(['lat', 'lon', 'SNAP'], as_index=False).sum()


### Multiply with 1000000000 to get from Tg to kg
filtered_data.loc[:,'CH4'] *= 1000000000
filtered_data.loc[:,'NOX'] *= 1000000000
filtered_data.loc[:,'NMVOC'] *= 1000000000
filtered_data.loc[:,'CO'] *= 1000000000
filtered_data.loc[:,'SO2'] *= 1000000000
filtered_data.loc[:,'NH3'] *= 1000000000
filtered_data.loc[:,'PM2_5'] *= 1000000000
filtered_data.loc[:,'PM10'] *= 1000000000
print(filtered_data)


df = filtered_data.fillna(0)
gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat))
gdf.crs = crs_wgs

projected_data = gdf.to_crs(crs_utm)
projected_data.to_file(OutFolder + "projected_data_area.shp", driver="ESRI Shapefile")
projected_data = gpd.read_file(OutFolder + "projected_data_area.shp")

croped_data = cl.clip_shp(projected_data, grid)
croped_data = projected_data.loc[projected_data['SNAP'] != 0]

### Create cols with xcor & ycor from projected data
croped_data['xcor'] = croped_data['geometry'].x
croped_data['ycor'] = croped_data['geometry'].y
croped_data['GNFR'] = 'test'
croped_data['GNFR'].loc[croped_data['SNAP'] == 1] = 'A_PublicPower'
croped_data['GNFR'].loc[croped_data['SNAP'] == 2] = 'C_OtherStationaryComb'
croped_data['GNFR'].loc[croped_data['SNAP'] == 3] = 'B_Industry'
croped_data['GNFR'].loc[croped_data['SNAP'] == 4] = 'B_Industry'
croped_data['GNFR'].loc[croped_data['SNAP'] == 5] = 'D_Fugitives'
croped_data['GNFR'].loc[croped_data['SNAP'] == 6] = 'E_Solvents'
croped_data['GNFR'].loc[croped_data['SNAP'] == 7] = 'F_RoadTransport'
croped_data['GNFR'].loc[croped_data['SNAP'] == 8] = 'G_Shipping'
croped_data['GNFR'].loc[croped_data['SNAP'] == 9] = 'J_Waste'
croped_data['GNFR'].loc[croped_data['SNAP'] == 10] = 'K_AgriLivestock'
croped_data['GNFR'].loc[croped_data['SNAP'] == 11] = 'H_Aviation'
croped_data['GNFR'].loc[croped_data['SNAP'] == 12] = 'I_OffRoad'

print(croped_data)

croped_data.to_file(OutFolder + "croped_data_area.shp", driver="ESRI Shapefile")
croped_data = gpd.read_file(OutFolder + "croped_data_area.shp")


col_data = croped_data.drop(['lon', 'lat', 'geometry'], axis=1)
print(col_data)
col_names = col_data.columns.unique()
print(col_data.SNAP.unique())


SNAP_sectors = col_data.SNAP.unique()






#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








#**************************************************************************************************************
#*************************************** Raster to Points Geodataframe ***************************************#
#**************************************************************************************************************



#*** Population Density ***#


pop_raster = rasterio.open(Proxy_Folder + "/pop_proj.tif")
print(pop_raster.crs)
print(pop_raster.count)

pop_pr_x_list = []
pop_pr_y_list = []
pop_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = pop_raster.index(x,y)
    val = pop_raster.read(1)[row,col]
    pop_pr_x_list.append(x)
    pop_pr_y_list.append(y)
    pop_pr_val_list.append(val)
    
pop_norm_post_proc_df = pd.DataFrame(list(zip(pop_pr_x_list, pop_pr_y_list, pop_pr_val_list)), columns =['xcor', 'ycor', 'pop_val']) 
pop_raster_to_points_gdf = gpd.GeoDataFrame(pop_norm_post_proc_df, geometry=gpd.points_from_xy(x=pop_norm_post_proc_df.xcor, y=pop_norm_post_proc_df.ycor), crs = crs_utm)
pop_raster_to_points_gdf['boolean'] = 0
pop_raster_to_points_gdf['boolean'].loc[pop_raster_to_points_gdf['pop_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
pop_grid_gdf = gpd.sjoin(grid, pop_raster_to_points_gdf, how='inner', op='contains')
pop_grid_norm_gdf = pop_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#








#*** Industry ***#


ind_raster = rasterio.open(Proxy_Folder + "/lu_industry.tif")
print(ind_raster.crs)
print(ind_raster.count)

ind_pr_x_list = []
ind_pr_y_list = []
ind_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = ind_raster.index(x,y)
    val = ind_raster.read(1)[row,col]
    ind_pr_x_list.append(x)
    ind_pr_y_list.append(y)
    ind_pr_val_list.append(val)
    
ind_norm_post_proc_df = pd.DataFrame(list(zip(ind_pr_x_list, ind_pr_y_list, ind_pr_val_list)), columns =['xcor', 'ycor', 'ind_val']) 
ind_raster_to_points_gdf = gpd.GeoDataFrame(ind_norm_post_proc_df, geometry=gpd.points_from_xy(x=ind_norm_post_proc_df.xcor, y=ind_norm_post_proc_df.ycor), crs = crs_utm)
ind_raster_to_points_gdf['boolean'] = 0
ind_raster_to_points_gdf['boolean'].loc[ind_raster_to_points_gdf['ind_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
ind_grid_gdf = gpd.sjoin(grid, ind_raster_to_points_gdf, how='inner', op='contains')
ind_grid_norm_gdf = ind_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#








#*** Agriculture ***#


agr_raster = rasterio.open(Proxy_Folder + "/lu_agriculture.tif")
print(agr_raster.crs)
print(agr_raster.count)

agr_pr_x_list = []
agr_pr_y_list = []
agr_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = agr_raster.index(x,y)
    val = agr_raster.read(1)[row,col]
    agr_pr_x_list.append(x)
    agr_pr_y_list.append(y)
    agr_pr_val_list.append(val)
    
agr_norm_post_proc_df = pd.DataFrame(list(zip(agr_pr_x_list, agr_pr_y_list, agr_pr_val_list)), columns =['xcor', 'ycor', 'agr_val']) 
agr_raster_to_points_gdf = gpd.GeoDataFrame(agr_norm_post_proc_df, geometry=gpd.points_from_xy(x=agr_norm_post_proc_df.xcor, y=agr_norm_post_proc_df.ycor), crs = crs_utm)
agr_raster_to_points_gdf['boolean'] = 0
agr_raster_to_points_gdf['boolean'].loc[agr_raster_to_points_gdf['agr_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
agr_grid_gdf = gpd.sjoin(grid, agr_raster_to_points_gdf, how='inner', op='contains')
agr_grid_norm_gdf = agr_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#








#*** Non Road Mobility Sources ***#


snap8_raster = rasterio.open(Proxy_Folder + "/lu_offroad.tif")
print(snap8_raster.crs)
print(snap8_raster.count)

snap8_pr_x_list = []
snap8_pr_y_list = []
snap8_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = snap8_raster.index(x,y)
    val = snap8_raster.read(1)[row,col]
    snap8_pr_x_list.append(x)
    snap8_pr_y_list.append(y)
    snap8_pr_val_list.append(val)
    
snap8_norm_post_proc_df = pd.DataFrame(list(zip(snap8_pr_x_list, snap8_pr_y_list, snap8_pr_val_list)), columns =['xcor', 'ycor', 'offroad_val']) 
snap8_raster_to_points_gdf = gpd.GeoDataFrame(snap8_norm_post_proc_df, geometry=gpd.points_from_xy(x=snap8_norm_post_proc_df.xcor, y=snap8_norm_post_proc_df.ycor), crs = crs_utm)
snap8_raster_to_points_gdf['boolean'] = 0
snap8_raster_to_points_gdf['boolean'].loc[snap8_raster_to_points_gdf['offroad_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
snap8_grid_gdf = gpd.sjoin(grid, snap8_raster_to_points_gdf, how='inner', op='contains')
snap8_grid_norm_gdf = snap8_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#








#*** Waste ***#


waste_raster = rasterio.open(Proxy_Folder + "/lu_waste.tif")
print(waste_raster.crs)
print(waste_raster.count)

waste_pr_x_list = []
waste_pr_y_list = []
waste_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = waste_raster.index(x,y)
    val = waste_raster.read(1)[row,col]
    waste_pr_x_list.append(x)
    waste_pr_y_list.append(y)
    waste_pr_val_list.append(val)
    
waste_norm_post_proc_df = pd.DataFrame(list(zip(waste_pr_x_list, waste_pr_y_list, waste_pr_val_list)), columns =['xcor', 'ycor', 'waste_val']) 
waste_raster_to_points_gdf = gpd.GeoDataFrame(waste_norm_post_proc_df, geometry=gpd.points_from_xy(x=waste_norm_post_proc_df.xcor, y=waste_norm_post_proc_df.ycor), crs = crs_utm)
waste_raster_to_points_gdf['boolean'] = 0
waste_raster_to_points_gdf['boolean'].loc[waste_raster_to_points_gdf['waste_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
waste_grid_gdf = gpd.sjoin(grid, waste_raster_to_points_gdf, how='inner', op='contains')
waste_grid_norm_gdf = waste_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#








#*** Snap 34 ***#


snap34_raster = rasterio.open(Proxy_Folder + "/lu_snap34.tif")
print(snap34_raster.crs)
print(snap34_raster.count)

snap34_pr_x_list = []
snap34_pr_y_list = []
snap34_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = snap34_raster.index(x,y)
    val = snap34_raster.read(1)[row,col]
    snap34_pr_x_list.append(x)
    snap34_pr_y_list.append(y)
    snap34_pr_val_list.append(val)
    
snap34_norm_post_proc_df = pd.DataFrame(list(zip(snap34_pr_x_list, snap34_pr_y_list, snap34_pr_val_list)), columns =['xcor', 'ycor', 'snap34_val']) 
snap34_raster_to_points_gdf = gpd.GeoDataFrame(snap34_norm_post_proc_df, geometry=gpd.points_from_xy(x=snap34_norm_post_proc_df.xcor, y=snap34_norm_post_proc_df.ycor), crs = crs_utm)
snap34_raster_to_points_gdf['boolean'] = 0
snap34_raster_to_points_gdf['boolean'].loc[snap34_raster_to_points_gdf['snap34_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
snap34_grid_gdf = gpd.sjoin(grid, snap34_raster_to_points_gdf, how='inner', op='contains')
snap34_grid_norm_gdf = snap34_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#









#*** Airport ***#


airport_raster = rasterio.open(Proxy_Folder + "/lu_airport.tif")
print(airport_raster.crs)
print(airport_raster.count)

airport_pr_x_list = []
airport_pr_y_list = []
airport_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = airport_raster.index(x,y)
    val = airport_raster.read(1)[row,col]
    airport_pr_x_list.append(x)
    airport_pr_y_list.append(y)
    airport_pr_val_list.append(val)
    
airport_norm_post_proc_df = pd.DataFrame(list(zip(airport_pr_x_list, airport_pr_y_list, airport_pr_val_list)), columns =['xcor', 'ycor', 'airport_val']) 
airport_raster_to_points_gdf = gpd.GeoDataFrame(airport_norm_post_proc_df, geometry=gpd.points_from_xy(x=airport_norm_post_proc_df.xcor, y=airport_norm_post_proc_df.ycor), crs = crs_utm)
airport_raster_to_points_gdf['boolean'] = 0
airport_raster_to_points_gdf['boolean'].loc[airport_raster_to_points_gdf['airport_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
airport_grid_gdf = gpd.sjoin(grid, airport_raster_to_points_gdf, how='inner', op='contains')
airport_grid_norm_gdf = airport_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#









#*** Ports ***#


##ports_raster = rasterio.open(Proxy_Folder + "/lu_ports.tif")
ports_raster = rasterio.open(Proxy_Folder + "/lu_shipping_with ports.tif")
print(ports_raster.crs)
print(ports_raster.count)

ports_pr_x_list = []
ports_pr_y_list = []
ports_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = ports_raster.index(x,y)
    val = ports_raster.read(1)[row,col]
    ports_pr_x_list.append(x)
    ports_pr_y_list.append(y)
    ports_pr_val_list.append(val)
    
ports_norm_post_proc_df = pd.DataFrame(list(zip(ports_pr_x_list, ports_pr_y_list, ports_pr_val_list)), columns =['xcor', 'ycor', 'ports_val']) 
ports_raster_to_points_gdf = gpd.GeoDataFrame(ports_norm_post_proc_df, geometry=gpd.points_from_xy(x=ports_norm_post_proc_df.xcor, y=ports_norm_post_proc_df.ycor), crs = crs_utm)
ports_raster_to_points_gdf['boolean'] = 0
ports_raster_to_points_gdf['boolean'].loc[ports_raster_to_points_gdf['ports_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
ports_grid_gdf = gpd.sjoin(grid, ports_raster_to_points_gdf, how='inner', op='contains')
ports_grid_norm_gdf = ports_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#









#*** Snap 1 ***#


snap1_raster = rasterio.open(Proxy_Folder + "/lu_snap1.tif")
print(snap1_raster.crs)
print(snap1_raster.count)

snap1_pr_x_list = []
snap1_pr_y_list = []
snap1_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = snap1_raster.index(x,y)
    val = snap1_raster.read(1)[row,col]
    snap1_pr_x_list.append(x)
    snap1_pr_y_list.append(y)
    snap1_pr_val_list.append(val)
    
snap1_norm_post_proc_df = pd.DataFrame(list(zip(snap1_pr_x_list, snap1_pr_y_list, snap1_pr_val_list)), columns =['xcor', 'ycor', 'snap1_val']) 
snap1_raster_to_points_gdf = gpd.GeoDataFrame(snap1_norm_post_proc_df, geometry=gpd.points_from_xy(x=snap1_norm_post_proc_df.xcor, y=snap1_norm_post_proc_df.ycor), crs = crs_utm)
snap1_raster_to_points_gdf['boolean'] = 0
snap1_raster_to_points_gdf['boolean'].loc[snap1_raster_to_points_gdf['snap1_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
snap1_grid_gdf = gpd.sjoin(grid, snap1_raster_to_points_gdf, how='inner', op='contains')
snap1_grid_norm_gdf = snap1_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean'], axis=1)

#**************************#

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








#**************************************************************************************************************#
#******************************** Proxy (raster) info from grid to CAMS cell ************************************#
#**************************************************************************************************************#








#*** Population Density ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_pop_comb = gpd.sjoin(cams_grid_utm, pop_raster_to_points_gdf, how='inner', op='contains')
cams_grid_pop_comb['cams_index'] = cams_grid_pop_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_pop_comb['center_x'] = cams_grid_pop_comb['geometry'].centroid.x
cams_grid_pop_comb['center_y'] = cams_grid_pop_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_pop_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_pop_group = pd.DataFrame(data=grouped)
cams_grid_pop_group['pop'] = 0
cams_grid_pop_group['pop'].loc[cams_grid_pop_group['boolean'] > 0] = 1

cams_grid_pop = cams_grid_pop_group.drop(['index_right', 'xcor', 'ycor', 'pop_val', 'boolean'], axis=1)
cams_grid_pop['cams_index'] = cams_grid_pop.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_pop = gpd.GeoDataFrame(cams_grid_pop, geometry=gpd.points_from_xy(x=cams_grid_pop.center_x, y=cams_grid_pop.center_y), crs = crs_utm)

#**************************#








#*** Industry ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_ind_comb = gpd.sjoin(cams_grid_utm, ind_raster_to_points_gdf, how='inner', op='contains')
cams_grid_ind_comb['cams_index'] = cams_grid_ind_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_ind_comb['center_x'] = cams_grid_ind_comb['geometry'].centroid.x
cams_grid_ind_comb['center_y'] = cams_grid_ind_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_ind_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_ind_group = pd.DataFrame(data=grouped)
cams_grid_ind_group['ind'] = 0
cams_grid_ind_group['ind'].loc[cams_grid_ind_group['boolean'] > 0] = 1

cams_grid_ind = cams_grid_ind_group.drop(['index_right', 'xcor', 'ycor', 'ind_val', 'boolean'], axis=1)
cams_grid_ind['cams_index'] = cams_grid_ind.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_ind = gpd.GeoDataFrame(cams_grid_ind, geometry=gpd.points_from_xy(x=cams_grid_ind.center_x, y=cams_grid_ind.center_y), crs = crs_utm)

#**************************#








#*** Agriculture ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_agr_comb = gpd.sjoin(cams_grid_utm, agr_raster_to_points_gdf, how='inner', op='contains')
cams_grid_agr_comb['cams_index'] = cams_grid_agr_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_agr_comb['center_x'] = cams_grid_agr_comb['geometry'].centroid.x
cams_grid_agr_comb['center_y'] = cams_grid_agr_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_agr_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_agr_group = pd.DataFrame(data=grouped)
cams_grid_agr_group['agr'] = 0
cams_grid_agr_group['agr'].loc[cams_grid_agr_group['boolean'] > 0] = 1

cams_grid_agr = cams_grid_agr_group.drop(['index_right', 'xcor', 'ycor', 'agr_val', 'boolean'], axis=1)
cams_grid_agr['cams_index'] = cams_grid_agr.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_agr = gpd.GeoDataFrame(cams_grid_agr, geometry=gpd.points_from_xy(x=cams_grid_agr.center_x, y=cams_grid_agr.center_y), crs = crs_utm)

#**************************#








#*** Non Road Mobility Sources ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_snap8_comb = gpd.sjoin(cams_grid_utm, snap8_raster_to_points_gdf, how='inner', op='contains')
cams_grid_snap8_comb['cams_index'] = cams_grid_snap8_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_snap8_comb['center_x'] = cams_grid_snap8_comb['geometry'].centroid.x
cams_grid_snap8_comb['center_y'] = cams_grid_snap8_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_snap8_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_snap8_group = pd.DataFrame(data=grouped)
cams_grid_snap8_group['offroad'] = 0
cams_grid_snap8_group['offroad'].loc[cams_grid_snap8_group['boolean'] > 0] = 1

cams_grid_snap8 = cams_grid_snap8_group.drop(['index_right', 'xcor', 'ycor', 'offroad_val', 'boolean'], axis=1)
cams_grid_snap8['cams_index'] = cams_grid_snap8.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_snap8 = gpd.GeoDataFrame(cams_grid_snap8, geometry=gpd.points_from_xy(x=cams_grid_snap8.center_x, y=cams_grid_snap8.center_y), crs = crs_utm)

#**************************#








#*** Waste ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_waste_comb = gpd.sjoin(cams_grid_utm, waste_raster_to_points_gdf, how='inner', op='contains')
cams_grid_waste_comb['cams_index'] = cams_grid_waste_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_waste_comb['center_x'] = cams_grid_waste_comb['geometry'].centroid.x
cams_grid_waste_comb['center_y'] = cams_grid_waste_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_waste_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_waste_group = pd.DataFrame(data=grouped)
cams_grid_waste_group['waste'] = 0
cams_grid_waste_group['waste'].loc[cams_grid_waste_group['boolean'] > 0] = 1

cams_grid_waste = cams_grid_waste_group.drop(['index_right', 'xcor', 'ycor', 'waste_val', 'boolean'], axis=1)
cams_grid_waste['cams_index'] = cams_grid_waste.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_waste = gpd.GeoDataFrame(cams_grid_waste, geometry=gpd.points_from_xy(x=cams_grid_waste.center_x, y=cams_grid_waste.center_y), crs = crs_utm)

#**************************#








#*** Snap 34 ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_snap34_comb = gpd.sjoin(cams_grid_utm, snap34_raster_to_points_gdf, how='inner', op='contains')
cams_grid_snap34_comb['cams_index'] = cams_grid_snap34_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_snap34_comb['center_x'] = cams_grid_snap34_comb['geometry'].centroid.x
cams_grid_snap34_comb['center_y'] = cams_grid_snap34_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_snap34_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_snap34_group = pd.DataFrame(data=grouped)
cams_grid_snap34_group['snap34'] = 0
cams_grid_snap34_group['snap34'].loc[cams_grid_snap34_group['boolean'] > 0] = 1

cams_grid_snap34 = cams_grid_snap34_group.drop(['index_right', 'xcor', 'ycor', 'snap34_val', 'boolean'], axis=1)
cams_grid_snap34['cams_index'] = cams_grid_snap34.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_snap34 = gpd.GeoDataFrame(cams_grid_snap34, geometry=gpd.points_from_xy(x=cams_grid_snap34.center_x, y=cams_grid_snap34.center_y), crs = crs_utm)

#**************************#








#*** Airport ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_airport_comb = gpd.sjoin(cams_grid_utm, airport_raster_to_points_gdf, how='inner', op='contains')
cams_grid_airport_comb['cams_index'] = cams_grid_airport_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_airport_comb['center_x'] = cams_grid_airport_comb['geometry'].centroid.x
cams_grid_airport_comb['center_y'] = cams_grid_airport_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_airport_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_airport_group = pd.DataFrame(data=grouped)
cams_grid_airport_group['airport'] = 0
cams_grid_airport_group['airport'].loc[cams_grid_airport_group['boolean'] > 0] = 1

cams_grid_airport = cams_grid_airport_group.drop(['index_right', 'xcor', 'ycor', 'airport_val', 'boolean'], axis=1)
cams_grid_airport['cams_index'] = cams_grid_airport.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_airport = gpd.GeoDataFrame(cams_grid_airport, geometry=gpd.points_from_xy(x=cams_grid_airport.center_x, y=cams_grid_airport.center_y), crs = crs_utm)

#**************************#








#*** Ports ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_ports_comb = gpd.sjoin(cams_grid_utm, ports_raster_to_points_gdf, how='inner', op='contains')
cams_grid_ports_comb['cams_index'] = cams_grid_ports_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_ports_comb['center_x'] = cams_grid_ports_comb['geometry'].centroid.x
cams_grid_ports_comb['center_y'] = cams_grid_ports_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_ports_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_ports_group = pd.DataFrame(data=grouped)
cams_grid_ports_group['ports'] = 0
cams_grid_ports_group['ports'].loc[cams_grid_ports_group['boolean'] > 0] = 1

cams_grid_ports = cams_grid_ports_group.drop(['index_right', 'xcor', 'ycor', 'ports_val', 'boolean'], axis=1)
cams_grid_ports['cams_index'] = cams_grid_ports.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_ports = gpd.GeoDataFrame(cams_grid_ports, geometry=gpd.points_from_xy(x=cams_grid_ports.center_x, y=cams_grid_ports.center_y), crs = crs_utm)

#**************************#








#*** Snap 1 ***#


#*** Spatial join CAMS with proxy grid  ***#

cams_grid_snap1_comb = gpd.sjoin(cams_grid_utm, snap1_raster_to_points_gdf, how='inner', op='contains')
cams_grid_snap1_comb['cams_index'] = cams_grid_snap1_comb.index


#*** Calculate centroid OR bounds of grid polygons - CAMS cells with multiple info (rows:2025)  ***#

cams_grid_snap1_comb['center_x'] = cams_grid_snap1_comb['geometry'].centroid.x
cams_grid_snap1_comb['center_y'] = cams_grid_snap1_comb['geometry'].centroid.y


#*** Group by common key (column values) - CAMS cells proxy exist - Boolean (rows:35)  ***#

grouped = cams_grid_snap1_comb.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()
cams_grid_snap1_group = pd.DataFrame(data=grouped)
cams_grid_snap1_group['snap1'] = 0
cams_grid_snap1_group['snap1'].loc[cams_grid_snap1_group['boolean'] > 0] = 1

cams_grid_snap1 = cams_grid_snap1_group.drop(['index_right', 'xcor', 'ycor', 'snap1_val', 'boolean'], axis=1)
cams_grid_snap1['cams_index'] = cams_grid_snap1.index


#*** Proxy Existance in CAMS POINTS  ***#

cams_grid_temp_snap1 = gpd.GeoDataFrame(cams_grid_snap1, geometry=gpd.points_from_xy(x=cams_grid_snap1.center_x, y=cams_grid_snap1.center_y), crs = crs_utm)

#**************************#

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








#**************************************************************************************************************
#********************************************CAMS Final Polygons******************************************************************
#**************************************************************************************************************



cams_grid_temp = cams_grid_temp_pop.append([cams_grid_temp_ind, cams_grid_temp_agr, cams_grid_temp_snap8, cams_grid_temp_waste, cams_grid_temp_snap34, cams_grid_temp_airport, cams_grid_temp_ports, cams_grid_temp_snap1])

cams_group = cams_grid_temp.groupby(['center_x', 'center_y','cams_index'], as_index=False).sum()

cams_group_points = gpd.GeoDataFrame(cams_group, geometry=gpd.points_from_xy(x=cams_group.center_x, y=cams_group.center_y), crs = crs_utm)

cams_grid_utm_sjoin = gpd.sjoin(cams_grid_utm, cams_group_points, how='left', op='contains')

cams_proxy_existance = cams_grid_utm_sjoin.drop(['center_x', 'center_y', 'index_right'], axis=1)

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








#**************************************************************************************************************
#********************************* Country regions for coastline clip *****************************************
#**************************************************************************************************************



if country_code == 'GRC':
    greece_regions_shp = gpd.read_file("X:/NOA/covid/maps/Reviewers/Regions_Greece_UTM34.shp")
    greece_regions_clip = cl.clip_shp(greece_regions_shp, grid)
    greece_regions_clip['dissolve'] = 1
    greece_regions_dis = greece_regions_clip.dissolve(by='dissolve')
    greece_regions_dis.crs = crs_utm

    domain_grid = grid.copy()
    domain_grid['dissolve'] = 1
    domain_dis = domain_grid.dissolve(by='dissolve')
    domain_dis.crs = crs_utm
    greece_coast_union = gpd.overlay(greece_regions_dis, domain_dis, how='union')
    greece_coast_poly = greece_coast_union[['geometry']]
    greece_coast_poly['sea'] = 0
    greece_coast_poly['land'] = 0

    if len(greece_coast_poly.geometry.unique()) > 1:
        greece_coast_poly['sea'].iloc[1] = 1
        greece_coast_poly['land'].iloc[0] = 1
        clip_land_for_sea = greece_coast_poly.loc[greece_coast_poly['sea'] == 1]
        clip_sea_for_land = greece_coast_poly.loc[greece_coast_poly['land'] == 1]
    else:
        greece_coast_poly['land'].iloc[0] = 1
        clip_sea_for_land = greece_coast_poly.loc[greece_coast_poly['land'] == 1]
        clip_land_for_sea = clip_sea_for_land.copy()


else:
    countries_eurostat_2020 = gpd.read_file(Eurostat_countries_shp)


    domain_grid = grid_sea_wgs.copy()
    domain_grid['dissolve'] = 1
    domain_dis_wgs = domain_grid.dissolve(by='dissolve')
    domain_dis_wgs.crs = crs_wgs
    domain_dis = domain_dis_wgs.to_crs(crs_utm)


    countries_eurostat_2020_clip_wgs = cl.clip_shp(countries_eurostat_2020, domain_dis_wgs)
    countries_eurostat_2020_clip_wgs.crs = crs_wgs
    countries_eurostat_2020_clip = countries_eurostat_2020_clip_wgs.to_crs(crs_utm)

    countries_eurostat_2020_clip['dissolve'] = 1
    countries_eurostat_2020_dis = countries_eurostat_2020_clip.dissolve(by='dissolve')
    countries_eurostat_2020_dis.crs = crs_utm
    countries_eurostat_2020_dis.to_file(OutFolder + "countries_eurostat_2020_dis.shp", driver="ESRI Shapefile")
    countries_eurostat_2020_union = gpd.overlay(countries_eurostat_2020_dis, domain_dis, how='union')
    countries_eurostat_2020_poly_temp = countries_eurostat_2020_union[['geometry']]
    countries_eurostat_2020_poly = cl.clip_shp(countries_eurostat_2020_poly_temp, grid)
    countries_eurostat_2020_poly.crs = crs_utm
    countries_eurostat_2020_poly['sea'] = 0
    countries_eurostat_2020_poly['land'] = 0

    if len(countries_eurostat_2020_poly.geometry.unique()) > 1:
        countries_eurostat_2020_poly['sea'].iloc[1] = 1
        countries_eurostat_2020_poly['land'].iloc[0] = 1
        clip_land_for_sea = countries_eurostat_2020_poly.loc[countries_eurostat_2020_poly['sea'] == 1]
        clip_sea_for_land = countries_eurostat_2020_poly.loc[countries_eurostat_2020_poly['land'] == 1]
    else:
        countries_eurostat_2020_poly['land'].iloc[0] = 1
        clip_sea_for_land = countries_eurostat_2020_poly.loc[countries_eurostat_2020_poly['land'] == 1]
        clip_land_for_sea = clip_sea_for_land.copy()
    

###**************************************************************************************************************
###**************************************************************************************************************
###**************************************************************************************************************









#**************************************************************************************************************
#********************************* CAMS points ******************************************************************
#**************************************************************************************************************



try:
    cams_point = gpd.read_file(OutFolder + "cams_points.shp")
    cams_points_to_grid_utm = gpd.sjoin(grid, cams_point, how='inner', op='contains')
    cams_points_to_grid_utm['SNAP'] = cams_points_to_grid_utm['snap']
    cams_points_to_grid_utm['SO2'] = cams_points_to_grid_utm['SOX']
    cams_points_to_grid_utm['PM2_5'] = cams_points_to_grid_utm['PM25']
    cams_points_to_grid_utm.to_file(OutFolder + "cams_points_to_grid_utm.shp", driver="ESRI Shapefile")
except:
    print("There are no cams points from E-PRTR")
    cams_point = pd.DataFrame(columns = ['SNAP', 'xcor', 'ycor', 'Hi', 'Vi', 'Ti', 'radi', 'CH4', 'MainIASect', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2'])
    cams_points_to_grid_utm = gpd.GeoDataFrame(cams_point, crs = crs_utm)

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








#**************************************************************************************************************
#*************************************** Raster to Points Geodataframe ***************************************#
#**************************************************************************************************************



#*** Urban Center ***#


uc_raster = rasterio.open(Proxy_Folder + "/ghs_incr_fact_" + str(uc_increase_factor) + ".tif")
print(uc_raster.crs)
print(uc_raster.count)

uc_pr_x_list = []
uc_pr_y_list = []
uc_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = uc_raster.index(x,y)
    val = uc_raster.read(1)[row,col]
    uc_pr_x_list.append(x)
    uc_pr_y_list.append(y)
    uc_pr_val_list.append(val)
    
uc_norm_post_proc_df = pd.DataFrame(list(zip(uc_pr_x_list, uc_pr_y_list, uc_pr_val_list)), columns =['xcor', 'ycor', 'uc_val']) 
uc_raster_to_points_gdf = gpd.GeoDataFrame(uc_norm_post_proc_df, geometry=gpd.points_from_xy(x=uc_norm_post_proc_df.xcor, y=uc_norm_post_proc_df.ycor), crs = crs_utm)
uc_raster_to_points_gdf['boolean'] = 0
uc_raster_to_points_gdf['boolean'].loc[uc_raster_to_points_gdf['uc_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
uc_grid_gdf = gpd.sjoin(grid, uc_raster_to_points_gdf, how='inner', op='contains')

#*** Keep only UC values ***#
uc_clip = uc_grid_gdf.loc[uc_grid_gdf['boolean'] == 1]
uc_clip['dissolve'] = 1
uc_4_clip = uc_clip.dissolve(by='dissolve')
uc_4_clip.crs = crs_utm
uc_4_clip.to_file(OutFolder + "uc_4_clip.shp", driver="ESRI Shapefile")

#**************************#

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








proxy_grid_list = []
cams_fine_list = []
fine_temp_list = []
small_cell_emission_list = []
cams_no_proxy_fine_list = []
small_cell_no_proxy_emission_list = []
urbem_merged_list = []
urbem_full_list = []
small_cell_no_proxy_redistr_list = []
urbem_merged_redistr_list = []
area_list = []
smooth_list = []
urbem_stat_list = []
cams_stat_list = []







for i in range(0,len(SNAP_sectors)):
    print("***************************** i ********************************")
    print("i", i)
    print("SNAP_sectors", SNAP_sectors[i])
    
    
    #*** Read CAMS emissions croped df from dbase data  ***#

    df = col_data[col_data['SNAP'] == SNAP_sectors[i]]

    size = (df.shape[0], df.shape[1])

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(x=df.xcor, y=df.ycor), crs = crs_utm)


    #*** Read CAMS point emissions ***#

    points_gdf = cams_points_to_grid_utm[cams_points_to_grid_utm['SNAP'] == SNAP_sectors[i]]


    #*** Join CAMS emissions gdf with CAMS polygons  ***#

    cams_emissions_gdf = gpd.sjoin(cams_proxy_existance, gdf, how='left', op='contains')

    #*** CAMS emissions polygons set column values ***#

    cams_emissions = cams_emissions_gdf.drop(['index_right', 'xcor', 'ycor'], axis=1)
    cams_emissions['ISO3'] = 'GRC'
    cams_emissions['Year'] = 2015
    cams_emissions['SNAP'] = SNAP_sectors[i]
    cams_emissions['SourceType'] = 'A'
    cams_emissions_temp = cams_emissions.fillna(0)

    #*** Calculate area in CAMS cells before clip with domian ***#

    cams_emissions_temp['area'] = cams_emissions_temp.area

    #*** Clip with the domain and Calculate area in cliped CAMS cells ***#

    cams_emissions_clip = cl.clip_shp(cams_emissions_temp, grid)
    cams_emissions_clip['clip_area'] = cams_emissions_clip.area

    #*** Calculate area difference and coefficient with clip ***#

    cams_emissions_clip['analogia'] = (cams_emissions_clip['clip_area'] - cams_emissions_clip['area'])/cams_emissions_clip['area']
    cams_emissions_clip['area_coef'] = 1 + cams_emissions_clip['analogia']

    #*** Calculate emissions based on area coefficient ***#

    cams_emissions_clip['CH4'] = cams_emissions_clip['CH4'] * cams_emissions_clip['area_coef']
    cams_emissions_clip['CO'] = cams_emissions_clip['CO'] * cams_emissions_clip['area_coef']
    cams_emissions_clip['NH3'] = cams_emissions_clip['NH3'] * cams_emissions_clip['area_coef']
    cams_emissions_clip['NMVOC'] = cams_emissions_clip['NMVOC'] * cams_emissions_clip['area_coef']
    cams_emissions_clip['NOX'] = cams_emissions_clip['NOX'] * cams_emissions_clip['area_coef']
    cams_emissions_clip['PM10'] = cams_emissions_clip['PM10'] * cams_emissions_clip['area_coef']
    cams_emissions_clip['PM2_5'] = cams_emissions_clip['PM2_5'] * cams_emissions_clip['area_coef']
    cams_emissions_clip['SO2'] = cams_emissions_clip['SO2'] * cams_emissions_clip['area_coef']

    
    CAMS_emissions = cams_emissions_clip.drop(['area', 'clip_area', 'analogia', 'area_coef'], axis=1)







###*******************************************************************#
###*******************************************************************#
###****************** Distribution per Snap sector *******************#
###*******************************************************************#
###*******************************************************************#








#*** Snap 1 ***#



    if (SNAP_sectors[i] == 1):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['snap1'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['snap1'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(snap1_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')

            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0



            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()


            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            snap1_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            snap1_cells_with_proxy = cl.clip_shp(snap1_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(snap1_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = snap1_cells_with_proxy.loc[snap1_cells_with_proxy['cams_index'] == cams_cells[big_cell]]

                ###*** Normalize snap1 val per CAMS cell ***#
                
                small_cell_per_cams_index_df['snap1_val_norm'] = 0
                small_cell_per_cams_index_df['snap1_val'].loc[small_cell_per_cams_index_df['snap1_val'] < 0] = 0
                small_cell_per_cams_index_df['snap1_val_norm'] = small_cell_per_cams_index_df['snap1_val']/sum(small_cell_per_cams_index_df['snap1_val'])


                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['snap1_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['snap1_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['snap1_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['snap1_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['snap1_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['snap1_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['snap1_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['snap1_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

            

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(snap1_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                
                ##print(cams_fine_list)

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            snap1_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            snap1_cells_with_no_proxy = cl.clip_shp(snap1_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(snap1_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = snap1_cells_with_no_proxy.loc[snap1_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['snap1_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)


                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['snap1'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['snap1'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)

                ##print(sum(small_cell_redistr['CH4']))

                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point) > 0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'snap1_val', 'snap1_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 2 / Snap 6 / Snap 7 ***#



    elif ((SNAP_sectors[i] == 2) | (SNAP_sectors[i] == 6) | (SNAP_sectors[i] == 7)):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['pop'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['pop'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(pop_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            ##print(proxy_cams_emissions_gdf)
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0

            ##proxy_cams_emissions_gdf.to_file(OutFolder + "proxy_cams_emissions_gdf_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            pop_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            pop_cells_with_proxy = cl.clip_shp(pop_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(pop_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = pop_cells_with_proxy.loc[pop_cells_with_proxy['cams_index'] == cams_cells[big_cell]]

                ###*** Normalize pop val per CAMS cell ***#
                
                small_cell_per_cams_index_df['pop_val_norm'] = 0
                small_cell_per_cams_index_df['pop_val'].loc[small_cell_per_cams_index_df['pop_val'] < 0] = 0
                small_cell_per_cams_index_df['pop_val_norm'] = small_cell_per_cams_index_df['pop_val']/sum(small_cell_per_cams_index_df['pop_val'])


                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['pop_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['pop_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['pop_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['pop_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['pop_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['pop_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['pop_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['pop_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

            

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(pop_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]
 
            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            pop_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            pop_cells_with_no_proxy = cl.clip_shp(pop_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(pop_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = pop_cells_with_no_proxy.loc[pop_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['pop_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)

                ##print(sum(small_cell_per_cams_index_df['CH4']))

                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)


        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['pop'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['pop'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)


                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])










        ### ~ Urban Center increase and save to shp



        if SNAP_sectors[i] == 7:
            ghs = gdal.Open(Proxy_Folder + "/ghs_incr_fact_" + str(uc_increase_factor) + ".tif")
            ghs_table = ghs.GetRasterBand(1)
            ghs_geo_transform = ghs.GetGeoTransform()
            ghs_np_table = ghs_table.ReadAsArray()

            ghs_np_table[ghs_np_table == 0]= 1
            ghs_np_fl = ghs_np_table.flatten()
            ghs_norm_df = pd.DataFrame(data=ghs_np_fl, columns=["boolean"])

            cams_ghs_join = pd.concat([urbem_final_shp, ghs_norm_df], axis=1, join='inner')


            cams_ghs_all_increase = cams_ghs_join.copy()
            cams_ghs_all_increase.loc[:,'CH4_km'] *= cams_ghs_all_increase.loc[:,'boolean']
            cams_ghs_all_increase.loc[:,'NOX_km'] *= cams_ghs_all_increase.loc[:,'boolean']
            cams_ghs_all_increase.loc[:,'NMVOC_km'] *= cams_ghs_all_increase.loc[:,'boolean']
            cams_ghs_all_increase.loc[:,'CO_km'] *= cams_ghs_all_increase.loc[:,'boolean']
            cams_ghs_all_increase.loc[:,'SO2_km'] *= cams_ghs_all_increase.loc[:,'boolean']
            cams_ghs_all_increase.loc[:,'NH3_km'] *= cams_ghs_all_increase.loc[:,'boolean']
            cams_ghs_all_increase.loc[:,'PM2_5_km'] *= cams_ghs_all_increase.loc[:,'boolean']
            cams_ghs_all_increase.loc[:,'PM10_km'] *= cams_ghs_all_increase.loc[:,'boolean']

            cams_ghs_all_increase_final = cams_ghs_all_increase#.drop(['boolean'], axis=1)
            print(cams_ghs_all_increase_final)
            urbem_final_uc_increase_shp = cams_ghs_all_increase_final.copy()
            urbem_final_uc_increase_shp.to_file(OutFolder + "urbem_final_uc_increase_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")









###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")
        ##print(CAMS_emissions)
        print(CAMS_emissions_shp['CH4'].sum())

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'pop_val', 'pop_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 5 ***#



    elif (SNAP_sectors[i] == 5):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['ind'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['ind'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(ind_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            ##proxy_cams_emissions_gdf = gpd.sjoin(grid_proxy_cliped, cams_proxy, how='inner', op='intersects')
            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0


            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):
            ##for j in range (0, 1):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                ##print(grid_index_dupl_cells)

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            ind_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            ind_cells_with_proxy = cl.clip_shp(ind_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(ind_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = ind_cells_with_proxy.loc[ind_cells_with_proxy['cams_index'] == cams_cells[big_cell]]

                ###*** Normalize ind val per CAMS cell ***#
                
                small_cell_per_cams_index_df['ind_val_norm'] = 0
                small_cell_per_cams_index_df['ind_val'].loc[small_cell_per_cams_index_df['ind_val'] < 0] = 0
                small_cell_per_cams_index_df['ind_val_norm'] = small_cell_per_cams_index_df['ind_val']/sum(small_cell_per_cams_index_df['ind_val'])

                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['ind_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['ind_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['ind_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['ind_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['ind_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['ind_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['ind_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['ind_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(ind_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')
            ##grid_no_proxy_full_emissions.fillna(0)


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            ind_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            ind_cells_with_no_proxy = cl.clip_shp(ind_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(ind_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = ind_cells_with_no_proxy.loc[ind_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['ind_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)

                ##print(sum(small_cell_per_cams_index_df['CH4']))

                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                    ##print(grid_index_dupl_cells)

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    ##print(grid_index_dupl_cells)
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['ind'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['ind'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):
            ##for big_cell in range(0, 1):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)

                ##print(sum(small_cell_redistr['CH4']))

                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")
        ##print(CAMS_emissions)
        print(CAMS_emissions_shp['CH4'].sum())

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'ind_val', 'ind_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 8 ***#



    elif (SNAP_sectors[i] == 8):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_land_for_sea)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['ports'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()
        ##print("**************************cams_proxy**************************")
        ##print(cams_proxy)

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['ports'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()
        ##print("**************************cams_no_proxy**************************")
        ##print(cams_no_proxy)








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(ports_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]
            ##grid_proxy_cliped.to_file(OutFolder + "grid_proxy_cliped_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            ##proxy_cams_emissions_gdf = gpd.sjoin(grid_proxy_cliped, cams_proxy, how='inner', op='intersects')
            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            ##print(proxy_cams_emissions_gdf)
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0

            ##proxy_cams_emissions_gdf.to_file(OutFolder + "proxy_cams_emissions_gdf_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()
            ##print(cams_index_cells_with_proxy)

            
            ##proxy_grid_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == 1324]
            ##print(proxy_grid_cells)

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):
            ##for j in range (0, 1):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                ##print(grid_index_dupl_cells)

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    ##print(pollutants_check)
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]
                    ##print(pollutants_with_value)

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        ##print(grid_index_dupl_cells)
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        ##print(grid_index_dupl_cells)
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                
                ##print(cams_fine_list)

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		
            ##print(grid_proxy_with_gaps_overlay)
            ##grid_proxy_with_gaps_overlay.to_file(OutFolder + "grid_proxy_with_gaps_overlay_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            ##grid_proxy_full_emissions.to_file(OutFolder + "grid_proxy_full_emissions_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")
            ##print(grid_proxy_full_emissions)
            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            ports_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            ports_cells_with_proxy = cl.clip_shp(ports_cells_with_proxy_temp, cams_proxy)
            ##ports_cells_with_proxy.to_file(OutFolder + "ports_cells_with_proxy_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(ports_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):
            ##for big_cell in range(23, 24):

                small_cell_per_cams_index_df = ports_cells_with_proxy.loc[ports_cells_with_proxy['cams_index'] == cams_cells[big_cell]]

                ###*** Normalize ports val per CAMS cell ***#
                
                small_cell_per_cams_index_df['ports_val_norm'] = 0
                small_cell_per_cams_index_df['ports_val'].loc[small_cell_per_cams_index_df['ports_val'] < 0] = 0
                small_cell_per_cams_index_df['ports_val_norm'] = small_cell_per_cams_index_df['ports_val']/sum(small_cell_per_cams_index_df['ports_val'])

                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['ports_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['ports_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['ports_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['ports_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['ports_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['ports_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['ports_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['ports_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(ports_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                
                ##print(cams_fine_list)

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            ports_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            ports_cells_with_no_proxy = cl.clip_shp(ports_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(ports_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = ports_cells_with_no_proxy.loc[ports_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['ports_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)


                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['ports'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['ports'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)

                ##print(sum(small_cell_redistr['CH4']))

                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_land_for_sea)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")
        print(CAMS_emissions_shp['CH4'].sum())

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'ports_val', 'ports_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 9 ***#



    elif (SNAP_sectors[i] == 9):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['waste'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['waste'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(waste_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0


            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            waste_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            waste_cells_with_proxy = cl.clip_shp(waste_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(waste_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = waste_cells_with_proxy.loc[waste_cells_with_proxy['cams_index'] == cams_cells[big_cell]]

                ###*** Normalize waste val per CAMS cell ***#
                
                small_cell_per_cams_index_df['waste_val_norm'] = 0
                small_cell_per_cams_index_df['waste_val'].loc[small_cell_per_cams_index_df['waste_val'] < 0] = 0
                small_cell_per_cams_index_df['waste_val_norm'] = small_cell_per_cams_index_df['waste_val']/sum(small_cell_per_cams_index_df['waste_val'])

                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['waste_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['waste_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['waste_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['waste_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['waste_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['waste_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['waste_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['waste_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(waste_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            waste_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            waste_cells_with_no_proxy = cl.clip_shp(waste_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(waste_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = waste_cells_with_no_proxy.loc[waste_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['waste_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)


                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['waste'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['waste'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)


                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)


        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'waste_val', 'waste_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 10 ***#



    elif (SNAP_sectors[i] == 10):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['agr'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['agr'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(agr_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0


            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            agr_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            agr_cells_with_proxy = cl.clip_shp(agr_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(agr_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = agr_cells_with_proxy.loc[agr_cells_with_proxy['cams_index'] == cams_cells[big_cell]]

                ###*** Normalize agr val per CAMS cell ***#
                
                small_cell_per_cams_index_df['agr_val_norm'] = 0
                small_cell_per_cams_index_df['agr_val'].loc[small_cell_per_cams_index_df['agr_val'] < 0] = 0
                small_cell_per_cams_index_df['agr_val_norm'] = small_cell_per_cams_index_df['agr_val']/sum(small_cell_per_cams_index_df['agr_val'])

                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['agr_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['agr_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['agr_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['agr_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['agr_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['agr_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['agr_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['agr_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(agr_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            agr_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            agr_cells_with_no_proxy = cl.clip_shp(agr_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(agr_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = agr_cells_with_no_proxy.loc[agr_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['agr_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)


                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                    ##print(grid_index_dupl_cells)

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['agr'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['agr'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)

                ##print(sum(small_cell_redistr['CH4']))

                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'agr_val', 'agr_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 11 ***#



    elif (SNAP_sectors[i] == 11):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['airport'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['airport'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(airport_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0


            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()


            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                
                ##print(cams_fine_list)

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            airport_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            airport_cells_with_proxy = cl.clip_shp(airport_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(airport_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = airport_cells_with_proxy.loc[airport_cells_with_proxy['cams_index'] == cams_cells[big_cell]]

                ###*** Normalize airport val per CAMS cell ***#
                
                small_cell_per_cams_index_df['airport_val_norm'] = 0
                small_cell_per_cams_index_df['airport_val'].loc[small_cell_per_cams_index_df['airport_val'] < 0] = 0
                small_cell_per_cams_index_df['airport_val_norm'] = small_cell_per_cams_index_df['airport_val']/sum(small_cell_per_cams_index_df['airport_val'])

                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['airport_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['airport_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['airport_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['airport_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['airport_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['airport_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['airport_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['airport_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(airport_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            airport_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            airport_cells_with_no_proxy = cl.clip_shp(airport_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(airport_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = airport_cells_with_no_proxy.loc[airport_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['airport_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)

                ##print(sum(small_cell_per_cams_index_df['CH4']))

                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['airport'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['airport'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)


                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)


        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'airport_val', 'airport_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 12 ***#



    elif (SNAP_sectors[i] == 12):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['offroad'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['offroad'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(snap8_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0


            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            snap8_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            snap8_cells_with_proxy = cl.clip_shp(snap8_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(snap8_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = snap8_cells_with_proxy.loc[snap8_cells_with_proxy['cams_index'] == cams_cells[big_cell]]
                
                small_cell_per_cams_index_df['offroad_val_norm'] = 0
                small_cell_per_cams_index_df['offroad_val'].loc[small_cell_per_cams_index_df['offroad_val'] < 0] = 0
                small_cell_per_cams_index_df['offroad_val_norm'] = small_cell_per_cams_index_df['offroad_val']/sum(small_cell_per_cams_index_df['offroad_val'])


                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['offroad_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['offroad_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['offroad_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['offroad_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['offroad_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['offroad_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['offroad_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['offroad_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(snap8_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]
                    ##print(pollutants_with_value)

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            snap8_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            snap8_cells_with_no_proxy = cl.clip_shp(snap8_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(snap8_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = snap8_cells_with_no_proxy.loc[snap8_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['offroad_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)

                ##print(sum(small_cell_per_cams_index_df['CH4']))

                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['offroad'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['offroad'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)

                ##print(sum(small_cell_redistr['CH4']))

                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'offroad_val', 'offroad_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###*** Snap 3 / Snap 4  / Snap 34 ***#



    elif ((SNAP_sectors[i] == 3) | (SNAP_sectors[i] == 4) | (SNAP_sectors[i] == 34)):

        #*** Clip CAMS polygons with sea (KEEP only LAND polygons) ***#

        if Domain_with_sea == 'YES':
            CAMS_emissions = cl.clip_shp(CAMS_emissions, clip_sea_for_land)
            CAMS_emissions = CAMS_emissions[~CAMS_emissions.is_empty]
        else:
            CAMS_emissions = CAMS_emissions.copy()

        
        #*** Proxy existance in CAMS cells - Separate CAMS emissions polygons ***#

        cams_proxy = CAMS_emissions.loc[CAMS_emissions['snap34'] == 1]
        cams_proxy_index = cams_proxy.cams_index.unique()

        cams_no_proxy = CAMS_emissions.loc[CAMS_emissions['snap34'] == 0]
        cams_no_proxy_index = cams_no_proxy.cams_index.unique()








###*******************************************************************#
###*************** Check if proxy EXISTS per grid cell ***************#
###*******************************************************************#








        if len(cams_proxy) > 0:

            #*** Clip Proxy polygons (contain proxy value) ***#

            grid_proxy_cliped = cl.clip_shp(snap34_grid_norm_gdf, cams_proxy)
            grid_proxy_cliped = grid_proxy_cliped[~grid_proxy_cliped.is_empty]

            #*** Spatial join grid (with proxy val) with CAMS cell emission values ***#

            proxy_cams_emissions_gdf = gpd.overlay(cams_proxy, grid_proxy_cliped, how='intersection')
            proxy_cams_emissions_gdf['area'] = proxy_cams_emissions_gdf.area
            proxy_cams_emissions_gdf['max_area'] = 0


            cams_index_cells_with_proxy = proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
            grid_index = proxy_cams_emissions_gdf['grid_index'].unique()
            proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):



                #*** Select unique num IN grid cell index ***#
                grid_index_dupl_cells = proxy_cams_emissions_gdf.loc[proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]

                grid_index_dupl_cells.fillna(-99)



                #*** Duplicate Cells with multiple AREA values ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_proxy_temp['nn_pixels'].loc[grid_proxy_temp['max_area'] == 1] = 1
                    cams_fine_list.append(grid_proxy_temp)

                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        
                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_fine_list.append(grid_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_proxy_with_gaps_overlay = pd.concat(cams_fine_list)		

            ###*** Whole grid cells with corrected value per cell ***#

            grid_proxy_full_emissions = gpd.sjoin(grid, grid_proxy_with_gaps_overlay, how='inner', op='contains')

            grid_proxy_full_emissions['grid_index'] = grid_proxy_full_emissions['grid_index_left']
            snap34_cells_with_proxy_temp = grid_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            snap34_cells_with_proxy = cl.clip_shp(snap34_cells_with_proxy_temp, cams_proxy)

        


            ###*** Analyze grid cells per CAMS cell ***#

            cams_cells = list(snap34_cells_with_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells)):

                small_cell_per_cams_index_df = snap34_cells_with_proxy.loc[snap34_cells_with_proxy['cams_index'] == cams_cells[big_cell]]
                
                small_cell_per_cams_index_df['snap34_val_norm'] = 0
                small_cell_per_cams_index_df['snap34_val'].loc[small_cell_per_cams_index_df['snap34_val'] < 0] = 0
                small_cell_per_cams_index_df['snap34_val_norm'] = small_cell_per_cams_index_df['snap34_val']/sum(small_cell_per_cams_index_df['snap34_val'])

                ###*** Calculate EMISSIONS per small cell ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())* small_cell_per_cams_index_df['snap34_val_norm']
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())* small_cell_per_cams_index_df['snap34_val_norm']
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())* small_cell_per_cams_index_df['snap34_val_norm']
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())* small_cell_per_cams_index_df['snap34_val_norm']
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())* small_cell_per_cams_index_df['snap34_val_norm']
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())* small_cell_per_cams_index_df['snap34_val_norm']
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())* small_cell_per_cams_index_df['snap34_val_norm']
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())* small_cell_per_cams_index_df['snap34_val_norm']


                small_cell_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_proxy = pd.concat(small_cell_emission_list)

        else:
            print("There are no proxy values...")

            urbem_with_proxy = cams_proxy.copy()








###****************************************************************************#
###*************** Check if proxy DOESN' T EXISTS per grid cell ***************#
###****************************************************************************#








        if len(cams_no_proxy) > 0:

            #*** Clip no_proxy polygons (contain no_proxy value) ***#

            grid_no_proxy_cliped = cl.clip_shp(snap34_grid_norm_gdf, cams_no_proxy)
            grid_no_proxy_cliped = grid_no_proxy_cliped[~grid_no_proxy_cliped.is_empty]

            #*** Spatial join grid (with no_proxy val) with CAMS cell emission values ***#

            no_proxy_cams_emissions_gdf = gpd.overlay(cams_no_proxy, grid_no_proxy_cliped, how='intersection')
            no_proxy_cams_emissions_gdf['area'] = no_proxy_cams_emissions_gdf.area
            no_proxy_cams_emissions_gdf['max_area'] = 0

            cams_index_cells_with_no_proxy = no_proxy_cams_emissions_gdf.cams_index.unique()

            

            #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#

            grid_index = no_proxy_cams_emissions_gdf['grid_index'].unique()
            no_proxy_cams_emissions_gdf['nn_pixels'] = 0





            for j in range (0, len(grid_index)):

                grid_index_dupl_cells = no_proxy_cams_emissions_gdf.loc[no_proxy_cams_emissions_gdf['grid_index'] == grid_index[j]]
                grid_index_dupl_cells.fillna(-99)


                #*** Cell with GREATER AREA value ***#
                if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                    grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    grid_no_proxy_temp['nn_pixels'].loc[grid_no_proxy_temp['max_area'] == 1] = 1
                    cams_no_proxy_fine_list.append(grid_no_proxy_temp)



                #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
                elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):
                    print(grid_index_dupl_cells)

                    #*** Select only EMISSION cols ***#
                    pollutants_check = grid_index_dupl_cells.iloc[:, 11:19]
                    


                    #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                    pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                    if len(pollutants_with_value.columns) > 0:

                        #*** From the remaining EMISSION cols select the first ***#
                        first_col_with_val = pollutants_with_value.iloc[:, :1]

                        #*** From this col select COL NAME ***#
                        col_name = list(first_col_with_val)[0]

                        #*** From this col select MAX ROW value ***#
                        max_pollutant_val = float(first_col_with_val.max())

                        #*** From duplicate cells select ROW containing this MAX value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1

                        #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)


                    else:

                        #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                        grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                        grid_index_dupl_cells['nn_pixels'].loc[grid_index_dupl_cells['max_area'] == 1] = 1
                        grid_no_proxy_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                        cams_no_proxy_fine_list.append(grid_no_proxy_temp)
                        
                else:
                    grid_index_dupl_cells['max_area'] = 1
                    cams_no_proxy_fine_list.append(grid_index_dupl_cells)
                

            #*** Cliped grid cells with corrected value per cell ***#

            grid_no_proxy_with_gaps_overlay = pd.concat(cams_no_proxy_fine_list)		

            #*** Whole grid cells with corrected value per cell ***#

            grid_no_proxy_full_emissions = gpd.sjoin(grid, grid_no_proxy_with_gaps_overlay, how='inner', op='contains')


            grid_no_proxy_full_emissions['grid_index'] = grid_no_proxy_full_emissions['grid_index_left']
            snap34_cells_with_no_proxy_temp = grid_no_proxy_full_emissions.drop(['index_right', 'grid_index_left', 'grid_index_right'], axis=1)
            snap34_cells_with_no_proxy = cl.clip_shp(snap34_cells_with_no_proxy_temp, cams_no_proxy)


            #*** Analyze grid cells per CAMS cell ***#

            cams_cells_no_proxy = list(snap34_cells_with_no_proxy.cams_index.unique())





            for big_cell in range(0, len(cams_cells_no_proxy)):


                small_cell_per_cams_index_df = snap34_cells_with_no_proxy.loc[snap34_cells_with_no_proxy['cams_index'] == cams_cells_no_proxy[big_cell]]

                small_cell_per_cams_index_df['snap34_val_norm'] = 0

                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_per_cams_index_df['CH4_km'] = float(small_cell_per_cams_index_df['CH4'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['CO_km'] = float(small_cell_per_cams_index_df['CO'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NH3_km'] = float(small_cell_per_cams_index_df['NH3'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NMVOC_km'] = float(small_cell_per_cams_index_df['NMVOC'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['NOX_km'] = float(small_cell_per_cams_index_df['NOX'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM10_km'] = float(small_cell_per_cams_index_df['PM10'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['PM2_5_km'] = float(small_cell_per_cams_index_df['PM2_5'].unique())/len(small_cell_per_cams_index_df)
                small_cell_per_cams_index_df['SO2_km'] = float(small_cell_per_cams_index_df['SO2'].unique())/len(small_cell_per_cams_index_df)

                ##print(sum(small_cell_per_cams_index_df['CH4']))

                small_cell_no_proxy_emission_list.append(small_cell_per_cams_index_df)

            urbem_with_no_proxy = pd.concat(small_cell_no_proxy_emission_list)

        else:
            print("There are no no_proxy values...")

            urbem_with_no_proxy = cams_no_proxy.copy()









###****************************************************************************#
###********************* Final Calculations per grid cell *********************#
###****************************************************************************#








        #*** Merge grid cells with and without proxies ***#
        urbem_merged_list.append(urbem_with_no_proxy)
        urbem_merged_list.append(urbem_with_proxy)

        urbem_full_temp = pd.concat(urbem_merged_list).reset_index(drop = True)
        urbem_full = urbem_full_temp.drop(['max_area'], axis=1)
        urbem_full['max_area'] = 0


        #*** Find duplicate index of grid cells and select max areas' value per grid cell ***#
        urbem_grid_index = urbem_full.grid_index.unique()





        for k in range (0, len(urbem_grid_index)):

            grid_index_dupl_cells = urbem_full.loc[urbem_full['grid_index'] == urbem_grid_index[k]]
            grid_index_dupl_cells.fillna(-99)


            #*** Cell with GREATER AREA value ***#
            if (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) != 1):
                grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['area'].idxmax()] = 1
                grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                urbem_full_list.append(grid_temp)



            #*** If AREA values are EUQAL - Select cell with greater EMISSION value ***#
            elif (len(grid_index_dupl_cells) > 1) & (len(grid_index_dupl_cells['area'].unique()) == 1):

                #*** Select only EMISSION cols ***#
                pollutants_check = grid_index_dupl_cells.iloc[:, 12:20]
                


                #*** If whole EMISSION cols are nan or 0 ignore - Keep only EMISSION cols with values ***#
                pollutants_with_value = pollutants_check.loc[:, ((pollutants_check != 0).any(axis=0) & (pollutants_check != -99).any(axis=0))]

                if len(pollutants_with_value.columns) > 0:

                    #*** From the remaining EMISSION cols select the first ***#
                    first_col_with_val = pollutants_with_value.iloc[:, :1]

                    #*** From this col select COL NAME ***#
                    col_name = list(first_col_with_val)[0]

                    #*** From this col select MAX ROW value ***#
                    max_pollutant_val = float(first_col_with_val.max())

                    #*** From duplicate cells select ROW containing this MAX value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells[col_name] == max_pollutant_val] = 1

                    #*** Select only specific duplicate cell ROW with max area value to 1 and append to list***#
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)


                else:

                    #*** If ALL EMISSION cols are nan or 0 - Keep EMISSION ROW with MAX cams index value ***#
                    grid_index_dupl_cells['max_area'].loc[grid_index_dupl_cells['cams_index'].idxmax()] = 1
                    grid_temp = grid_index_dupl_cells.loc[grid_index_dupl_cells['max_area'] == 1]
                    urbem_full_list.append(grid_temp)
                    
            else:
                grid_index_dupl_cells['max_area'] = 1
                urbem_full_list.append(grid_index_dupl_cells)
            

        #*** Cliped grid cells with corrected value per cell ***#

        urbem_full_with_gaps_overlay = pd.concat(urbem_full_list)		

        #*** Whole grid cells with corrected value per cell ***#
        urbem_temp = gpd.sjoin(grid, urbem_full_with_gaps_overlay, how='inner', op='contains')
        urbem_temp['grid_index'] = urbem_temp['grid_index_left']
        urbem_final = urbem_temp.drop(['grid_index_left', 'index_right', 'grid_index_right'], axis=1)












        #*** Redistribute emissions only in small cells without proxies after final grid index selection ***#
        urbem_redistr_proxy = urbem_final.loc[urbem_final['snap34'] == 1]
        urbem_redistr_no_proxy = urbem_final.loc[urbem_final['snap34'] == 0]

        urbem_redistr_cams_index = urbem_redistr_no_proxy.cams_index.unique()



        #*** Avoid missing emissions from no proxy cells ***#
        if len(urbem_redistr_no_proxy) == 0:
            urbem_with_no_proxy_redistr = urbem_redistr_no_proxy.copy()
        elif len(urbem_redistr_no_proxy) != 0:
            
            for big_cell_redistr in range(0, len(urbem_redistr_cams_index)):


                small_cell_redistr = urbem_redistr_no_proxy.loc[urbem_redistr_no_proxy['cams_index'] == urbem_redistr_cams_index[big_cell_redistr]]


                #*** Calculate EMISSIONS per small cell based on COUNT ***#

                small_cell_redistr['CH4_km'] = float(small_cell_redistr['CH4'].unique())/len(small_cell_redistr)
                small_cell_redistr['CO_km'] = float(small_cell_redistr['CO'].unique())/len(small_cell_redistr)
                small_cell_redistr['NH3_km'] = float(small_cell_redistr['NH3'].unique())/len(small_cell_redistr)
                small_cell_redistr['NMVOC_km'] = float(small_cell_redistr['NMVOC'].unique())/len(small_cell_redistr)
                small_cell_redistr['NOX_km'] = float(small_cell_redistr['NOX'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM10_km'] = float(small_cell_redistr['PM10'].unique())/len(small_cell_redistr)
                small_cell_redistr['PM2_5_km'] = float(small_cell_redistr['PM2_5'].unique())/len(small_cell_redistr)
                small_cell_redistr['SO2_km'] = float(small_cell_redistr['SO2'].unique())/len(small_cell_redistr)

                ##print(sum(small_cell_redistr['CH4']))

                small_cell_no_proxy_redistr_list.append(small_cell_redistr)

            urbem_with_no_proxy_redistr = pd.concat(small_cell_no_proxy_redistr_list)

        #*** Merge grid cells with and without proxies ***#
        urbem_merged_redistr_list.append(urbem_with_no_proxy_redistr)
        urbem_merged_redistr_list.append(urbem_redistr_proxy)

        urbem_final = pd.concat(urbem_merged_redistr_list).reset_index(drop = True)
        print(urbem_final['CH4_km'].sum())

        
        #*** Avoid duplication due to point emissions ***#
        area = urbem_final.copy()

        point = points_gdf.replace(-999, 0)

        if len(point)>0:
            pollutants = point[['CH4', 'CO', 'NH3', 'NMVOC', 'NOX', 'PM10', 'PM2_5', 'SO2']]


            point_grid_cell = point.grid_index.unique()

            for p in pollutants:

                if ((area[str(p) + '_km'].sum()) - (point[str(p)].sum()))<=0:
                    for cell in point_grid_cell:
                        cell_value = 0
                        area[p].loc[area['grid_index'] == cell] = cell_value
                else:
                    for cell in point_grid_cell:
                        cell_value = float(area[p].loc[area['grid_index'] == cell])/(area[str(p) + '_km'].sum())*(area[str(p) + '_km'].sum()) - (point[str(p)].sum())
                        area[p].loc[area['grid_index'] == cell] = cell_value
        elif len(point) == 0:
            area = urbem_final.copy()
            

        ### ~ Save to shp as tn/km*year
        urbem_final_shp = cl.clip_shp(area, clip_sea_for_land)
        urbem_final_shp['CH4_km'] = urbem_final_shp['CH4_km']/1000
        urbem_final_shp['NOX_km'] = urbem_final_shp['NOX_km']/1000
        urbem_final_shp['NMVOC_km'] = urbem_final_shp['NMVOC_km']/1000
        urbem_final_shp['CO_km'] = urbem_final_shp['CO_km']/1000
        urbem_final_shp['SO2_km'] = urbem_final_shp['SO2_km']/1000
        urbem_final_shp['NH3_km'] = urbem_final_shp['NH3_km']/1000
        urbem_final_shp['PM2_5_km'] = urbem_final_shp['PM2_5_km']/1000
        urbem_final_shp['PM10_km'] = urbem_final_shp['PM10_km']/1000
        urbem_final_shp.to_file(OutFolder + "urbem_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        urbem_stat_list.append([SNAP_sectors[i], urbem_final_shp['CH4_km'].sum(), urbem_final_shp['NOX_km'].sum(), urbem_final_shp['NMVOC_km'].sum(), urbem_final_shp['CO_km'].sum(), urbem_final_shp['SO2_km'].sum(), urbem_final_shp['NH3_km'].sum(), urbem_final_shp['PM2_5_km'].sum(), urbem_final_shp['PM10_km'].sum(), len(urbem_final_shp)])








###****************************************************************************#
###*************************** CAMS Emissios per km ***************************#
###****************************************************************************#








        #*** Select unique cams index from grid ***#
        urbem_cams_index = urbem_final.cams_index.unique()

        CAMS_emissions['grid_pixel_count'] = 0





        #*** Select PIXEL COUNT per cams index from grid ***#
        for z in range (0, len(urbem_cams_index)):
            pixel_count = len(urbem_final.loc[urbem_final['cams_index'] == urbem_cams_index[z]])
            CAMS_emissions['grid_pixel_count'].loc[CAMS_emissions['cams_index'] == urbem_cams_index[z]] = pixel_count

        #*** Calculate CAMS EMISSIONS per km ***#
        CAMS_emissions['CH4_km'] = CAMS_emissions['CH4']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['CO_km'] = CAMS_emissions['CO']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NH3_km'] = CAMS_emissions['NH3']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NMVOC_km'] = CAMS_emissions['NMVOC']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['NOX_km'] = CAMS_emissions['NOX']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM10_km'] = CAMS_emissions['PM10']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['PM25_km'] = CAMS_emissions['PM2_5']/CAMS_emissions['grid_pixel_count']
        CAMS_emissions['SO2_km'] = CAMS_emissions['SO2']/CAMS_emissions['grid_pixel_count']

        ### ~ Save to shp as tn/km*year
        CAMS_emissions_shp = CAMS_emissions.copy()
        CAMS_emissions_shp['CH4_km'] = CAMS_emissions_shp['CH4_km']/1000
        CAMS_emissions_shp['NOX_km'] = CAMS_emissions_shp['NOX_km']/1000
        CAMS_emissions_shp['NMVOC_km'] = CAMS_emissions_shp['NMVOC_km']/1000
        CAMS_emissions_shp['CO_km'] = CAMS_emissions_shp['CO_km']/1000
        CAMS_emissions_shp['SO2_km'] = CAMS_emissions_shp['SO2_km']/1000
        CAMS_emissions_shp['NH3_km'] = CAMS_emissions_shp['NH3_km']/1000
        CAMS_emissions_shp['PM25_km'] = CAMS_emissions_shp['PM25_km']/1000
        CAMS_emissions_shp['PM10_km'] = CAMS_emissions_shp['PM10_km']/1000
        CAMS_emissions_shp.to_file(OutFolder + "CAMS_emissions_final_snap_" + str(SNAP_sectors[i]) + ".shp", driver="ESRI Shapefile")

        ### ~ Statistics tn/km*year

        cams_stat_list.append([SNAP_sectors[i], (CAMS_emissions_shp['CH4']/1000).sum(), (CAMS_emissions_shp['NOX']/1000).sum(), (CAMS_emissions_shp['NMVOC']/1000).sum(), (CAMS_emissions_shp['CO']/1000).sum(), (CAMS_emissions_shp['SO2']/1000).sum(), (CAMS_emissions_shp['NH3']/1000).sum(), (CAMS_emissions_shp['PM2_5']/1000).sum(), (CAMS_emissions_shp['PM10']/1000).sum()])








###****************************************************************************#
###**************************** Set final dataframe ***************************#
###****************************************************************************#








        #*** Replace nan with 0 ***#
        urbem_final_zero_nan = urbem_final.fillna(0)
        
        #*** Calculate the north - east coordinates per grid cell ***#

        urbem_final_zero_nan['xcor_sw'] = urbem_final_zero_nan.bounds['minx']
        urbem_final_zero_nan['xcor_ne'] = urbem_final_zero_nan.bounds['maxx']
        urbem_final_zero_nan['ycor_sw'] = urbem_final_zero_nan.bounds['miny']
        urbem_final_zero_nan['ycor_ne'] = urbem_final_zero_nan.bounds['maxy']

        #*** Add dataframe cols for final form ***#
        urbem_final_zero_nan['snap'] = urbem_final_zero_nan['SNAP']
        urbem_final_zero_nan['zcor_sw'] = 10
        urbem_final_zero_nan['zcor_ne'] = 10

        #*** Filter rows that have only zeroes for CH4, NOx, NMVOC, CO, SO2, NH3, PM2_5, PM10 ***#
        filter_zero_rows = urbem_final_zero_nan[(urbem_final_zero_nan['CH4_km'] != 0) | (urbem_final_zero_nan['NOX_km'] != 0) | (urbem_final_zero_nan['NMVOC_km'] != 0) | (urbem_final_zero_nan['CO_km'] != 0) | (urbem_final_zero_nan['SO2_km'] != 0) | (urbem_final_zero_nan['NH3_km'] != 0) | (urbem_final_zero_nan['PM2_5_km'] != 0) | (urbem_final_zero_nan['PM10_km'] != 0)]

        
        #*** Set decimal precision in each column - preparation for final table ***#

        filter_zero_rows['xcor_sw'] = filter_zero_rows['xcor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_sw'] = filter_zero_rows['ycor_sw'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['xcor_ne'] = filter_zero_rows['xcor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['ycor_ne'] = filter_zero_rows['ycor_ne'].map(lambda x: '%.0f' % x if not pd.isna(x) else '')
        filter_zero_rows['CH4'] = filter_zero_rows['CH4_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NOx'] = filter_zero_rows['NOX_km'].map(lambda x: '%.14f' % x if not pd.isna(x) else '')
        filter_zero_rows['NMVOC'] = filter_zero_rows['NMVOC_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['CO'] = filter_zero_rows['CO_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['SO2'] = filter_zero_rows['SO2_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['NH3'] = filter_zero_rows['NH3_km'].map(lambda x: '%.16f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM2.5'] = filter_zero_rows['PM2_5_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')
        filter_zero_rows['PM10'] = filter_zero_rows['PM10_km'].map(lambda x: '%.15f' % x if not pd.isna(x) else '')


        #*** Drop extra cols - Geodataframe to df ***#
        filtered_pollutants_df = filter_zero_rows.drop(['cams_index', 'pop',  'ind', 'agr', 'offroad', 'waste', 'snap34', 'ISO3', 'Year', 'SNAP', 'SourceType', 'NOX', 'PM2_5', 'area', 'snap34_val', 'snap34_val_norm', 'max_area'], axis=1)


        #*** Set col order ***#
        filtered_pollutants_df = filtered_pollutants_df[['geometry', 'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne', 'zcor_ne', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10', 'grid_index']]
        print(filtered_pollutants_df)


        #*** Set nan to -999 in final dataframe ***#
        filtered_pollutants_df.fillna(-999)

        #*** Append df to list ***#
        area_list.append(filtered_pollutants_df)

        proxy_grid_list.clear()
        cams_fine_list.clear()
        fine_temp_list.clear()
        small_cell_emission_list.clear()
        cams_no_proxy_fine_list.clear()
        small_cell_no_proxy_emission_list.clear()
        urbem_merged_list.clear()
        urbem_full_list.clear()
        small_cell_no_proxy_redistr_list.clear()
        urbem_merged_redistr_list.clear()
        smooth_list.clear()











###**************************************************************************************************************
###**************************************************************************************************************
###**************************************************************************************************************
###**************************************************************************************************************
###**************************************************************************************************************








areas_final = pd.concat(area_list)		
print(areas_final)

total_stat_list= []

areas_final_temp = areas_final.loc[areas_final['snap'] != 7]





##### ~ Statistics in Total areas tn/km*year
##areas_final_temp = areas_final.copy()
areas_final_temp = areas_final.loc[areas_final['snap'] != 7]

areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 1] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 2] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 3] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 4] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 5] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 6] = 0
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 7] = 0
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 8] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 9] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 10] = 0
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 11] = 10
areas_final_temp['zcor_sw'].loc[areas_final_temp['snap'] == 12] = 0
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 1] =1
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 2] =2
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 3] =3
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 4] =4
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 5] =5
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 6] =6
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 7] =7
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 8] =8
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 9] =1
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 10] =10
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 11] =5
areas_final_temp['snap'].loc[areas_final_temp['snap'] == 12] =10



### ~ Statistics in Total areas tn/km*year

areas_df = areas_final_temp.drop(['geometry'], axis=1)
total_areas_temp = areas_df.astype(float)
total_areas = total_areas_temp.rename(columns = {'PM2.5': 'PM25'}, inplace = False)
total_stat_list.append(['Total_area (Areas)', (total_areas['CH4']/1000).sum(), (total_areas['NOx']/1000).sum(), (total_areas['NMVOC']/1000).sum(), (total_areas['CO']/1000).sum(), (total_areas['SO2']/1000).sum(), (total_areas['NH3']/1000).sum(), (total_areas['PM25']/1000).sum(), (total_areas['PM10']/1000).sum(), len(total_areas)])





### ~ Statistics in Total areas per snap sector tn/km*year

total_stat_list.append(['Public Power', (total_areas['CH4'].loc[total_areas['snap'] == 1]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 1]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 1]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 1]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 1]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 1]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 1]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 1]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 1])])
total_stat_list.append(['Other Stationary Combustion', (total_areas['CH4'].loc[total_areas['snap'] == 2]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 2]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 2]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 2]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 2]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 2]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 2]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 2]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 2])])
total_stat_list.append(['Industry (2/3) (SNAP 3)', (total_areas['CH4'].loc[total_areas['snap'] == 3]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 3]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 3]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 3]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 3]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 3]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 3]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 3]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 3])])
total_stat_list.append(['Industry (1/3) (SNAP 4)', (total_areas['CH4'].loc[total_areas['snap'] == 4]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 4]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 4]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 4]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 4]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 4]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 4]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 4]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 4])])
total_stat_list.append(['Fugitives (SNAP 5)', (total_areas['CH4'].loc[total_areas['snap'] == 5]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 5]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 5]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 5]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 5]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 5]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 5]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 5]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 5])])
total_stat_list.append(['Solvents', (total_areas['CH4'].loc[total_areas['snap'] == 6]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 6]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 6]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 6]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 6]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 6]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 6]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 6]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 6])])
total_stat_list.append(['Shipping', (total_areas['CH4'].loc[total_areas['snap'] == 8]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 8]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 8]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 8]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 8]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 8]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 8]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 8]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 8])])
total_stat_list.append(['Waste', (total_areas['CH4'].loc[total_areas['snap'] == 9]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 9]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 9]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 9]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 9]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 9]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 9]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 9]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 9])])
total_stat_list.append(['Agriculture', (total_areas['CH4'].loc[total_areas['snap'] == 10]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 10]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 10]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 10]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 10]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 10]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 10]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 10]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 10])])
total_stat_list.append(['Aviation', (total_areas['CH4'].loc[total_areas['snap'] == 11]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 11]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 11]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 11]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 11]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 11]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 11]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 11]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 11])])
total_stat_list.append(['OffRoad', (total_areas['CH4'].loc[total_areas['snap'] == 12]/1000).sum(), (total_areas['NOx'].loc[total_areas['snap'] == 12]/1000).sum(), (total_areas['NMVOC'].loc[total_areas['snap'] == 12]/1000).sum(), (total_areas['CO'].loc[total_areas['snap'] == 12]/1000).sum(), (total_areas['SO2'].loc[total_areas['snap'] == 12]/1000).sum(), (total_areas['NH3'].loc[total_areas['snap'] == 12]/1000).sum(), (total_areas['PM25'].loc[total_areas['snap'] == 12]/1000).sum(), (total_areas['PM10'].loc[total_areas['snap'] == 12]/1000).sum(), len(total_areas.loc[total_areas['snap'] == 12])])





### ~ Clip by UC
areas_gdf = areas_final_temp.copy()
uc_cliped = cl.clip_shp(uc_4_clip, areas_gdf)
uc =  uc_cliped.drop(['grid_index', 'index_right', 'xcor', 'ycor', 'uc_val','boolean'], axis=1)

### ~ Save to shp as tn/km*year
urbem_uc = gpd.sjoin(areas_gdf, uc, how='inner', op='intersects')






### ~ Statistics in UC tn/km*year

urbem_uc_temp = urbem_uc.drop(['geometry'], axis=1)
urbem_uc_temp_fl = urbem_uc_temp.astype(float)
urbem_uc_shp = urbem_uc_temp_fl.rename(columns = {'PM2.5': 'PM25'}, inplace = False)
total_stat_list.append(['Urban_center (Areas)', (urbem_uc_shp['CH4']/1000).sum(), (urbem_uc_shp['NOx']/1000).sum(), (urbem_uc_shp['NMVOC']/1000).sum(), (urbem_uc_shp['CO']/1000).sum(), (urbem_uc_shp['SO2']/1000).sum(), (urbem_uc_shp['NH3']/1000).sum(), (urbem_uc_shp['PM25']/1000).sum(), (urbem_uc_shp['PM10']/1000).sum(), len(urbem_uc_shp)])





### ~ Statistics in UC per snap sector tn/km*year

total_stat_list.append(['Public Power (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 1]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 1])])
total_stat_list.append(['Other Stationary Combustion (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 2]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 2])])
total_stat_list.append(['Industry (2/3) (UC) (SNAP 3)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 3]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 3])])
total_stat_list.append(['Industry (1/3) (UC) (SNAP 4)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 4]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 4])])
total_stat_list.append(['Fugitives (UC) (SNAP 5)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 5]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 5])])
total_stat_list.append(['Solvents (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 6]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 6])])
total_stat_list.append(['Shipping (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 8]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 8])])
total_stat_list.append(['Waste (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 9]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 9])])
total_stat_list.append(['Agriculture (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 10]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 10])])
total_stat_list.append(['Aviation (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 11]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 11])])
total_stat_list.append(['OffRoad (UC)', (urbem_uc_shp['CH4'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), (urbem_uc_shp['NOx'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), (urbem_uc_shp['NMVOC'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), (urbem_uc_shp['CO'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), (urbem_uc_shp['SO2'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), (urbem_uc_shp['NH3'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), (urbem_uc_shp['PM25'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), (urbem_uc_shp['PM10'].loc[urbem_uc_shp['snap'] == 12]/1000).sum(), len(urbem_uc_shp.loc[urbem_uc_shp['snap'] == 12])])










### ~ Final csv tn/km*year

areas_df_4_csv = areas_final_temp.drop(['geometry','CH4', 'grid_index'], axis=1)

print(areas_df_4_csv.snap.unique())





### Check for EMPTY cells

print(np.where(areas_df_4_csv.applymap(lambda x: x == '')))
if {'snap', 'xcor_sw', 'ycor_sw', 'zcor_sw', 'xcor_ne', 'ycor_ne','zcor_ne', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10'}.issubset(areas_df_4_csv.columns):
    print("")
    print("**********************************************************************************")
    print("**********************************************************************************")
    print("WRITE csv... ")
    print("**********************************************************************************")
    print("**********************************************************************************")
    print("")
    if len(areas_df_4_csv.loc[(areas_df_4_csv['NOx'] == '') | (areas_df_4_csv['NMVOC'] == '') | (areas_df_4_csv['CO'] == '') | (areas_df_4_csv['SO2'] == '') | (areas_df_4_csv['NH3'] == '') | (areas_df_4_csv['PM2.5'] == '') | (areas_df_4_csv['PM10'] == '')]) == 0:
        ### Write csv

        areas_df_4_csv.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_areas_sources_" + str(Year) + ".csv", sep = ',', mode='w', header=True, index=False)

        ### Calculate stats

        urbem_stat_df = DataFrame(urbem_stat_list,columns=['SNAP', 'CH4(tn/year)', 'NOX(tn/year)', 'NMVOC(tn/year)', 'CO(tn/year)', 'SO2(tn/year)', 'NH3(tn/year)', 'PM25(tn/year)', 'PM10(tn/year)', 'Count'])
        cams_stat_df = DataFrame(cams_stat_list,columns=['SNAP', 'CH4(tn/year)', 'NOX(tn/year)', 'NMVOC(tn/year)', 'CO(tn/year)', 'SO2(tn/year)', 'NH3(tn/year)', 'PM25(tn/year)', 'PM10(tn/year)'])
        total_stat_df = DataFrame(total_stat_list,columns=['Snap/Region', 'CH4(tn/year)', 'NOx(tn/year)', 'NMVOC(tn/year)', 'CO(tn/year)', 'SO2(tn/year)', 'NH3(tn/year)', 'PM25(tn/year)', 'PM10(tn/year)', 'Count'])
        print(urbem_stat_df)
        print(cams_stat_df)
        print(total_stat_df)

        ### Write stats csv

        urbem_stat_df.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_urbem_stat_areas_sources_" + str(Year) + ".csv", sep = ';', mode='w', header=True, index=False)
        cams_stat_df.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_cams_stat_areas_sources_" + str(Year) + ".csv", sep = ';', mode='w', header=True, index=False)
        total_stat_df.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_total_stat_areas_sources_" + str(Year) + ".csv", sep = ';', mode='w', header=True, index=False)
else:
    print("NOT the write format: 1) Check cols, 2)Check for empty '' cells 3) Check if you drop the CH4 col")
















#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
















#**************************************************************************************************************
#**************************************************************************************************************
#*************************************** Calculate Line Source Emissions ***************************************#
#**************************************************************************************************************
#**************************************************************************************************************
















#*** URBEM ***#


areas_final = pd.concat(area_list)		
print(areas_final)

areas_final_temp = areas_final.loc[areas_final['snap'] == 7]





urbem_stat_df = DataFrame(urbem_stat_list,columns=['SNAP', 'CH4(tn/year)', 'NOX(tn/year)', 'NMVOC(tn/year)', 'CO(tn/year)', 'SO2(tn/year)', 'NH3(tn/year)', 'PM25(tn/year)', 'PM10(tn/year)', 'Count'])
cams_stat_df = DataFrame(cams_stat_list,columns=['SNAP', 'CH4(tn/year)', 'NOX(tn/year)', 'NMVOC(tn/year)', 'CO(tn/year)', 'SO2(tn/year)', 'NH3(tn/year)', 'PM25(tn/year)', 'PM10(tn/year)'])
print(urbem_stat_df)
print(cams_stat_df)










#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************











#**************************************************************************************************************
#*************************************** GHS urban core proxy increase ***************************************#
#**************************************************************************************************************



ghs = gdal.Open(Proxy_Folder + "/ghs_incr_fact_" + str(uc_increase_factor) + ".tif")
ghs_table = ghs.GetRasterBand(1)
ghs_geo_transform = ghs.GetGeoTransform()
ghs_np_table = ghs_table.ReadAsArray()

ghs_np_table[ghs_np_table == 0]= 1
ghs_np_fl = ghs_np_table.flatten()
ghs_norm_df = pd.DataFrame(data=ghs_np_fl, columns=["boolean"])



areas_df = areas_final_temp.drop(['geometry'], axis=1)
cams_all_traffic = areas_df.astype(float)



total_stat_list = []
areas_final_temp2 = cams_all_traffic.copy()
areas_final_temp2['zcor_sw'].loc[areas_final_temp2['snap'] == 7] = 0
areas_final_temp2['snap'].loc[areas_final_temp2['snap'] == 7] =7



### ~ Statistics in Total areas tn/km*year

total_areas_temp = areas_df.astype(float)
total_areas = total_areas_temp.rename(columns = {'PM2.5': 'PM25'}, inplace = False)
total_stat_list.append(['Total_area  - Before Increase (Lines)', (total_areas['CH4']/1000).sum(), (total_areas['NOx']/1000).sum(), (total_areas['NMVOC']/1000).sum(), (total_areas['CO']/1000).sum(), (total_areas['SO2']/1000).sum(), (total_areas['NH3']/1000).sum(), (total_areas['PM25']/1000).sum(), (total_areas['PM10']/1000).sum(), len(total_areas)])





### ~ Clip by UC
uc_areas_final_temp_list = []



uc_grid_index = uc_clip.grid_index.unique()



for index in range(0, len(uc_grid_index)):
    uc_areas_final_temp_list.append(total_areas.loc[total_areas['grid_index'] == uc_grid_index[index]])
urbem_uc = pd.concat(uc_areas_final_temp_list)
print(urbem_uc)





urbem_uc_shp = urbem_uc.rename(columns = {'PM2.5': 'PM25'}, inplace = False)
total_stat_list.append(['Urban_center  - Before Increase (Lines)', (urbem_uc_shp['CH4']/1000).sum(), (urbem_uc_shp['NOx']/1000).sum(), (urbem_uc_shp['NMVOC']/1000).sum(), (urbem_uc_shp['CO']/1000).sum(), (urbem_uc_shp['SO2']/1000).sum(), (urbem_uc_shp['NH3']/1000).sum(), (urbem_uc_shp['PM25']/1000).sum(), (urbem_uc_shp['PM10']/1000).sum(), len(urbem_uc_shp)])





print("cams_all_traffic NOx sum in kt/year", cams_all_traffic.NOx.sum()/1000/1000)





cams_nox_increase_list = []
cams_all_increase_list = []
for i in range (0,len(SNAP_sectors)):
    cams_all_traffic_per_snap = cams_all_traffic.loc[cams_all_traffic['snap'] == SNAP_sectors[i]]

    cams_ghs_join = pd.concat([cams_all_traffic_per_snap, ghs_norm_df], axis=1, join='inner')

    cams_ghs_nox_increase = cams_ghs_join.copy()

    cams_ghs_nox_increase.loc[:,'NOx'] *= cams_ghs_nox_increase.loc[:,'boolean']

    cams_ghs_nox_increase_final = cams_ghs_nox_increase#.drop(['boolean'], axis=1)
    cams_nox_increase_list.append(cams_ghs_nox_increase)

    cams_ghs_all_increase = cams_ghs_join.copy()
    cams_ghs_all_increase.loc[:,'CH4'] *= cams_ghs_all_increase.loc[:,'boolean']
    cams_ghs_all_increase.loc[:,'NOx'] *= cams_ghs_all_increase.loc[:,'boolean']
    cams_ghs_all_increase.loc[:,'NMVOC'] *= cams_ghs_all_increase.loc[:,'boolean']
    cams_ghs_all_increase.loc[:,'CO'] *= cams_ghs_all_increase.loc[:,'boolean']
    cams_ghs_all_increase.loc[:,'SO2'] *= cams_ghs_all_increase.loc[:,'boolean']
    cams_ghs_all_increase.loc[:,'NH3'] *= cams_ghs_all_increase.loc[:,'boolean']
    cams_ghs_all_increase.loc[:,'PM2.5'] *= cams_ghs_all_increase.loc[:,'boolean']
    cams_ghs_all_increase.loc[:,'PM10'] *= cams_ghs_all_increase.loc[:,'boolean']

    cams_ghs_all_increase_final = cams_ghs_all_increase#.drop(['boolean'], axis=1)
    cams_all_increase_list.append(cams_ghs_all_increase)

cams_nox_increase = pd.concat(cams_nox_increase_list)
cams_all_traffic_nox_increase = cams_nox_increase.drop(['boolean'], axis=1)

cams_all_increase = pd.concat(cams_all_increase_list)
cams_all_traffic_all_increase = cams_all_increase.drop(['boolean'], axis=1)
print(cams_all_traffic_all_increase)





print("cams traffic NOx increase sum in kt/year", cams_all_traffic_nox_increase.NOx.sum()/1000/1000)
print("cams traffic all increase sum in kt/year", cams_all_traffic_all_increase.NOx.sum()/1000/1000)





emis = cams_all_traffic.groupby(['xcor_sw','ycor_sw','xcor_ne', 'ycor_ne', 'grid_index'], as_index=False).sum()
emis_nox = cams_all_traffic_nox_increase.groupby(['xcor_sw','ycor_sw','xcor_ne', 'ycor_ne',  'grid_index'], as_index=False).sum()
emis_all = cams_all_traffic_all_increase.groupby(['xcor_sw','ycor_sw','xcor_ne', 'ycor_ne',  'grid_index'], as_index=False).sum()




emis_all_temp = emis_all.rename(columns = {'PM2.5': 'PM25'}, inplace = False)
total_stat_list.append(['Total_area  - After Increase (Lines)', (emis_all_temp['CH4']/1000).sum(), (emis_all_temp['NOx']/1000).sum(), (emis_all_temp['NMVOC']/1000).sum(), (emis_all_temp['CO']/1000).sum(), (emis_all_temp['SO2']/1000).sum(), (emis_all_temp['NH3']/1000).sum(), (emis_all_temp['PM25']/1000).sum(), (emis_all_temp['PM10']/1000).sum(), len(emis_all_temp)])





uc_emis_all_list = []

for index in range(0, len(uc_grid_index)):
    uc_emis_all_list.append(emis_all_temp.loc[emis_all_temp['grid_index'] == uc_grid_index[index]])
uc_emis_all = pd.concat(uc_emis_all_list)
print(uc_emis_all)



total_stat_list.append(['Urban_center - After Increase (Lines)', (uc_emis_all['CH4']/1000).sum(), (uc_emis_all['NOx']/1000).sum(), (uc_emis_all['NMVOC']/1000).sum(), (uc_emis_all['CO']/1000).sum(), (uc_emis_all['SO2']/1000).sum(), (uc_emis_all['NH3']/1000).sum(), (uc_emis_all['PM25']/1000).sum(), (uc_emis_all['PM10']/1000).sum(), len(uc_emis_all)])



total_stat_df = DataFrame(total_stat_list,columns=['Snap/Region', 'CH4(tn/year)', 'NOx(tn/year)', 'NMVOC(tn/year)', 'CO(tn/year)', 'SO2(tn/year)', 'NH3(tn/year)', 'PM25(tn/year)', 'PM10(tn/year)', 'Count'])
print(total_stat_df)


#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************











#**************************************************************************************************************
#************************************************* OSM proxy ***************************************************#
#**************************************************************************************************************



osm = gpd.read_file(Proxy_Folder + "/osm_utm34.shp")

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************











#**************************************************************************************************************
#****************************************** Distribute Emissions to line **********************************************#
#**************************************************************************************************************



grid['grid_index'] = grid.index

grid_index = grid.grid_index.unique()

osm_grid_sjoin = gpd.sjoin(osm, grid, how='inner', op='intersects')
osm_grid_sjoin['length'] = osm_grid_sjoin.geometry.length
osm_grid_sjoin.to_file(OutFolder + "osm_grid_sjoin.shp", driver="ESRI Shapefile")

osm_grid_index = osm_grid_sjoin.grid_index.unique()






osm_cell_all_list = []


for i in range (0,len(grid)):
    print(i)
    print(int((i / len(grid))*100), "% done")


    osm_cell = osm_grid_sjoin.loc[osm_grid_sjoin["grid_index"] == i]
    osm_cell['Highway'] = osm_cell['typeOfRoad']
    osm_cell.loc[(osm_cell['Highway'] == 'primary') | (osm_cell['Highway'] == 'primary_link'), 'Highway'] = 'primary' 
    osm_cell.loc[(osm_cell['Highway'] == 'secondary') | (osm_cell['Highway'] == 'secondary_link'), 'Highway'] = 'secondary' 
    osm_cell.loc[(osm_cell['Highway'] == 'trunk') | (osm_cell['Highway'] == 'trunk_link'), 'Highway'] = 'trunk'
    osm_cell.loc[(osm_cell['Highway'] == 'motorway') | (osm_cell['Highway'] == 'motorway_link'), 'Highway'] = 'motorway'



    st = osm_cell.Highway.unique()

    st_df = pd.DataFrame(st, columns=["roadtype"])


    road_weights_df = pd.DataFrame(0, index=np.arange(4), columns=["weights"])

    if len(st) > 0:
        for j in range (0,len(st_df)):
            r_type = st_df['roadtype'].iloc[j]

            if r_type == "motorway":
                road_weights_df["weights"].iloc[0] = 10
            elif r_type == "trunk":
                road_weights_df["weights"].iloc[1] = 5
            elif r_type == "primary":
                road_weights_df["weights"].iloc[2] = 2
            elif r_type == "secondary":
                road_weights_df["weights"].iloc[3] = 2

        street_sum = road_weights_df['weights'].sum()

        road_weights_df['norm_weights'] = road_weights_df['weights']/street_sum

        normalized_road_weights = road_weights_df#.drop(['weights'], axis=1)

        osm_cell_final = osm_cell[['length','Highway']]

        if len(emis_all.loc[emis_all['grid_index'] == i]) == 0:
            print("passsssssss...")

        else:
            CH4 = float(emis_all['CH4'].loc[emis_all['grid_index'] == i])
            NOx = float(emis_all['NOx'].loc[emis_all['grid_index'] == i])
            NMVOC = float(emis_all['NMVOC'].loc[emis_all['grid_index'] == i])
            CO = float(emis_all['CO'].loc[emis_all['grid_index'] == i])
            SO2 = float(emis_all['SO2'].loc[emis_all['grid_index'] == i])
            NH3 = float(emis_all['NH3'].loc[emis_all['grid_index'] == i])
            PM2_5 = float(emis_all['PM2.5'].loc[emis_all['grid_index'] == i])
            PM10 = float(emis_all['PM10'].loc[emis_all['grid_index'] == i])
            print(CH4)



            osm_cell_m = osm_cell.loc[osm_cell['Highway'] == "motorway"]
            norm_weight_m = normalized_road_weights['norm_weights'].iloc[0]
            
            osm_cell_m['CH4'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * CH4
            osm_cell_m['NOx'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * NOx
            osm_cell_m['NMVOC'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * NMVOC
            osm_cell_m['CO'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * CO
            osm_cell_m['SO2'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * SO2
            osm_cell_m['NH3'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * NH3
            osm_cell_m['PM2.5'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * PM2_5
            osm_cell_m['PM10'] = osm_cell_m['length']/osm_cell_m['length'].sum() * norm_weight_m * PM10
            osm_cell_all_list.append(osm_cell_m)


            osm_cell_t = osm_cell.loc[osm_cell['Highway'] == "trunk"]
            norm_weight_t = normalized_road_weights['norm_weights'].iloc[1]
            osm_cell_t['CH4'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * CH4
            osm_cell_t['NOx'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * NOx
            osm_cell_t['NMVOC'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * NMVOC
            osm_cell_t['CO'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * CO
            osm_cell_t['SO2'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * SO2
            osm_cell_t['NH3'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * NH3
            osm_cell_t['PM2.5'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * PM2_5
            osm_cell_t['PM10'] = osm_cell_t['length']/osm_cell_t['length'].sum() * norm_weight_t * PM10
            osm_cell_all_list.append(osm_cell_t)


            osm_cell_p = osm_cell.loc[osm_cell['Highway'] == "primary"]
            norm_weight_p = normalized_road_weights['norm_weights'].iloc[2]
            osm_cell_p['CH4'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * CH4
            osm_cell_p['NOx'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * NOx
            osm_cell_p['NMVOC'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * NMVOC
            osm_cell_p['CO'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * CO
            osm_cell_p['SO2'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * SO2
            osm_cell_p['NH3'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * NH3
            osm_cell_p['PM2.5'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * PM2_5
            osm_cell_p['PM10'] = osm_cell_p['length']/osm_cell_p['length'].sum() * norm_weight_p * PM10
            osm_cell_all_list.append(osm_cell_p)


            osm_cell_s = osm_cell.loc[osm_cell['Highway'] == "secondary"]
            norm_weight_s = normalized_road_weights['norm_weights'].iloc[3]
            osm_cell_s['CH4'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * CH4
            osm_cell_s['NOx'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * NOx
            osm_cell_s['NMVOC'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * NMVOC
            osm_cell_s['CO'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * CO
            osm_cell_s['SO2'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * SO2
            osm_cell_s['NH3'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * NH3
            osm_cell_s['PM2.5'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * PM2_5
            osm_cell_s['PM10'] = osm_cell_s['length']/osm_cell_s['length'].sum() * norm_weight_s * PM10
            osm_cell_all_list.append(osm_cell_s)
    else:
        print("pass:")








#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************











#**************************************************************************************************************
#****************************************** Final Calculations **********************************************************#
#**************************************************************************************************************



osm_cell_temp = pd.concat(osm_cell_all_list)		

osm_cell_all = osm_cell_temp.drop([('ID'), ('key'), ('Highway'), ('index_right'), ('grid_index')], axis=1)


out = osm_cell_all.copy()
print("Sum NOx after process: ", out['NOx'].sum())
print("CAMS sum NOx before process: ", emis_all['NOx'].sum())

out['CH4'] = out['CH4']*emis_all['CH4'].sum()/out['CH4'].sum()
out['NOx'] = out['NOx']*emis_all['NOx'].sum()/out['NOx'].sum()
out['NMVOC'] = out['NMVOC']*emis_all['NMVOC'].sum()/out['NMVOC'].sum()
out['CO'] = out['CO']*emis_all['CO'].sum()/out['CO'].sum()
out['SO2'] = out['SO2']*emis_all['SO2'].sum()/out['SO2'].sum()
out['NH3'] = out['NH3']*emis_all['NH3'].sum()/out['NH3'].sum()
out['PM2.5'] = out['PM2.5']*emis_all['PM2.5'].sum()/out['PM2.5'].sum()
out['PM10'] = out['PM10']*emis_all['PM10'].sum()/out['PM10'].sum()
print("Correct lost emissions - Sum NOx after process: ", out['NOx'].sum())
print("CAMS sum NOx before process: ", emis_all['NOx'].sum())
print("Emissions length == CAMS length: ", out['length'].sum() == osm_grid_sjoin['length'].sum())
print(out[['CH4', 'CO', 'PM10']])

out['width'] = 0
out['width'].loc[out['typeOfRoad'] == 'motorway'] = 20
out['width'].loc[out['typeOfRoad'] == 'motorway_link'] = 20
out['width'].loc[out['typeOfRoad'] == 'trunk'] = 16
out['width'].loc[out['typeOfRoad'] == 'trunk_link'] = 16
out['width'].loc[out['typeOfRoad'] == 'primary'] = 12
out['width'].loc[out['typeOfRoad'] == 'primary_link'] = 12
out['width'].loc[out['typeOfRoad'] == 'secondary'] = 12
out['width'].loc[out['typeOfRoad'] == 'secondary_link'] = 12
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************











#**************************************************************************************************************
#****************************************** Save lines to shp **********************************************************#
#**************************************************************************************************************



print(out[['CH4', 'CO', 'PM10']])


out = out.loc[out['geometry'] != None].reset_index(drop=True)
filter_zero_rows_final = out[(out['CH4'] != 0) | (out['NOx'] != 0) | (out['NMVOC'] != 0) | (out['CO'] != 0) | (out['SO2'] != 0) | (out['NH3'] != 0) | (out['PM2.5'] != 0) | (out['PM10'] != 0)]
out = filter_zero_rows_final.copy()


### ~ Save to shp as tn/km*year
out_lines = out.copy()
out_lines['CH4'] = out_lines['CH4']/1000
out_lines['NOx'] = out_lines['NOx']/1000
out_lines['NMVOC'] = out_lines['NMVOC']/1000
out_lines['CO'] = out_lines['CO']/1000
out_lines['SO2'] = out_lines['SO2']/1000
out_lines['NH3'] = out_lines['NH3']/1000
out_lines['PM2.5'] = out_lines['PM2.5']/1000
out_lines['PM10'] = out_lines['PM10']/1000
lines_shp = out_lines.rename(columns = {'PM2.5': 'PM25'}, inplace = False)
lines_shp.to_file(OutFolder + "lines.shp", driver="ESRI Shapefile")
print(len(out))









out['xcor_start'] = 0
out['ycor_start'] = 0
out['xcor_end'] = 0
out['ycor_end'] = 0

for k in range (0,len(out)):
    out['xcor_start'].iloc[k] = int(out['geometry'].iloc[k].bounds[0])
    out['ycor_start'].iloc[k] = int(out['geometry'].iloc[k].bounds[1])
    out['xcor_end'].iloc[k] = int(out['geometry'].iloc[k].bounds[2])
    out['ycor_end'].iloc[k] = int(out['geometry'].iloc[k].bounds[3])


print(out[['CH4', 'CO', 'PM10']])





#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************











#**************************************************************************************************************
#************************************************* Final df for csv *********************************************
#**************************************************************************************************************



out['snap'] = 7
out['elevation'] = 0
out['CH4'] = out['CH4'] * 1000/(365 * 24 * 3600)
out['NOx'] = out['NOx'] * 1000/(365 * 24 * 3600)
out['NMVOC'] = out['NMVOC'] * 1000/(365 * 24 * 3600)
out['CO'] = out['CO'] * 1000/(365 * 24 * 3600)
out['SO2'] = out['SO2'] * 1000/(365 * 24 * 3600)
out['NH3'] = out['NH3'] * 1000/(365 * 24 * 3600)
out['PM2.5'] = out['PM2.5'] * 1000/(365 * 24 * 3600)
out['PM10'] = out['PM10'] * 1000/(365 * 24 * 3600)
out.fillna(-999)
out_final = out[['snap', 'xcor_start', 'ycor_start', 'xcor_end', 'ycor_end', 'elevation', 'width', 'CH4', 'NOx', 'NMVOC', 'CO', 'SO2', 'NH3', 'PM2.5', 'PM10']]
print(out_final)

out_final_no_ch4 = out_final.drop(['CH4'], axis=1)

### Check for EMPTY cells

print(np.where(out_final_no_ch4.applymap(lambda x: x == '')))

### Save to csv

out_final_no_ch4.to_csv(csvFolder + "Nasia_" +Region + "_CAMS_v3_1_lines_sources_" + str(Year) + "_all_increase.csv", sep = ',', mode='w', header=True, index=False)
urbem_stat_df.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_urbem_stat_lines_sources_" + str(Year) + ".csv", sep = ';', mode='w', header=True, index=False)
cams_stat_df.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_cams_stat_lines_sources_" + str(Year) + ".csv", sep = ';', mode='w', header=True, index=False)
total_stat_df.to_csv(csvFolder + "Nasia_" + Region + "_CAMS_v3_1_total_stat_lines_sources_" + str(Year) + ".csv", sep = ';', mode='w', header=True, index=False)








#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************






elapsed_time = time.time() - start_time
hours, rem = divmod(elapsed_time, 3600)
minutes, seconds = divmod(rem, 60)
print("minutes", minutes)
print("seconds", elapsed_time)

##print("--- %s seconds ---" % (time.time() - start_time))

