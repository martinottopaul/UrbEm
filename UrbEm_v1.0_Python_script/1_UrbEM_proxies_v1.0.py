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
import shutil
from OSMPythonTools.nominatim import Nominatim
from OSMPythonTools.overpass import Overpass, overpassQueryBuilder
from OSMPythonTools.data import Data, dictRangeYears, ALL
from collections import OrderedDict
import shapefile
from rasterio import features


start_time = time.time()


### Set working dir & create folders

theWD = os.path.dirname(inspect.getfile(inspect.currentframe())) + "\\"
theWD = str(theWD.replace("\\", "/"))
print(theWD)

##################################################
##InFolder = sys.argv[1]
##OutFolder = sys.argv[2] + "\\"
##xmin = float(sys.argv[3])
##xmax = float(sys.argv[4])
##ymin = float(sys.argv[5])
##ymax = float(sys.argv[6])
##cell_size = int(sys.argv[7])
##crs_utm = 'EPSG:' + sys.argv[8]
##crs_wgs = 'EPSG:4326'
##crs_corine = 'EPSG:3035'
##uc_increase_factor = int(sys.argv[7])
##################################################

################################################
InFolder = theWD + "Input_Data"
OutFolder = theWD + "Athens/Results_Folder/Proxies/"
xmin = 716397
xmax = 761397
ymin = 4191261
ymax = 4236261
cell_size = 1000
epsg_code = 32634
crs_utm = 'EPSG:32634'
crs_wgs = 'EPSG:4326'
crs_corine = 'EPSG:3035'
uc_increase_factor = 1
Population_Density_Raster = InFolder + "/GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif"
Corine_Raster = InFolder + "/83684d24c50f069b613e0dc8e12529b893dc172f/u2018_clc2018_v2020_20u1_raster100m/DATA/U2018_CLC2018_V2020_20u1.tif"
Corine_GDB = InFolder + "/clc2018_clc2018_v2018_20_fgdb/CLC2018_CLC2018_V2018_20.gdb"
Corine_GDB_layer = fiona.listlayers(InFolder + '/clc2018_clc2018_v2018_20_fgdb/CLC2018_CLC2018_V2018_20.gdb')
E_PRTR_kmz = InFolder + '/E_PRTR/kmz/EPRTR_facilities_v17.kmz'
GHS_UC_shp = InFolder + "/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0/GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_0.shp"
Shipping_Routes_shp = InFolder + "/Shipping_Routes/shipping_routes_wgs.shp"

################################################


if not os.path.exists(OutFolder):
    print("Data Folder doesn't exist... Creating Folder...")
    os.makedirs(OutFolder)

elif os.path.exists(OutFolder):
    creation_time = os.path.getctime(OutFolder)

    if (start_time - creation_time) // (24 * 3600) >= 31:
        print("Data Folder has been created > 1 month ago... Deleting Folder...")
        shutil.rmtree(OutFolder)
        shutil.rmtree(theWD + "cache/")
        os.makedirs(OutFolder)
    elif len(os.listdir(OutFolder)) == 0:
        print("Data Folder is empty...")
    else:
        print("Data Folder exists or has been created < 1 month ago... Creating emissions...")
##        sys.exit()


print("Starting process...")








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#
















# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#
















# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#








#**************************************************************************************************************
#*************************************** Create small grid for domain ***************************************#
#**************************************************************************************************************



### Set cell size and resolution - Coordinates and Grid Parameters
# For proxies

### Hamburg
# ~ ymin = 5918656
# ~ ymax = 5948656
# ~ xmin = 551750
# ~ xmax = 581750

lat_point_list = [ymin, ymin, ymax, ymax]
lon_point_list = [xmin, xmax, xmax, xmin]

rows = int(np.ceil((ymax-ymin) /  cell_size))
cols = int(np.ceil((xmax-xmin) / cell_size))

proxy_width = rows
proxy_height = cols
print(proxy_width)
print(proxy_height)
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
##grid_wgs.to_file(OutFolder + "grid_wgs.shp", driver="ESRI Shapefile")


# Create point form polygon - grid centroid coordinates
grid_points = grid.copy()

grid_points['geometry'] = grid_points['geometry'].centroid
grid_points['xcor'] = grid_points['geometry'].x
grid_points['ycor'] = grid_points['geometry'].y
##grid_points.to_file(OutFolder + "grid_points.shp", driver="ESRI Shapefile")
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



grid_corine_bounds = grid_corine.bounds
xmin_corine = grid_corine_bounds['minx'].min()
xmax_corine = grid_corine_bounds['maxx'].max()
ymin_corine = grid_corine_bounds['miny'].min()
ymax_corine = grid_corine_bounds['maxy'].max()

print(xmin_corine)
print(xmax_corine)
print(ymin_corine)
print(ymax_corine)

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








###**************************************************************************************************************
###*************************************** Preparation of area proxies ***************************************#
###**************************************************************************************************************



# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#


print("Creating Population Density proxy...")

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#








###**************************************************************************************************************
###*************************************** 1) Population Density proxy ***************************************#
###**************************************************************************************************************



pop = gdal.Open(Population_Density_Raster)

# Project Raster
pop_proj = OutFolder + "/pop_proj.tif"
gdal.Warp(pop_proj, pop, format= "GTiff", outputBounds=[xmin, ymin, xmax, ymax], xRes= cell_size, yRes= cell_size, dstSRS=crs_utm, resampleAlg = gdal.GRA_Mode , srcNodata= -200, dstNodata = -9999, cutlineLayer= grid, cropToCutline= True, copyMetadata=True)

# Read croped raster and reclassify values < 0  to -9999
# Read croped raster

filehandle = gdal.Open(OutFolder + "pop_proj.tif")
attribute_table = filehandle.GetRasterBand(1)
meta = filehandle.GetMetadata()
geotransform = filehandle.GetGeoTransform()
geoproj = filehandle.GetProjection()
np_table = attribute_table.ReadAsArray()
xsize = filehandle.RasterXSize
ysize = filehandle.RasterYSize



###**************************************************************************************************************
###**************************************************************************************************************
###**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#








#**************************************************************************************************************
#*************************************** Reclassify Corine Raster 2018 ***************************************#
#**************************************************************************************************************



print("Creating Corine proxies...")




driver = gdal.GetDriverByName('GTiff')
Corine = gdal.Open(Corine_Raster)


band = Corine.GetRasterBand(1)
lista = band.ReadAsArray()

corine_np = lista.copy()

corine_utm = OutFolder + "corine_utm.tif"

gdal.Warp(corine_utm, Corine, format= "GTiff", outputBounds=[xmin, ymin, xmax, ymax], xRes= cell_size, yRes= cell_size, dstSRS=crs_utm, dstNodata = -9999, cutlineLayer= grid, cropToCutline= True, copyMetadata=True)
driver = gdal.GetDriverByName('GTiff')
Corine_masked = gdal.Open(corine_utm)
band_masked = Corine_masked.GetRasterBand(1)
list_masked = band_masked.ReadAsArray()

corine_np = list_masked.copy()







# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#






###### reclassification based on raster value

corine_np[np.where((list_masked == 1) | (list_masked == 2) | (list_masked == 10) | (list_masked == 11))] = 1
corine_np[np.where((list_masked == 3))] = 2
corine_np[np.where((list_masked == 4))] = 3
corine_np[np.where((list_masked == 5))] = 4
corine_np[np.where((list_masked == 6))] = 5
corine_np[np.where((list_masked == 12) | (list_masked == 13) | (list_masked == 14) | (list_masked == 15) | (list_masked == 16) | (list_masked == 17) | (list_masked == 18) | (list_masked == 19) | (list_masked == 20) | (list_masked == 21) | (list_masked == 22))] = 6
corine_np[np.where((list_masked == 23) | (list_masked == 24) | (list_masked == 25) | (list_masked == 26) | (list_masked == 27) | (list_masked == 28) | (list_masked == 29))] = 7
corine_np[np.where((list_masked == 30) | (list_masked == 31) | (list_masked == 32) | (list_masked == 33))] = 8
corine_np[np.where((list_masked == 34))] = 9
corine_np[np.where((list_masked == 35) | (list_masked == 36) | (list_masked == 40) | (list_masked == 41))] = 10
corine_np[np.where((list_masked == 37) | (list_masked == 38) | (list_masked == 39))] = 11
corine_np[np.where((list_masked == 42) | (list_masked == 43) | (list_masked == 44))] = 12
corine_np[np.where((list_masked == 7) | (list_masked == 8) | (list_masked == 9))] = 13
##corine_np[np.where((list_masked == 7) | (list_masked == 8))] = 13
##corine_np[np.where((list_masked == 9))] = 14
# ~ #corine_np[np.where((list_masked == -32768))] = np.nan






# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#







#*** Airport ***#


Airport = corine_np.copy()
Airport[Airport == 5] = 99
Airport[Airport != 99] = 0
Airport[Airport == 99] = 1




air_filename = OutFolder + 'lu_airport.tif'
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg_code)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(air_filename,proxy_height,proxy_width,1,gdal.GDT_Float64)
dataset.SetGeoTransform((xmin,cell_size,0,ymax,0,-cell_size))
dataset.SetProjection(srs.ExportToWkt())
dataset.GetRasterBand(1).WriteArray(Airport)
dataset.FlushCache()  








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#






#*** Mineral Dump/Dump/Construction Sites ***#


MineralDump_Dump_Construction = corine_np.copy()
MineralDump_Dump_Construction[MineralDump_Dump_Construction == 13] = 99
MineralDump_Dump_Construction[MineralDump_Dump_Construction != 99] = 0
MineralDump_Dump_Construction[MineralDump_Dump_Construction == 99] = 1








#*** Industry ***#


Industry = corine_np.copy()
Industry[Industry == 2] = 99
Industry[Industry != 99] = 0
Industry[Industry == 99] = 1



ind_filename = OutFolder + 'lu_industry.tif'
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg_code)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(ind_filename,proxy_height,proxy_width,1,gdal.GDT_Float64)
dataset.SetGeoTransform((xmin,cell_size,0,ymax,0,-cell_size))
dataset.SetProjection(srs.ExportToWkt())
dataset.GetRasterBand(1).WriteArray(Industry)
dataset.FlushCache()  






# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#







#*** Agriculture ***#


Agriculture = corine_np.copy()
Agriculture[Agriculture == 6] = 99
Agriculture[Agriculture != 99] = 0
Agriculture[Agriculture == 99] = 1



agr_filename = OutFolder + 'lu_agriculture.tif'
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg_code)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(agr_filename,proxy_height,proxy_width,1,gdal.GDT_Float64)
dataset.SetGeoTransform((xmin,cell_size,0,ymax,0,-cell_size))
dataset.SetProjection(srs.ExportToWkt())
dataset.GetRasterBand(1).WriteArray(Agriculture)
dataset.FlushCache()  








#*** Ports ***#


Ports = corine_np.copy()
Ports[Ports == 4] = 99
Ports[Ports != 99] = 0
Ports[Ports == 99] = 1



ports_filename = OutFolder + 'lu_ports.tif'
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg_code)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(ports_filename,proxy_height,proxy_width,1,gdal.GDT_Float64)
dataset.SetGeoTransform((xmin,cell_size,0,ymax,0,-cell_size))
dataset.SetProjection(srs.ExportToWkt())
dataset.GetRasterBand(1).WriteArray(Ports)
dataset.FlushCache()  








#*** Sea/Ocean ***#


Sea_Ocean = corine_np.copy()
Sea_Ocean[Sea_Ocean == 12] = 99
Sea_Ocean[Sea_Ocean != 99] = 0
Sea_Ocean[Sea_Ocean == 99] = 1








#*** Salines ***#


Salines = corine_np.copy()
Salines[Salines == 11] = 99
Salines[Salines != 99] = 0
Salines[Salines == 99] = 1








#*** Bare Soil/Beaches/Burnt Areas ***#


BareSoil_Beaches_BurntAreas = corine_np.copy()
BareSoil_Beaches_BurntAreas[BareSoil_Beaches_BurntAreas == 8] = 99
BareSoil_Beaches_BurntAreas[BareSoil_Beaches_BurntAreas != 99] = 0
BareSoil_Beaches_BurntAreas[BareSoil_Beaches_BurntAreas == 99] = 1








#*** Urban Areas ***#


Urban = corine_np.copy()
Urban[Urban == 1] = 99
Urban[Urban != 99] = 0
Urban[Urban == 99] = 1



urban_filename = OutFolder + 'lu_urban.tif'
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg_code)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(urban_filename,proxy_height,proxy_width,1,gdal.GDT_Float64)
dataset.SetGeoTransform((xmin,cell_size,0,ymax,0,-cell_size))
dataset.SetProjection(srs.ExportToWkt())
dataset.GetRasterBand(1).WriteArray(Urban)
dataset.FlushCache()  








#*** Water ***#


Water = corine_np.copy()
Water[Water == 10] = 99
Water[Water != 99] = 0
Water[Water == 99] = 1








#*** Glaciers ***#


Glaciers = corine_np.copy()
Glaciers[Glaciers == 9] = 99
Glaciers[Glaciers != 99] = 0
Glaciers[Glaciers == 99] = 1








#*** Forest ***#


Forest = corine_np.copy()
Forest[Forest == 7] = 99
Forest[Forest != 99] = 0
Forest[Forest == 99] = 1








#*** Snap 8 from CORINE ***#


corine_snap8 = corine_np.copy()
corine_snap8[(corine_snap8 == 2) | (corine_snap8 == 4) | (corine_snap8 == 5) | (corine_snap8 == 6) | (corine_snap8 == 13)] = 99
corine_snap8[corine_snap8 != 99] = 0
corine_snap8[corine_snap8 == 99] = 1



snap8_filename = OutFolder + 'lu_snap8_temp.tif'
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg_code)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(snap8_filename,proxy_height,proxy_width,1,gdal.GDT_Float64)
dataset.SetGeoTransform((xmin,cell_size,0,ymax,0,-cell_size))
dataset.SetProjection(srs.ExportToWkt())
dataset.GetRasterBand(1).WriteArray(corine_snap8)
dataset.FlushCache()  







# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#






#*** Offroad from CORINE ***#


offroad = corine_np.copy()
offroad[(offroad == 2) | (offroad == 6) | (offroad == 13)] = 99
offroad[offroad != 99] = 0
offroad[offroad == 99] = 1



offroad_filename = OutFolder + 'lu_offroad.tif'
srs = osr.SpatialReference()
srs.ImportFromEPSG(epsg_code)
driver = gdal.GetDriverByName('GTiff')
dataset = driver.Create(offroad_filename,proxy_height,proxy_width,1,gdal.GDT_Float64)
dataset.SetGeoTransform((xmin,cell_size,0,ymax,0,-cell_size))
dataset.SetProjection(srs.ExportToWkt())
dataset.GetRasterBand(1).WriteArray(offroad)
dataset.FlushCache()  







# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#






#*** Ports - Raster to polygon ***#


port_raster = rasterio.open(OutFolder + "lu_ports.tif")
print(port_raster.crs)
print(port_raster.count)

port_pr_x_list = []
port_pr_y_list = []
port_pr_val_list = []

#extract point value from raster
for point in grid_points['geometry']:
    x = point.xy[0][0]
    y = point.xy[1][0]
    row, col = port_raster.index(x,y)
    val = port_raster.read(1)[row,col]
    port_pr_x_list.append(x)
    port_pr_y_list.append(y)
    port_pr_val_list.append(val)
    
port_norm_post_proc_df = pd.DataFrame(list(zip(port_pr_x_list, port_pr_y_list, port_pr_val_list)), columns =['xcor', 'ycor', 'port_val']) 
port_raster_to_points_gdf = gpd.GeoDataFrame(port_norm_post_proc_df, geometry=gpd.points_from_xy(x=port_norm_post_proc_df.xcor, y=port_norm_post_proc_df.ycor), crs = crs_utm)
port_raster_to_points_gdf['boolean'] = 0
port_raster_to_points_gdf['boolean'].loc[port_raster_to_points_gdf['port_val'] > 0] = 1


#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
port_grid_gdf = gpd.sjoin(grid, port_raster_to_points_gdf, how='inner', op='contains')
port_grid_norm_gdf = port_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean', 'grid_index'], axis=1)

#**************************#







# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#






#*** Shipping Routes ***#


bbox= (xmin_wgs, ymin_wgs, xmax_wgs, ymax_wgs)
ships_gdf = gpd.read_file(Shipping_Routes_shp, bbox=bbox)

ships_gdf_utm = ships_gdf.to_crs(crs_utm)
if len(ships_gdf_utm) > 0 :
    grid_ships_join = gpd.sjoin(grid, ships_gdf_utm, how="inner", op='intersects')
    grid_ships = grid_ships_join.drop(['grid_index', 'index_right', 'OBJECTID', 'Id', 'Routes', 'Shape_Leng'], axis=1)
    grid_ships['ship_val'] = 1
else:
    grid_ships = pd.DataFrame()

#**************************#






# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#














# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#






#*** Snap 8 - Raster to polygon ***#


snap8_raster = rasterio.open(OutFolder + "lu_snap8_temp.tif")
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
    
snap8_norm_post_proc_df = pd.DataFrame(list(zip(snap8_pr_x_list, snap8_pr_y_list, snap8_pr_val_list)), columns =['xcor', 'ycor', 'snap8_val']) 
snap8_raster_to_points_gdf = gpd.GeoDataFrame(snap8_norm_post_proc_df, geometry=gpd.points_from_xy(x=snap8_norm_post_proc_df.xcor, y=snap8_norm_post_proc_df.ycor), crs = crs_utm)
snap8_raster_to_points_gdf['boolean'] = 0
snap8_raster_to_points_gdf['boolean'].loc[snap8_raster_to_points_gdf['snap8_val'] > 0] = 1

#*** Raster values to polygons - Spatial join of points including raster values with grid ***#
snap8_grid_gdf = gpd.sjoin(grid, snap8_raster_to_points_gdf, how='inner', op='contains')
snap8_grid_norm_gdf = snap8_grid_gdf.drop(['index_right', 'xcor', 'ycor', 'boolean', 'grid_index'], axis=1)

#**************************#







# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#






# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#







#*** Ports Final Claculation - Polygon to raster with shipping routes ***#


if len(grid_ships) > 0:
    port_final = gpd.sjoin(port_grid_norm_gdf, grid_ships, how="left", op='contains')
    port_final['ship_val'].loc[port_final['ship_val'] != 1] = 0
    port_final['port_val'].loc[port_final['ship_val'] == 1] = 1
else:
    port_final = port_grid_norm_gdf.copy()


#*** Snap 8 Final Claculation - Polygon to raster with shipping routes ***#


if len(grid_ships) > 0:
    snap8_final = gpd.sjoin(snap8_grid_norm_gdf, grid_ships, how="left", op='contains')
    snap8_final['ship_val'].loc[snap8_final['ship_val'] != 1] = 0
    snap8_final['snap8_val'].loc[snap8_final['ship_val'] == 1] = 1
else:
    snap8_final = snap8_grid_norm_gdf.copy()


#*** Raster we are going to use as example for cell size, bounds and crs ***#

if len(grid_ships) > 0:
    rst_fn = OutFolder + "/lu_industry.tif"

    ## Output Rasters

    ##snap8_final_raster = OutFolder + "/lu_snap8_full.tif"
    snap8_final_raster = OutFolder + "/lu_snap8.tif"
    ship_raster = OutFolder + "/lu_shipping.tif"
    ship_with_ports_raster = OutFolder + "/lu_shipping_with ports.tif"

    ## Start process to burn geodataframe to geotiff

    rst = rasterio.open(rst_fn)

    meta = rst.meta.copy()
    meta.update(compress='lzw')



    with rasterio.open(snap8_final_raster, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(snap8_final.geometry, snap8_final.snap8_val))

        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)



    with rasterio.open(ship_raster, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(grid_ships.geometry, grid_ships.ship_val))

        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)



    with rasterio.open(ship_with_ports_raster, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(port_final.geometry, port_final.port_val))

        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)



else:

    rst_fn = OutFolder + "/lu_industry.tif"

    ## Output Rasters

    ##snap8_final_raster = OutFolder + "/lu_snap8_full.tif"
    snap8_final_raster = OutFolder + "/lu_snap8.tif"
    ship_with_ports_raster = OutFolder + "/lu_shipping_with ports.tif"

    ## Start process to burn geodataframe to geotiff

    rst = rasterio.open(rst_fn)

    meta = rst.meta.copy()
    meta.update(compress='lzw')



    with rasterio.open(snap8_final_raster, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(snap8_final.geometry, snap8_final.snap8_val))

        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)



    with rasterio.open(ship_with_ports_raster, 'w+', **meta) as out:
        out_arr = out.read(1)

        # this is where we create a generator of geom, value pairs to use in rasterizing
        shapes = ((geom,value) for geom, value in zip(port_final.geometry, port_final.port_val))

        burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
        out.write_band(1, burned)

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************







# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#





# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#









###**************************************************************************************************************
###*************************************** Reclassify Corine Polygons 2018 ***************************************#
###**************************************************************************************************************



print("Creating E-PRTR - Corine proxies...")




bbox= (xmin_corine - 100, ymin_corine - 100, xmax_corine + 100, ymax_corine + 100)

corine_polygons = gpd.read_file(Corine_GDB, driver='FileGDB', layer = Corine_GDB_layer[0], bbox=bbox)

corine_poly_croped = cl.clip_shp(corine_polygons, grid_corine)

corine_poly_croped.crs = crs_corine
corine_poly_wgs = corine_poly_croped.to_crs(crs_wgs)

corine_poly_final = corine_poly_wgs.drop(['Remark', 'PRTR_5', 'PRTR_34', 'PRTR_1', 'c18'], axis=1)

corine_poly_final.to_file(OutFolder + "/corine_poly_final.shp", driver='ESRI Shapefile')
corine_poly_final = gpd.read_file(OutFolder + "/corine_poly_final.shp")

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#





# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#









#**************************************************************************************************************
#*************************************** Read KMZ coordinates & save to df ***************************************#
#**************************************************************************************************************



##eprtr_df = geotable.load(InFolder + '/E_PRTR/kmz/eprtr_v17_athens.kmz')
eprtr_df = geotable.load(E_PRTR_kmz)




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




sectors = eprtr_df.geometry_layer.unique()
proj = eprtr_df.geometry_proj4.unique()

eprtr_df['Snap1'] = 0
eprtr_df['Snap6'] = 0
eprtr_df['Snap9'] = 0
eprtr_df['Snap34'] = 0
eprtr_df['Snap10'] = 0

eprtr_df["geometry_layer"]= eprtr_df["geometry_layer"].astype(str)

eprtr_df["Sector"]= eprtr_df["geometry_layer"].str.get(0)
eprtr_df["SubSector"]= eprtr_df["geometry_layer"].str.get(3)

eprtr_df["Sector"] = eprtr_df["Sector"].astype(str)
eprtr_df["SubSector"] = eprtr_df["SubSector"].astype(str)


# ~ *** ~ 1.(a) Mineral oil and gas refineries ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "1") & (eprtr_df['SubSector'] == "a"), 'Snap1'] = 1

# ~ *** ~ 1.(b) Gasification and liquefaction ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "1") & (eprtr_df['SubSector'] == "b"), 'Snap1'] = 1

# ~ *** ~ 1.(c) Thermal power stations and other combustion installations ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "1") & (eprtr_df['SubSector'] == "c"), 'Snap1'] = 1

# ~ *** ~ 1.(d) Coke ovens ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "1") & (eprtr_df['SubSector'] == "d"), 'Snap34'] = 1

# ~ *** ~ 1.(e) Coal rolling mills ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "1") & (eprtr_df['SubSector'] == "e"), 'Snap1'] = 1

# ~ *** ~ 1.(f) Manufacture of coal products and solid smokeless fuel ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "1") & (eprtr_df['SubSector'] == "f"), 'Snap1'] = 1

# ~ *** ~ 2.(a) Metal ore (including sulphide ore) roasting or sintering installations ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "2") & (eprtr_df['SubSector'] == "a"), 'Snap34'] = 1

# ~ *** ~ 2.(b) Installations for the production of pig iron or steel (primary or secondary ,elting) including continuous casting   ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "2") & (eprtr_df['SubSector'] == "b"), 'Snap34'] = 1

# ~ *** ~ 2.(c) Processing of ferrous metals ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "2") & (eprtr_df['SubSector'] == "c"), 'Snap34'] = 1

# ~ *** ~ 2.(d) Ferrous metal foundries ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "2") & (eprtr_df['SubSector'] == "d"), 'Snap34'] = 1

# ~ *** ~ 2.(e.1) Installations for the production of the non-ferrous crude metals from ore, concentrates or secondary raw materials by metallugical, chemical or electrolytic processes ~ *** ~ #
# ~ *** ~ 2.(e.2) Installation for the smelting, including the alloying, of the non-ferrous metals, including recovered products (refining, foundry casting,etc.) ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "2") & (eprtr_df['SubSector'] == "e"), 'Snap34'] = 1

# ~ *** ~ 2.(f) Surface treatment of metals and plastics using electrolytic or chemical processes ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "2") & (eprtr_df['SubSector'] == "f"), 'Snap34'] = 1

# ~ *** ~ 3.(a) Underground mining and related operations ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "3") & (eprtr_df['SubSector'] == "a"), 'Snap34'] = 1

# ~ *** ~ 3.(b) Opencast mining and quarrying ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "3") & (eprtr_df['SubSector'] == "b"), 'Snap34'] = 1

# ~ *** ~ 3.(c) Production of cement clinker or lime in rotary kilns or other furnaces ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "3") & (eprtr_df['SubSector'] == "c"), 'Snap34'] = 1

# ~ *** ~ 3.(d) Installations for the production of the asbestos and the manifracture of asbestos-based products ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "3") & (eprtr_df['SubSector'] == "d"), 'Snap34'] = 1

# ~ *** ~ 3.(e) Installations for the manifracture of the glass, including glass fibre ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "3") & (eprtr_df['SubSector'] == "e"), 'Snap34'] = 1

# ~ *** ~ 3.(f) Melting mineral substances, including the production of mineral fibres ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "3") & (eprtr_df['SubSector'] == "f"), 'Snap34'] = 1

# ~ *** ~ 3.(g) Manufacture of ceramic products including tiles, bricks, stoneware or porcelain ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "3") & (eprtr_df['SubSector'] == "g"), 'Snap34'] = 1

# ~ *** ~ 4.(a) Industrial scale production of basic organic chemicals ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "4") & (eprtr_df['SubSector'] == "a"), 'Snap34'] = 1

# ~ *** ~ 4.(b) Industrial scale production of basic inorganic chemicals ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "4") & (eprtr_df['SubSector'] == "b"), 'Snap34'] = 1

# ~ *** ~ 4.(c) Industrial scale production of phosphorous, nitrogen or potassium based fertilizers ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "4") & (eprtr_df['SubSector'] == "c"), 'Snap34'] = 1

# ~ *** ~ 4.(d) Industrial scale production of basic plant health products and of b ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "4") & (eprtr_df['SubSector'] == "d"), 'Snap34'] = 1

# ~ *** ~ 4.(e) Industrial scale production of basic pharmaceutical products ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "4") & (eprtr_df['SubSector'] == "e"), 'Snap34'] = 1

# ~ *** ~ 4.(f) Industrial scale production of explosives and pyrotechnic products ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "4") & (eprtr_df['SubSector'] == "f"), 'Snap34'] = 1

# ~ *** ~ 5.(a) Instalation and recovery or disposal of hazardous waste ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "5") & (eprtr_df['SubSector'] == "a"), 'Snap9'] = 1

# ~ *** ~ 5.(b) Incineration of non-hazardous waste included in Directive 2000/76/EC - waste incineration ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "5") & (eprtr_df['SubSector'] == "b"), 'Snap9'] = 1

# ~ *** ~ 5.(c) Disposal of non-hazardous waste ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "5") & (eprtr_df['SubSector'] == "c"), 'Snap9'] = 1

# ~ *** ~ 5.(d) Landfills (excluding landfills closed before the 16.7.2001) ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "5") & (eprtr_df['SubSector'] == "d"), 'Snap9'] = 1

# ~ *** ~ 5.(e) Disposal or recycling of animal carcasses and animal waste ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "5") & (eprtr_df['SubSector'] == "e"), 'Snap9'] = 1

# ~ *** ~ 5.(f) Urban waste-water treatment plants ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "5") & (eprtr_df['SubSector'] == "f"), 'Snap9'] = 1

# ~ *** ~ 5.(g) Independently operated industrial waste-water treatment plants serving a listed activity ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "5") & (eprtr_df['SubSector'] == "g"), 'Snap9'] = 1

# ~ *** ~ 6.(a) Production of pulp from timber or similar fibrous materials ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "6") & (eprtr_df['SubSector'] == "a"), 'Snap34'] = 1

# ~ *** ~ 6.(b) Production of paper and board and other primary wood products ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "6") & (eprtr_df['SubSector'] == "b"), 'Snap34'] = 1

# ~ *** ~ 6.(c) Preservation of wood and wood products with chemicals ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "6") & (eprtr_df['SubSector'] == "c"), 'Snap34'] = 1

# ~ *** ~ 7.(a) Intensive rearing of poultry or pigs ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "7") & (eprtr_df['SubSector'] == "a"), 'Snap10'] = 1

# ~ *** ~ 7.(b) Intensive aquaculture ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "7") & (eprtr_df['SubSector'] == "b"), 'Snap10'] = 1

# ~ *** ~ 8.(a) Slaughterhouses ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "8") & (eprtr_df['SubSector'] == "a"), 'Snap34'] = 1

# ~ *** ~ 8.(b) Treatment and processing of animal and vegetable materials in food and drink production ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "8") & (eprtr_df['SubSector'] == "b"), 'Snap34'] = 1

# ~ *** ~ 8.(c) Treatment and processing of milk ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "8") & (eprtr_df['SubSector'] == "c"), 'Snap34'] = 1

# ~ *** ~ 9.(a) Pretreatment or dyeing of fibres or textiles ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "9") & (eprtr_df['SubSector'] == "a"), 'Snap6'] = 1

# ~ *** ~ 9.(b) Tanning of hides and skins ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "9") & (eprtr_df['SubSector'] == "b"), 'Snap6'] = 1

# ~ *** ~ 9.(c) Surface treatment of substances, objects or products using organic solvents ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "9") & (eprtr_df['SubSector'] == "c"), 'Snap6'] = 1

# ~ *** ~ 9.(d) Production of carbon or electro-graphite through incineration or graphitization ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "9") & (eprtr_df['SubSector'] == "d"), 'Snap6'] = 1

# ~ *** ~ 9.(e) Building of, painting or removal of paint from ships ~ *** ~ #
eprtr_df.loc[(eprtr_df['Sector'] == "9") & (eprtr_df['SubSector'] == "e"), 'Snap6'] = 1




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#





# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




# ~ #**************** ~ Create Snap9 Proxy from E-PRTR ~ ****************#

print("Calculate Snap34...")

# ~ *** ~ Geotable select only snap9 / drop NA from geometry col / drop extra cols ~ *** #

eprtr_snap34_df = eprtr_df.loc[eprtr_df['Snap34'] == 1]
eprtr_snap34_df_no_nan = eprtr_snap34_df.dropna(subset=['geometry_object'])
eprtr_snap34_df_final = eprtr_snap34_df_no_nan.drop(['Name', 'geometry_proj4', 'Snap1', 'Snap6', 'Snap9', 'Snap10'], axis=1)

# ~ *** ~ From Geotable to Geodataframe via list and np.array ~ *** ~ #

point_list = eprtr_snap34_df_final.to_numpy().tolist()
point_snap34_df = DataFrame(point_list,columns=['geometry_object', 'geometry_layer', 'Snap34', 'Sector', 'SubSector'])
point_snap34_gdf = gpd.GeoDataFrame(point_snap34_df, geometry='geometry_object')

try:
    # ~ *** ~ Crop and project points with grid ~ *** ~ #

    eprtr_croped_snap34 = cl.clip_shp(point_snap34_gdf, grid_wgs)
    eprtr_croped_snap34.crs = crs_wgs

    # ~ *** ~ Intersect points with grid = boolean df T/F ~ *** ~ #

    grid_wgs_inter_point_34 = gpd.sjoin(corine_poly_final, eprtr_croped_snap34, how='inner', op='contains')
    ##grid_wgs_inter_point_34 = gpd.sjoin(grid_wgs, eprtr_croped_snap34, how='inner', op='intersects')

    grid_wgs_inter_point_34['index'] = grid_wgs_inter_point_34.index

    corine_poly_final['index'] = corine_poly_final.index
    corine_poly_final['Inter'] = 0

    grouped_inter = grid_wgs_inter_point_34.groupby(['index'], as_index=False).sum()
    grouped_inter['Inter'] = 1

    inter_point_merged_34 = pd.merge(corine_poly_final, grouped_inter, on = ["index","Inter"], how='outer')

    snap34_grouped = inter_point_merged_34.groupby(['index'], as_index=False).sum()

    corine_poly_final['snap34'] = snap34_grouped['Inter']
    grid_drop = corine_poly_final.drop(['index', 'Inter'], axis=1)
    snap34_df_proxy = grid_drop.to_crs(crs_utm)

except:
    print("There are no E-PRTR points in the selected domain")
    
    # ~ *** ~ Intersect points with grid = boolean df T/F ~ *** ~ #

    grid_wgs_inter_point_34 = corine_poly_final.copy()

    grid_wgs_inter_point_34['index'] = grid_wgs_inter_point_34.index

    corine_poly_final['index'] = corine_poly_final.index
    corine_poly_final['Inter'] = 0

    inter_point_merged_34 = corine_poly_final.copy()

    snap34_grouped = inter_point_merged_34.groupby(['index'], as_index=False).sum()

    corine_poly_final['snap34'] = snap34_grouped['Inter']
    grid_drop = corine_poly_final.drop(['index', 'Inter'], axis=1)
    snap34_df_proxy = grid_drop.to_crs(crs_utm)


#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#





# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#





# ~ #**************** ~ Create Snap1 Proxy from E-PRTR ~ ****************#

print("Calculate Snap1...")

# ~ *** ~ Geotable select only snap1 / drop NA from geometry col / drop extra cols ~ *** #

eprtr_snap1_df = eprtr_df.loc[eprtr_df['Snap1'] == 1]
eprtr_snap1_df_no_nan = eprtr_snap1_df.dropna(subset=['geometry_object'])
eprtr_snap1_df_final = eprtr_snap1_df_no_nan.drop(['Name', 'geometry_proj4', 'Snap34', 'Snap6', 'Snap9', 'Snap10'], axis=1)

# ~ *** ~ From Geotable to Geodataframe via list and np.array ~ *** ~ #

point_list = eprtr_snap1_df_final.to_numpy().tolist()
point_snap1_df = DataFrame(point_list,columns=['geometry_object', 'geometry_layer', 'snap1', 'Sector', 'SubSector'])
point_snap1_gdf = gpd.GeoDataFrame(point_snap1_df, geometry='geometry_object')

try:
    # ~ *** ~ Crop and project points with grid ~ *** ~ #

    eprtr_croped_snap1 = cl.clip_shp(point_snap1_gdf, grid_wgs)
    eprtr_croped_snap1.crs = crs_wgs

    # ~ *** ~ Intersect points with grid = boolean df T/F ~ *** ~ #

    grid_wgs_inter_point_1 = gpd.sjoin(corine_poly_final, eprtr_croped_snap1, how='inner', op='contains')
    ##grid_wgs_inter_point_34 = gpd.sjoin(grid_wgs, eprtr_croped_snap1, how='inner', op='intersects')

    grid_wgs_inter_point_1['index'] = grid_wgs_inter_point_1.index

    corine_poly_final['index'] = corine_poly_final.index
    corine_poly_final['Inter'] = 0

    grouped_inter = grid_wgs_inter_point_1.groupby(['index'], as_index=False).sum()
    grouped_inter['Inter'] = 1

    inter_point_merged_1 = pd.merge(corine_poly_final, grouped_inter, on = ["index","Inter"], how='outer')

    snap1_grouped = inter_point_merged_1.groupby(['index'], as_index=False).sum()

    corine_poly_final['snap1'] = snap1_grouped['Inter']
    grid_drop = corine_poly_final.drop(['index', 'Inter'], axis=1)
    snap1_df_proxy = grid_drop.to_crs(crs_utm)

except:

    print("There are no E-PRTR points in the selected domain")
    
    # ~ *** ~ Intersect points with grid = boolean df T/F ~ *** ~ #

    grid_wgs_inter_point_1 = corine_poly_final.copy()

    grid_wgs_inter_point_1['index'] = grid_wgs_inter_point_1.index

    corine_poly_final['index'] = corine_poly_final.index
    corine_poly_final['Inter'] = 0

    inter_point_merged_1 = corine_poly_final.copy()

    snap1_grouped = inter_point_merged_1.groupby(['index'], as_index=False).sum()

    corine_poly_final['snap1'] = snap1_grouped['Inter']
    grid_drop = corine_poly_final.drop(['index', 'Inter'], axis=1)
    snap1_df_proxy = grid_drop.to_crs(crs_utm)


    
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#





# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#









#**************************************************************************************************************
#*************************************** Raster to Points Geodataframe ***************************************#
#**************************************************************************************************************



#*** Create Snap9 Proxy from E-PRTR ***#


print("Calculate Snap9...")

# ~ *** ~ Geotable select only snap9 / drop NA from geometry col / drop extra cols ~ *** #

eprtr_snap9_df = eprtr_df.loc[eprtr_df['Snap9'] == 1]
eprtr_snap9_df_no_nan = eprtr_snap9_df.dropna(subset=['geometry_object'])
eprtr_snap9_df_final = eprtr_snap9_df_no_nan.drop(['Name', 'geometry_proj4', 'Snap1', 'Snap6', 'Snap34', 'Snap10'], axis=1)

# ~ *** ~ From Geotable to Geodataframe via list and np.array ~ *** ~ #

point_list_snap9 = eprtr_snap9_df_final.to_numpy().tolist()
point_snap9_df = DataFrame(point_list_snap9,columns=['geometry_object', 'geometry_layer', 'Snap9', 'Sector', 'SubSector'])
point_snap9_gdf = gpd.GeoDataFrame(point_snap9_df, geometry='geometry_object')
point_snap9_gdf.crs = crs_wgs

try:
    # ~ *** ~ Crop and project points with grid ~ *** ~ #

    eprtr_croped_snap9 = cl.clip_shp(point_snap9_gdf, grid_wgs)
    eprtr_croped_snap9.crs = crs_wgs
    eprtr_croped_snap9.to_file(OutFolder + "eprtr_croped_snap9.shp", driver="ESRI Shapefile")



    # ~ #**************** ~ Clear Memory ~ ****************#
    gc.collect()
    # ~ #**************************************************#



    # ~ *** ~ Intersect points with grid = boolean df T/F ~ *** ~ #

    grid_wgs_inter_point_9 = gpd.sjoin(corine_poly_final, eprtr_croped_snap9, how='inner', op='contains')

    grid_wgs_inter_point_9['index'] = grid_wgs_inter_point_9.index

    corine_poly_final['index'] = corine_poly_final.index
    corine_poly_final['Inter'] = 0

    grouped_inter = grid_wgs_inter_point_9.groupby(['index'], as_index=False).sum()
    grouped_inter['Inter'] = 1

    inter_point_merged_9 = pd.merge(corine_poly_final, grouped_inter, on = ["index","Inter"], how='outer')

    snap9_grouped = inter_point_merged_9.groupby(['index'], as_index=False).sum()

    corine_poly_final['snap9'] = snap9_grouped['Inter']
    grid_drop = corine_poly_final.drop(['index', 'Inter'], axis=1)
    snap9_df_proxy = grid_drop.to_crs(crs_utm)

except:
    print("There are no E-PRTR points in the selected domain")
    
    # ~ *** ~ Intersect points with grid = boolean df T/F ~ *** ~ #

    grid_wgs_inter_point_9 = corine_poly_final.copy()

    grid_wgs_inter_point_9['index'] = grid_wgs_inter_point_9.index

    corine_poly_final['index'] = corine_poly_final.index
    corine_poly_final['Inter'] = 0

    inter_point_merged_9 = corine_poly_final.copy()

    snap9_grouped = inter_point_merged_9.groupby(['index'], as_index=False).sum()

    corine_poly_final['snap9'] = snap9_grouped['Inter']
    grid_drop = corine_poly_final.drop(['index', 'Inter'], axis=1)
    snap9_df_proxy = grid_drop.to_crs(crs_utm)


    
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#




# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#





# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#









#*** Snap 9 FULL ***#


snap9_df_proxy['snap9'].loc[(snap9_df_proxy['Code_18'] == '211')|(snap9_df_proxy['Code_18'] == '212')|(snap9_df_proxy['Code_18'] == '213')|(snap9_df_proxy['Code_18'] == '221')|(snap9_df_proxy['Code_18'] == '222')|(snap9_df_proxy['Code_18'] == '223')|(snap9_df_proxy['Code_18'] == '231')|(snap9_df_proxy['Code_18'] == '241')|(snap9_df_proxy['Code_18'] == '242')|(snap9_df_proxy['Code_18'] == '243')|(snap9_df_proxy['Code_18'] == '244')] = 1


#*** Snap 34 FULL ***#


snap9_df_proxy['snap34'].loc[(snap9_df_proxy['Code_18'] == '121')|(snap9_df_proxy['Code_18'] == '131')|(snap9_df_proxy['Code_18'] == '132')|(snap9_df_proxy['Code_18'] == '133')] = 1


#*** Snap 1 FULL ***#


snap9_df_proxy['snap1'].loc[(snap9_df_proxy['Code_18'] == '121')|(snap9_df_proxy['Code_18'] == '131')|(snap9_df_proxy['Code_18'] == '132')|(snap9_df_proxy['Code_18'] == '133')] = 1


snap9_df_proxy.to_file(OutFolder + "eprtr_df_proxy.shp", driver="ESRI Shapefile")



# ~ Convert shp to geotiff - Create snap9 & snap34 raster proxies

## Read geodataframe - shp from folder

shp_fn = OutFolder + "eprtr_df_proxy.shp"

## Raster we are going to use as example for cell size, bounds and crs

rst_fn = OutFolder + "/lu_industry.tif"

## Output Rasters

waste_raster = OutFolder + "/lu_waste.tif"
snap34_raster = OutFolder + "/lu_snap34.tif"
snap1_raster = OutFolder + "/lu_snap1.tif"

## Start process to burn geodataframe to geotiff

eprtr_df_proxy_gdf = gpd.read_file(shp_fn)

rst = rasterio.open(rst_fn)

meta = rst.meta.copy()
meta.update(compress='lzw')


#*** Write Snap 9 Raster ***#


with rasterio.open(waste_raster, 'w+', **meta) as out:
    out_arr = out.read(1)

    # this is where we create a generator of geom, value pairs to use in rasterizing
    shapes = ((geom,value) for geom, value in zip(eprtr_df_proxy_gdf.geometry, eprtr_df_proxy_gdf.snap9))

    burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
    out.write_band(1, burned)


#*** Write Snap 34 Raster ***#


with rasterio.open(snap34_raster, 'w+', **meta) as out:
    out_arr = out.read(1)

    # this is where we create a generator of geom, value pairs to use in rasterizing
    shapes = ((geom,value) for geom, value in zip(eprtr_df_proxy_gdf.geometry, eprtr_df_proxy_gdf.snap34))

    burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
    out.write_band(1, burned)


#*** Write Snap 1 Raster ***#


with rasterio.open(snap1_raster, 'w+', **meta) as out:
    out_arr = out.read(1)

    # this is where we create a generator of geom, value pairs to use in rasterizing
    shapes = ((geom,value) for geom, value in zip(eprtr_df_proxy_gdf.geometry, eprtr_df_proxy_gdf.snap1))

    burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
    out.write_band(1, burned)
    
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#








#**************************************************************************************************************
#*************************************** Download OSM Data *****************************************************#
#**************************************************************************************************************



print("Creating OSM proxy...")




typeOfRoad_len = []

bbox=[ymin_wgs, xmin_wgs, ymax_wgs, xmax_wgs]


nominatim = Nominatim()
overpass = Overpass()


query_motor = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="motorway"', includeGeometry=True)
query_motor_link = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="motorway_link"', includeGeometry=True)
query_primary = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="primary"', includeGeometry=True)
query_primary_link = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="primary_link"', includeGeometry=True)
query_secondary = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="secondary"', includeGeometry=True)
query_secondary_link = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="secondary_link"', includeGeometry=True)
query_trunk = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="trunk"', includeGeometry=True)
query_trunk_link = overpassQueryBuilder(bbox=bbox, elementType='way', selector='"highway"="trunk_link"', includeGeometry=True)





filename = OutFolder + '\\osm_roads_domain_wgs.shp'
w = shapefile.Writer(filename, shapeType = 3)
w.field('ID', 'N')
w.field('key', 'C')
w.field('typeOfRoad', 'C')
fid = 0








#*** Motorway ***#


result_motor = overpass.query(query_motor, timeout=90000)
elements_motor = result_motor.elements()
firstElement_motor = result_motor.elements()[0]
geometry_motor = firstElement_motor.geometry()
print(len(elements_motor))

for k in range(0,len(elements_motor)):
    geometry_motor = elements_motor[k].geometry()
    coordinates_motor = geometry_motor.get("coordinates")
    final_coordinates = []
    final_coordinates.append(coordinates_motor)
    w.line(final_coordinates)
    w.record(fid, 'Highway', 'motorway')
    fid = fid + 1








#*** Motorway link ***#


result_motor_link = overpass.query(query_motor_link, timeout=90000)
elements_motor_link = result_motor_link.elements()
firstElement_motor_link = result_motor_link.elements()[0]
geometry_motor_link = firstElement_motor_link.geometry()
print(len(elements_motor_link))

for k in range(0,len(elements_motor_link)):
    geometry_motor_link = elements_motor_link[k].geometry()
    coordinates_motor_link = geometry_motor_link.get("coordinates")
    final_coordinates = []
    final_coordinates.append(coordinates_motor_link)
    w.line(final_coordinates)
    w.record(fid, 'Highway', 'motorway_link')
    fid = fid + 1








#*** Primary ***#


result_primary = overpass.query(query_primary, timeout=90000)
elements_primary = result_primary.elements()
firstElement_primary = result_primary.elements()[0]
geometry_primary = firstElement_primary.geometry()
print(len(elements_primary))

for k in range(0,len(elements_primary)):
##for k in range(0,2):
    geometry_primary = elements_primary[k].geometry()
    coordinates_primary = geometry_primary.get("coordinates")
    if len(coordinates_primary) > 1:
        final_coordinates = []
        final_coordinates.append(coordinates_primary)
    else:
        final_coordinates = coordinates_primary
        

    w.line(final_coordinates)
    w.record(fid, 'Highway', 'primary')
    fid = fid + 1







#*** Primary link ***#


result_primary_link = overpass.query(query_primary_link, timeout=90000)
elements_primary_link = result_primary_link.elements()
firstElement_primary_link = result_primary_link.elements()[0]
geometry_primary_link = firstElement_primary_link.geometry()
print(len(elements_primary_link))

for k in range(0,len(elements_primary_link)):
    geometry_primary_link = elements_primary_link[k].geometry()
    coordinates_primary_link = geometry_primary_link.get("coordinates")
    final_coordinates = []
    final_coordinates.append(coordinates_primary_link)
    w.line(final_coordinates)
    w.record(fid, 'Highway', 'primary_link')
    fid = fid + 1








#*** Secondary ***#


result_secondary = overpass.query(query_secondary)
elements_secondary = result_secondary.elements()
firstElement_secondary = result_secondary.elements()[0]
geometry_secondary = firstElement_secondary.geometry()
print(len(elements_secondary))

for k in range(0,len(elements_secondary)):
    geometry_secondary = elements_secondary[k].geometry()
    coordinates_secondary = geometry_secondary.get("coordinates")
    if len(coordinates_secondary) > 1:
        final_coordinates = []
        final_coordinates.append(coordinates_secondary)
    else:
        final_coordinates = coordinates_secondary
    w.line(final_coordinates)
    w.record(fid, 'Highway', 'secondary')
    fid = fid + 1








#*** Secondary link ***#


result_secondary_link = overpass.query(query_secondary_link, timeout=90000)
elements_secondary_link = result_secondary_link.elements()
firstElement_secondary_link = result_secondary_link.elements()[0]
geometry_secondary_link = firstElement_secondary_link.geometry()
print(len(elements_secondary_link))

for k in range(0,len(elements_secondary_link)):
    geometry_secondary_link = elements_secondary_link[k].geometry()
    coordinates_secondary_link = geometry_secondary_link.get("coordinates")
    final_coordinates = []
    final_coordinates.append(coordinates_secondary_link)
    w.line(final_coordinates)
    w.record(fid, 'Highway', 'secondary_link')
    fid = fid + 1








#*** Trunk ***#


result_trunk = overpass.query(query_trunk, timeout=90000)
elements_trunk = result_trunk.elements()
firstElement_trunk = result_trunk.elements()[0]
geometry_trunk = firstElement_trunk.geometry()
print(len(elements_trunk))

for k in range(0,len(elements_trunk)):
    geometry_trunk = elements_trunk[k].geometry()
    coordinates_trunk = geometry_trunk.get("coordinates")
    final_coordinates = []
    final_coordinates.append(coordinates_trunk)
    w.line(final_coordinates)
    w.record(fid, 'Highway', 'trunk')
    fid = fid + 1







#*** Trunk link ***#


result_trunk_link = overpass.query(query_trunk_link, timeout=90000)
elements_trunk_link = result_trunk_link.elements()
firstElement_trunk_link = result_trunk_link.elements()[0]
geometry_trunk_link = firstElement_trunk_link.geometry()
print(len(elements_trunk_link))

for k in range(0,len(elements_trunk_link)):
    geometry_trunk_link = elements_trunk_link[k].geometry()
    coordinates_trunk_link = geometry_trunk_link.get("coordinates")
    final_coordinates = []
    final_coordinates.append(coordinates_trunk_link)
    w.line(final_coordinates)
    w.record(fid, 'Highway', 'trunk_link')
    fid = fid + 1




# create the PRJ file
prj = open(OutFolder + '\\osm_roads_domain_wgs.prj', "w")
epsg = 'GEOGCS["GCS_WGS_1984",'
epsg += 'DATUM["D_WGS_1984",'
epsg += 'SPHEROID["WGS 84",6378137,298.257223563]]'
epsg += ',PRIMEM["Greenwich",0],'
epsg += 'UNIT["degree",0.0174532925199433]]'
prj.write(epsg)
prj.close()


w.close()

osm_wgs = gpd.read_file(filename)
osm_utm34 = osm_wgs.to_crs(crs_utm)
osm_utm34.to_file(OutFolder + "osm_utm34.shp", driver="ESRI Shapefile")
for fname in os.listdir(OutFolder):
    if fname.startswith("osm_roads_domain_wgs"):
        os.remove(os.path.join(OutFolder, fname))

#********************************************************************************#
#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#








#**************************************************************************************************************
#*************************************** Create uni proxy *******************************************************#
#**************************************************************************************************************



print("Creating uni proxy...")





# ~ Convert shp to geotiff - Create snap9 & snap34 raster proxies

## Read geodataframe - shp from folder
grid['uni']=1

## Raster we are going to use as example for cell size, bounds and crs

rst_fn = OutFolder + "/lu_industry.tif"

## Output Rasters

uni_raster = OutFolder + "/lu_uni.tif"

## Start process to burn geodataframe to geotiff

rst = rasterio.open(rst_fn)

meta = rst.meta.copy()
meta.update(compress='lzw')

with rasterio.open(uni_raster, 'w+', **meta) as out:
    out_arr = out.read(1)

    # this is where we create a generator of geom, value pairs to use in rasterizing
    shapes = ((geom,value) for geom, value in zip(grid.geometry, grid.uni))

    burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
    out.write_band(1, burned)

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








# ~ #**************** ~ Clear Memory ~ ****************#
gc.collect()
# ~ #**************************************************#








###**************************************************************************************************************
###*************************************** GHS UC proxy ************************************************************#
###**************************************************************************************************************



print("Creating ghs proxy...")




ghs = gpd.read_file(GHS_UC_shp)

ghs_cliped = gpd.clip(ghs, grid_wgs)

ghs_drop = ghs_cliped.filter(['CTR_MN_ISO', 'UC_NM_MN', 'geometry'])

ghs_utm = ghs_drop.to_crs(crs_utm)
ghs_utm['boolean'] = uc_increase_factor


# ~ Convert shp to geotiff - Create snap9 & snap34 raster proxies

## Raster we are going to use as example for cell size, bounds and crs

rst_fn = OutFolder + "/lu_industry.tif"

## Output Rasters

ghs_raster = OutFolder + "/ghs_incr_fact_" + str(uc_increase_factor) + ".tif"

## Start process to burn geodataframe to geotiff

rst = rasterio.open(rst_fn)

meta = rst.meta.copy()
meta.update(compress='lzw')

with rasterio.open(ghs_raster, 'w+', **meta) as out:
    out_arr = out.read(1)

    # this is where we create a generator of geom, value pairs to use in rasterizing
    shapes = ((geom,value) for geom, value in zip(ghs_utm.geometry, ghs_utm.boolean))

    burned = features.rasterize(shapes=shapes, fill=0, out=out_arr, transform=out.transform)
    out.write_band(1, burned)

#**************************************************************************************************************
#**************************************************************************************************************
#**************************************************************************************************************








elapsed_time = time.time() - start_time
hours, rem = divmod(elapsed_time, 3600)
minutes, seconds = divmod(rem, 60)
print("minutes", minutes)
print("seconds", elapsed_time)

##print("--- %s seconds ---" % (time.time() - start_time))



