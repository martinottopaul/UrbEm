# UrbEm v1.0
The Urban Emission downscaling model (UrbEm)for air quality modeling

0. Introduction

The use of regional emission inventories can be challenging for urban-scale AQ applications and air quality management in cities. Nevertheless, their exploitation through disaggregation by utilizing spatial proxies is a credible solution for European cities that lack bottom-up emission inventories. 

To this end, we developed the UrbEm approach, which enables in a modular manner downscaling of gridded regional emissions with specific spatial proxies based on a variety of open access, robust, sustainable and frequently updated sources. UrbEm can be applied to any urban area in Europe and provides methodological homogeneity between different cities.

To demonstrate the general applicability and performance of the developed method and tool, we introduced the method, and compared the spatial distribution of uniformly disaggregated regional emissions with emissions downscaled with the UrbEm approach for the differing cities of Athens and Hamburg (manuscript submitted for publication, pre-print accessible on request). 

The UrbEm downscaling approach is completely free of cost and open source. Its application is realized in (1) a series of R scripts and (2) as Pyhton script. Both applications rely on e.g. CAMS-REG emission inventories as well as a set (maps) of spatial proxies, which need to be downloaded before using UrbEm.

1. Access necessary input data

1.1. Emission datasets

UrbEm v1.0 is configured to read CAMS-REG-AP v3.1 regional emissions. After registration, these can be downloaded free of cost at: https://eccad.aeris-data.fr
After registration in the "Access Data" section, the CAMS-REG-AP dataset should be selected and downloaded for all pollutants and sectors.

The second emission database applied is the European Pollutant Release and Transfer Register E-PRTR, which can be downloaded without registration at: https://www.eea.europa.eu/data-and-maps/data/member-states-reporting-art-7-under-the-european-pollutant-release-and-transfer-register-e-prtr-regulation-23/european-pollutant-release-and-transfer-register-e-prtr-data-base

Both datasets should be placed in separate folders.

1.2. Spatial proxies
Besides a collection of spatial proxies (download here: https://doi.org/10.5281/zenodo.5508739), which have been specifically prepared for application in UrbEm, the Global Human Settlement Layer (download here: https://ghsl.jrc.ec.europa.eu/ghs_pop2019.php) Population density product "GHS_POP_E2015_GLOBE_R2019A_4326_30SS_V1_0" needs to be downloaded.

Make sure all proxies are placed in one folder.


2. Apply UrbEm v1.0

Both solutions (R and Python) are configured to 
(1) read CAMS-REG-AP v3.1 and E-PRTR emission input files, 
(2) downscale these gridded and point emissions with the downloaded set of spatial proxies,
(3) to arrive at annual total emissions for a selected year and a selected urban domain,
(4) and write these as area, line and/or point source emission files,
(5) in a *.csv file format for the EPISODE-CityChem preprocessor UECT (Karl et al. 2019).

Although UrbEm v1.0 delivers only UECT/EPISODE-CityChem file format as output, it is generally possible to change to the desired output format by code modification. Nevertheless, we promote to use the EPISODE-CityChem air quality model for urban-scales due to its efficiency, performance and ongoing development. EPISODE-CityChem can be downloaded free of cost at https://doi.org/10.5281/zenodo.1116173.

2.1. UrbEm Rscripts

The UrbEm v1.0 Rscripts are separated in three main scripts:
1_UrbEm_pointsources_v1.R
2_UrbEm_areasources_v2.R
3_UrbEm_linesources_v3.

These scripts need to be run sequentially to create point, area and line emission files. 

Before running the scripts, make sure the R libraries raster, sf, rgdal, osmdata and lwgeom are installed in your R environment. 

Addtionally the following auxiliary functions (scripts distributed with UrbEm v1.0) are necessary and will be sourced at the beginning of some main scripts:
areasources_e-prtr_pointsource_correction.R
areasources_to_osm_linesources.R
proxy_distribution.R
proxy_preparation.R

While there are no changes in the auxiliary scripts necessary to run UrbEm, there need to be changes made in the main scripts. Each of the main scripts has an input section at the beginning, which needs to be adjusted, for e.g.:
- setting input folders of emission files and proxies
- setting output folders and output text strings
- setting a domain definition
- setting downscaling options

The input section of each main script, as well as the code itself delivers a documentation of each step made in the script.

Feel free to adjust the code for your purposes.



2.2. UrbEm Python scripts

The UrbEm v1.0 Python scripts are separated in two main scripts: 
- 1_UrbEM_proxies_v1.py 
- 2_UrbEM_emissions_v1.py

These scripts need to be run sequentially to 1. create proxy, point, area and line emission files.

Before running the scripts, make sure the python libraries pandas, numpy, gdal, geopandas, os, sys, fnmatch, inspect, rasterio, rasterio.mask, earthpy.spatial, shapely.geometry, earthpy, fiona, osgeo, gc, geotable, pyproj, shapely, time, shutil, OSMPythonTools.nominatim, OSMPythonTools.overpass, OSMPythonTools.data, collections and shapefile are installed in your Python v3 environment.

Spatial datasets: In order to be able to run the python scripts the user should also download:
- Population density data (Global dataset/ 2015 / WGS84 / 30 arcsec): https://ghsl.jrc.ec.europa.eu/download.php?ds=pop 
- CORINE raster and GDB files: https://land.copernicus.eu/pan-european/corine-land-cover/clc2018 
- E-PRTR kmz data: 
https://www.eea.europa.eu/data-and-maps/data/member-states-reporting-art-7-under-the-european-pollutant-release-and-transfer-register-e-prtr-regulation-23/e-prtr-facilities-kmz-format/eprtr_facilities_v9.kmz 
- E-PRTR csv data: 
https://www.eea.europa.eu/data-and-maps/data/member-states-reporting-art-7-under-the-european-pollutant-release-and-transfer-register-e-prtr-regulation-23/european-pollutant-release-and-transfer-register-e-prtr-data-base/eprtr_v9_csv.zip 
- Urban Center data: https://ghsl.jrc.ec.europa.eu/ghs_stat_ucdb2015mt_r2019a.php 
- Eurostat countries: https://ec.europa.eu/eurostat/web/gisco/geodata/reference-data/administrative-units-statistical-units/countries 
- Shipping Routes

Each of the main scripts has an input section at the beginning, which needs to be adjusted, for e.g.:
- setting input folders of emission files and proxies
- setting a domain definition
- setting downscaling options

The input section of each script, as well as the code itself delivers a documentation of each step made in the script.

Feel free to adjust the code for your purposes.
