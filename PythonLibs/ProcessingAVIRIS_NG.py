# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 15:25:01 2019

@author: wayne
"""


import json
import math
import affine
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from osgeo import gdal,ogr,osr


'''
Print ENVI header metadata; Open ENVI image with GDAL
'''

# open the ENVI file
img = gdal.Open("ang20180814t224053_rfl_v2r2/ang20180814t224053_corr_v2r2_img")

nbands = img.RasterCount
nrows = img.RasterYSize
ncols = img.RasterXSize

print("\n".join(["Bands:\t"+str(nbands),"Rows:\t"+str(nrows),"Cols:\t"+str(ncols)]))

'''
Get band information 
The header file also contains a list of bands and wavelengths. 
You can access a dictionary containing the wavelength of the center of each band using img.GetMetadata(). 
Make a table describing each of the bands:
'''

# band descriptions. we reference this dictionary throughout the tutorial.
band_dictionary = {
    "visible-violet": {'lower': 375, 'upper': 450, 'color': 'violet'},
    "visible-blue": {'lower': 450, 'upper': 485, 'color': 'blue'},
    "visible-cyan": {'lower': 485, 'upper': 500, 'color': 'cyan'},
    "visible-green": {'lower': 500, 'upper': 565, 'color': 'green'},
    "visible-yellow": {'lower': 565, 'upper': 590, 'color': 'yellow'},
    "visible-orange": {'lower': 590, 'upper': 625, 'color': 'orange'},
    "visible-red": {'lower': 625, 'upper': 740, 'color': 'red'},
    "near-infrared": {'lower': 740, 'upper': 1100, 'color': 'gray'},
    "shortwave-infrared": {'lower': 1100, 'upper': 2500, 'color': 'white'}
}

# function to classify bands
between = lambda wavelength, region: region['lower'] < wavelength <= region['upper']
def classifier(band):
    for region, limits in band_dictionary.items():
        if between(band, limits):
            return(region)

# lists of band numbers, band centers, and em classes
band_numbers = [int(b.split("_")[1]) for b in img.GetMetadata().keys() if b != "wavelength_units"]
band_centers = [float(b.split(" ")[0]) for b in img.GetMetadata().values() if b != "Nanometers"]
em_regions = [classifier(b) for b in band_centers]

# data frame describing bands
bands = pd.DataFrame({ 
    "Band number": band_numbers, 
    "Band center (nm)": band_centers, 
    "EM region": em_regions }, index = band_numbers).sort_index()

# print the first ten rows
bands.head(10)


'''
2. Extract image data at specific locations
Get geographic and image (pixel,line) coordinates for NGEE Arctic permafrost monitoring sites
'''
sites = "ngee-data/Teller_Permafrost_Monitoring/TL_Permafrost_Monitoring.shp"

# open with ogr
driver = ogr.GetDriverByName("ESRI Shapefile")
permafrost_sites = driver.Open(sites, 0)

# get the first feature in the shapefile as JSON
site1 = json.loads(permafrost_sites[0].GetFeature(0).ExportToJson())
site1


# shapefiles have a nested structure: layer(s) -> feature(s) -> geometry
lyr = permafrost_sites.GetLayer() # get the only layer in the shapefile
feat = lyr.GetFeature(1)          # get the first feature in the layer (1 feature per site)
geom = feat.GetGeometryRef()      # get the feature's geometry

# get transform for decimal degrees
from_srs = lyr.GetSpatialRef()                                         # get shapefile srs def
to_srs = osr.SpatialReference()                                        # init ogr srs object
to_srs.ImportFromEPSG(4326)                                            # import wgs84 srs def
xytransform = osr.CoordinateTransformation(from_srs,to_srs)            # get transform object

# get UTM and lat/long coordinates for each of the sites
utm_coordinate_pairs = {}
ll_coordinate_pairs = {}
for feature in lyr:
    geom = feature.GetGeometryRef()                                    # get site geometry
    utm_coordinate_pairs[feature['Name']] = (geom.GetX(), geom.GetY()) # get x,y utm coordinates 
    geom.Transform(xytransform)                                        # to wgs84
    ll_coordinate_pairs[feature['Name']] = (geom.GetX(), geom.GetY())  # get lon, lat
    
# get the x and y UTM coordinates for the first site
x, y = utm_coordinate_pairs['TL_IS_2']
affine_transform = affine.Affine.from_gdal(*img.GetGeoTransform())     # affine forward transform
inverse_transform = ~affine_transform                                  # invert transform
px, py = inverse_transform * (x, y)                                    # apply to x,y coordinates
px, py = int(px + 0.5), int(py + 0.5)                                  # get new x,y as integers

# print the three coordinates (UTM, geographic, image)
print( "\n".join(["Site 1 UTM coordinates (x,y): "+"\t"*4+str((x,y)),
       " are equal to geographic coordinates (lng,lat): \t"+str(ll_coordinate_pairs['TL_IS_2']),
       " and fall within image coordinates (pixel,line):\t"+str((px,py))]) )


'''
Once we have the inverse transform we can use it repeatedly for this raster dataset. 
Print the image coordinates for each of the sites:
'''
# get image coordinates for each site
image_coordinate_pairs = {
    site: inverse_transform * pair for site,pair in utm_coordinate_pairs.items()}

print("site id: (col, row)")
for site,coord in image_coordinate_pairs.items():  # must convert to int to raster array  
    print(site + ": (" + str(int(coord[0] + 0.5)) + ", " + str(int(coord[1] + 0.5)) +")")
    
    
'''
Plot the spectral curves for the first two permafrost monitoring sites
'''    
#site 1
site1name = list(image_coordinate_pairs.keys())[0]    # get name of site1
site1xy = list(image_coordinate_pairs.values())[0]    # get image coordinates of site1
px, py = int(site1xy[0] + 0.5), int(site1xy[1] + 0.5) # convert image coordinates to integers

band1_array = img.GetRasterBand(1).ReadAsArray()
print("Band 1 reflectance at site 1: "+str(band1_array[py,px]))


# function gets value at input xy from input band
get_pixel = lambda img,band,y,x: img.GetRasterBand(band).ReadAsArray()[y,x]

# make a copy of the bands data frame and add reflectance column for site 1
_bands = bands
_bands[site1name+" reflectance"] = [get_pixel(img,b,py,px) for b in range(1,nbands+1)]
_bands.head(10)


##################
#site 2
site2name = list(image_coordinate_pairs.keys())[1]    # get name of site1
site2xy = list(image_coordinate_pairs.values())[1]    # get image coordinates of site1
px, py = int(site2xy[0] + 0.5), int(site2xy[1] + 0.5) # convert image coordinates to int

_bands[site2name+" reflectance"] = [get_pixel(img,b,py,px) for b in range(1,nbands+1)]
_bands.head(10)




















