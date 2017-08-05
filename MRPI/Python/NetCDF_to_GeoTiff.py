#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 12:07:10 2017

@author: svimal
"""


from osgeo import gdal, osr
import xarray as xr
import numpy as np
import numpy

#def create_geotiff(netcdf_filename, var_name, soil=False,):

# Inputs
netcdf_filename = "/home2/svimal/Data/VIC_Fluxes/1920s/fluxes.1924-11.nc"
var_name = "Soil_liquid" #"Tair" #
var_name =  "Tair" #

xr_ensemble = xr.open_dataset(netcdf_filename)
data = xr_ensemble[var_name]
nodatavalue = -9999
data = np.ma.masked_array(data, mask=data==nodatavalue, fill_value=nodatavalue)
gdal_datatype = gdal.GDT_UInt16
np_datatype = numpy.uint16
driver = gdal.GetDriverByName( "GTiff" )

NCOLS = 310
NROWS = 350
XLLCORNER = -124.875
YLLCORNER = 31.125
CELLSIZE = 0.0625
NODATA_VALUE = -9999

nbands=365
dst_filename = "/home/svimal/test8.tif"
dst_ds = driver.Create(dst_filename, NCOLS, NROWS, nbands, gdal_datatype )

#Cellsize in same units as spatial reference
cellsize=0.0625

dst_ds.SetGeoTransform( [ XLLCORNER, CELLSIZE, 0, YLLCORNER, 0, CELLSIZE ] )
srs = osr.SpatialReference()
srs.SetWellKnownGeogCS("WGS84")
dst_ds.SetProjection( srs.ExportToWkt() )
for band, t in enumerate(xr_ensemble.time):
    #flipped_data = np.flipud(data[band].data)
    #dst_ds.GetRasterBand(band+1).WriteArray(flipped_data)
    dst_ds.GetRasterBand(band+1).WriteArray(data[band].data)

# Once we're done, close properly the dataset
dst_ds = None

#min & max random values of the output raster
#zmin=0; zmax=12345
#for band in range(nbands):
#    dst_ds.GetRasterBand(band+1).WriteArray( raster[band, :, :] )
#raster = numpy.random.randint(zmin,zmax, (nbands, nrows, ncols)).astype(np_datatype )
## See http://gdal.org/python/osgeo.gdal_array-module.html#codes
## for mapping between gdal and numpy data types
## These are only required if you wish to georeference (http://en.wikipedia.org/wiki/Georeference)
## your output geotiff, you need to know what values to input, don't just use the ones below
#Coordinates of the lower left corner of the image
#in same units as spatial reference
#xllcorner=147.2
#yllcorner=-34.54

#Cellsize in same units as spatial reference
#cellsize=0.01

#dst_ds.SetGeoTransform( [ xllcorner, cellsize, 0, yllcorner, 0, -cellsize ] )
#srs = osr.SpatialReference()
#srs.SetWellKnownGeogCS("WGS84")
#dst_ds.SetProjection( srs.ExportToWkt() )

###########################



#
#
#    if var_name is not "Soil_liquid":
#        for i, t in enumerate(data.time):
#            print i, t
#            flipped_data = np.flipud(data[i].data)
#    
#    
#    if soil=True:
#        for layer in [0,1,2]:
#            data = xr_ensemble["Soil_liquid"]
