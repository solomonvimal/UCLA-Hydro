#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 16:20:52 2017

@author: svimal
"""

import geopandas
from rasterio import features
from affine import Affine
import numpy as np
import xarray as xray

def transform_from_latlon(lat, lon):
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(shapes, coords, latitude='latitude', longitude='longitude',
              fill=np.nan, **kwargs):
    """Rasterize a list of (geometry, fill_value) tuples onto the given
    xray coordinates. This only works for 1d latitude and longitude
    arrays.
    """
    transform = transform_from_latlon(coords[latitude], coords[longitude])
    out_shape = (len(coords[latitude]), len(coords[longitude]))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
    return xray.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))

# this shapefile is from natural earth data
# http://www.naturalearthdata.com/downloads/10m-cultural-vectors/10m-admin-1-states-provinces/

counties = geopandas.read_file("/home2/svimal/Github/CA_drought/data/Spatial/CA_counties/CA_counties_WGS.shp")
county_ids = {k: i for i, k in enumerate(counties.NAME)}
shapes = zip(counties.geometry, range(len(counties)))
ds = xray.open_dataset("/home2/svimal/Data/VIC_Fluxes/1920s/fluxes.1924-11.nc")
ds["counties"] = rasterize(shapes, ds.coords, longitude="lon", latitude="lat")

# Plot of all counties
(ds.counties.sel(lat=slice(32.5, 42.5), lon=slice(-125.2, -114.0)).plot())
  
# Plot of only one county
(ds.counties == county_ids['Los Angeles']).plot()

# Make plots
(ds.Tair
 .isel(time=slice(2))
 .where(ds.counties == county_ids['Orange'])
 .sel(lat=slice(33.4, 34.0), lon=slice(-118.2, -117.4))
 .plot.imshow(col='time', col_wrap=4))

# Plot of time series 
(ds.Tair
 .isel(time=slice(30))
 .where(ds.counties == county_ids['Orange'])
 .sel(lat=slice(32, 43), lon=slice(-122, -112))
 .mean(['lat', 'lon'])
 .plot())

series = (ds.Tair
 .isel(time=slice(30))
 .where(ds.counties == county_ids['Orange'])
 .sel(lat=slice(32, 43), lon=slice(-122, -112))
 .mean(['lat', 'lon']))
