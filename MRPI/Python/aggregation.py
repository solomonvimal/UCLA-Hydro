#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 18:06:40 2017
Data Aggregation 
@author: svimal
"""

import glob, os
import numpy as np
import xarray as xray
from affine import Affine
from rasterio import features
import geopandas
os.getcwd()

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

fluxes = glob.glob("/home2/svimal/Data/VIC_Fluxes/*/*.nc") # Monthly for 86 years
fnames_yearly = sorted(list(set([fname[0:-5]+"??.nc" for fname in fluxes])))
counties = geopandas.read_file("/home2/svimal/Github/CA_drought/data/Spatial/CA_counties/CA_counties_WGS.shp")
county_ids = {k: i for i, k in enumerate(counties.NAME)}
shapes = zip(counties.geometry, range(len(counties)))


variables = ["Prec", "Evap", "Runoff", "Baseflow", "Tair", "SWE", "Soil_liquid"]

    
def worker(fname):
    year =  fname[-10:-6]
    try:
        os.mkdir("/home2/svimal/Github/UCLA-Hydro/MRPI/Aggregated_Data/"+year)
    except:
        pass
    ds = xray.open_mfdataset(fname) # Open multiple files at a time (12 files for each year)
    #ds.time
    ds["counties"] = rasterize(shapes, ds.coords, longitude="lon", latitude="lat")
    #index = range(len(ds.time))
    #data = {'index': index}
    ds = ds.resample("W", "time", how="mean")
    for week in range(len(ds.time)):
        filename = "/home2/svimal/Github/UCLA-Hydro/MRPI/Aggregated_Data/" + year + "/"+ year +"_week_"+str(week+1) + ".csv"
        # Write the header for the CSV file
        with open(filename, "w") as f:
            f.write("County, Precipitation, Evapotranspiration, Runoff, Baseflow, Air_Temperature, Snow_Water_Equivalent, Soil_Liquid1, Soil_Liquid2, Soil_Liquid3 \n")
        for name in counties.NAME:
            line = []
            for v in variables:
                if v != "Soil_liquid":
                    series = eval("(ds." + v + ".where(ds.counties == county_ids['" + name + "']).mean(['lat', 'lon']))")
                    #.isel(time=slice(30))#.sel(lat=slice(32, 43), lon=slice(-122, -112)).
                    series = list(np.array(series))
                    #data[fname[-10:-6]+"_"+v] = series
                    try:
                        line.append(round(float(series[week])))    
                    except:
                        line.append(series[week])
                else:
                    soil_liquid = []
                    series = eval("(ds." + v + ".where(ds.counties == county_ids['" + name + "']).mean(['lat', 'lon']))")
                    for sl in [0,1,2]:
                        soil_liquid.append(round(float(series.T[sl][week]),3))
                    line.extend(soil_liquid)
            line = str(line).strip("[").strip("]")
            with open(filename, "a+") as f:
                f.write(name + ", " + str(line) + "\n")
    return

for fname in list(reversed(fnames_yearly)):
    worker(fname)
   

'''
from multiprocessing.pool import ThreadPool as Pool
# from multiprocessing import Pool

pool_size = 7  # "parallelness"
pool = Pool(pool_size)

for fname in list(reversed(fnames_yearly)):
    pool.apply_async(worker, (fname,))

pool.close()
pool.join()


os.getcwd()

def write_asc(fname, data_list, filename):
    
    """
    input: netCDF file, and list containing the x,y,z data, and filename
    output: filename
    creates an asc raster file and clips it to california size
    """
    data =  netCDF4.Dataset(fname, "r")
    
    ncols = len(data.variables["lon"][:])
    nrows = len(data.variables["lat"][:])
    xllcorner = float(min(data.variables["lon"][:])) + (-0.0625/2.0)
    yllcorner = float(min(data.variables["lat"][:])) + (-0.0625/2.0)
    cell_size =  abs(data.variables["lat"][1] - data.variables["lat"][0])
    #nodata_value = -9999

    with open(filename, "w") as f:
        f.write("NCOLS " + str(ncols) + "\n")
        f.write("NROWS " + str(nrows) + "\n")
        f.write("XLLCORNER " + str(xllcorner) + "\n") 
        f.write("YLLCORNER " + str(yllcorner) + "\n")
        f.write("CELLSIZE " + str(cell_size) + "\n")
        f.write("NODATA_VALUE = -9999\n")
        data_2d =  np.reshape(data_list, (nrows, ncols))
        for line in data_2d:
            f.write(str(list(line)).strip("[").strip("]")[1:].strip("\n") + "\n")
    return filename


# Path to folders 
forcings = glob.glob("/home2/svimal/Data/VIC_Forcing/*/*.nc") # 2 datasets because its divided into cali1 and cali2
fluxes = glob.glob("/home2/svimal/Data/VIC_Fluxes/*/*.nc") # Monthly for 86 years

# No. of files each
len(forcings), len(fluxes) # Forcing files are stored as two parts for california

# Start with fluxes 

# Select all the files of one year and work on that. 
years_sorted=[]
for year in range(1920,2015+1):
    #print(year)
    for flux_file in fluxes:
        year_c = str(flux_file.split(".")[1].split("-")[0])
        #print(year_c)
        if str(year)==year_c:
            years_sorted.append(flux_file)

year1 = years_sorted[0:12]
year1_sorted = sorted(year1, key=lambda x: x.split(".")[1].split("-")[1])

day = 1
data_lists = [] # data.variables["Soil_liquid"][29][2][349][309] # #data.variables["Soil_liquid"][time][lat][lon]

for fname in year1_sorted:
    data = netCDF4.Dataset(fname, "r")
    for soil_layer in [0,1,2]:
        for time in range(len(data["time"][:])):
            data_list = []
            for lat in list(reversed(range(len(data["lat"][:])))): # reversed to get order of tiff format lat-lon correct
                for lon in range(len(data["lon"][:])):
                    d = float(data.variables["Soil_liquid"][time][soil_layer][lat][lon])
                    if numpy.isnan(float(d)) == True:
                        data_list.append(0)
                    else:
                        data_list.append(d)
            data_lists.append(data_list)
            day = day+1
            if day == 366:
                filename = "/home2/svimal/Data/MRPI_results/" + "Soil_Layer" + str(soil_layer) + "_Week53.asc"
                #data_sum=np.sum(data_lists, axis=0)
                write_asc(fname, data_list, filename)
                data_lists=[]
            elif day/7.0==1.0:
                filename = "/home2/svimal/Data/MRPI_results/" + "Soil_Layer" + str(soil_layer) + "_Week" + str(day/7)+".asc"
                data_sum = np.sum(data_lists, axis=0)
                write_asc(fname, data_sum, filename)
                data_lists = []
            else:
                pass
'''
