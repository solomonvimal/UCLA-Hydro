#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  2 18:06:40 2017

@author: svimal
"""

import glob, netCDF4, numpy, os
import numpy as np
# Data aggregation

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