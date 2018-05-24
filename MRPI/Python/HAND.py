#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Code to create HAND raster

Non native Python dependency: 
TauDEM - install from https://github.com/dtarb/TauDEM 
After installing TauDEM, to test if it works fine, try typing the following two commands on the terminal 
1. pitremove
2. dinfdistdown

Usage instructions:
On terminal type the following
python HAND.py path/to/dem_file_1.tif path/to/dem_file_2.tif ... path/to/dem_file_n.tif

"""

import os
import subprocess
import time


dem_path = "/home/Downloads/8135001_DTM_actual.tif"



def HAND_from_DEM(dem_path):
    '''
    Args:
      dem: Path the the raw DEM
    Returns:
      The raster files of: 
      - filled elevation rater
      - D8 Flow Direction
      - DInf Flow Direction
      - Flow Accumulation (D8 Specific Catchment Area)
      - Stream raster for Thresholds 100, 500, 1000
      - Height Above nearest Drainage (HAND) raster 
    '''

    file_name = dem_path
    fel = PitRemove(file_name)
    
    try:
        if (fel == file_name[:-4] + "_fel.tif"):
            pass
    except:
        print("Flow direction not done, waiting for a few seconds")
        time.sleep(3)
    d8 = d8FlowDir(fel, file_name);
    dInf = dInfFlowDir(fel, file_name);
    AreaD8 = d8ContribArea(d8, file_name);
    for thresh in [1155]: # 500, 1000, etc.
        stream, threshold = Threshold(AreaD8, thresh, file_name)
        hand_raster = HAND(dInf, fel, stream, threshold, file_name)
    del fel, stream,threshold, AreaD8, thresh, d8, file_name, dInf

def PitRemove(file_name):
    '''
    Takes the raw DEM and runs TauDEM PitRemove function for it and returns fel (filled elevation raster) as output.
    '''
    fel = file_name[:-4] + "_fel.tif"
    command = ['pitremove','-z',file_name,'-fel',fel]
    subprocess.check_output(command)
    print("pitremove is done!")
    return fel

def d8FlowDir(fel, file_name):
    '''
    Creates the D8 flow direction raster from input - filled DEM
    '''
    d8 = file_name[:-4] + "_d8.tif"
    command = ['d8flowdir','-fel',fel,'-p',d8]
    subprocess.check_output(command)
    print("D8 done")
    return d8

def dInfFlowDir(fel, file_name):
    '''
    Creates the D Infinity flow direction raster from input - filled DEM
    '''
    dInf = file_name[:-4] + "_dInf.tif"
    command =  ['dinfflowdir','-fel',fel,"-ang",dInf]
    subprocess.check_output(command)
    print("DInf done")
    return dInf

def d8ContribArea(d8, file_name):
    '''
    Creates the d8 Contributing Area (Flow Accumulation) raster from input D8
    '''
    AreaD8 =  file_name[:-4] + "_AreaD8.tif"
    command =  ['aread8', '-p', d8, '-ad8', AreaD8, '-nc']
    subprocess.check_output(command)
    print("d8Contrib Area is done!")
    return AreaD8

# Stream Definition By Threshold
def Threshold(AreaD8, threshold, file_name):
    '''
    Takes AreaD8 and input thresholds to create a stream raster 
    '''
    stream = file_name[:-4] + "_stream_" + str(threshold) + ".tif"
    command =  ['threshold','-ssa',AreaD8,'-src',stream,'-thresh',str(threshold)]
    subprocess.check_output(command)
    print("thresholding is done")
    return stream, str(threshold)

# DInf Vertical Distance Down
def HAND(dInf, fel, stream, threshold, file_name):
    '''
    Takes as input Inf, fel, stream, threshold to create the HAND raster.
    '''
    HAND_raster = file_name[:-4] + "_HAND_" + str(threshold) + ".tif"
    command =  ['dinfdistdown','-ang',dInf,'-fel',fel,'-src',stream,'-dd',HAND_raster,'-m ave v','-nc']
    subprocess.check_output(command)
    print("hand is done!")
    print(command)
    return HAND_raster



def main():
    # print command line arguments
    for arg in sys.argv[1:]:
        HAND_from_DEM(arg)

if __name__ == "__main__":
    main()
