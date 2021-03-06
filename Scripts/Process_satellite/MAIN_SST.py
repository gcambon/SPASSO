##############################################################
##                   SPASSO                  vs 0.0.0beta  ###
##############################################################
##  This  file  is  free software ;  it  is distributed  in  #
##  the hope that  it  will  be  useful,  but  without  any  #
##  warranty.  You  can  redistribute  it and/or modify  it  #
##  under  the  terms  of  the  GNU  General Public License  #
##  as  published  by  the   Free  Software  Foundation  at  #
##       http://www.gnu.org/copyleft/gpl.html                #
##############################################################
#                  MAIN_SST.py
#
# Python script for convertion of Seas Surface Temperature data
# Input: CMEMS *_SST*.nc & JPL ourocean *_SST*.nc files
#
# Output: .mat file containing Longitude (-180 to 180), Latitude and SST[degC] 
#
##############################################################
# Credits: F.d'Ovidio
#          F.Nencioli
#          L.Rousselet
#          A.Doglioli
#          Alice Della Penna
#          Stephanie Barrillon
#          Anne Petrenko
#          Anais Ricout
###############################################################
import os

import sys #to have argument in the script call
from netCDF4 import Dataset
import cmocean as cm_oc
import colormaps as cmaps
import csv
import matplotlib.rcsetup as rcsetup
import numpy as np
import scipy.io as sio #to read .mat files
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl #to plot graphics
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
import matplotlib 
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import glob #to obtain file paths
from mpl_toolkits.basemap import Basemap #to plot maps
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from netcdftime import utime
import cftime
import time
import datetime 
#############################################
args=sys.argv
### TEST MODE 
if len(args)==1:
    dir_wrk='/home/SPASSO/Cruises/DEMO/Wrk'
    dir_cruise='/home/SPASSO/Cruises/DEMO'
### OPERATIONAL MODE
else:
    dir_wrk=sys.argv[1]
    dir_cruise=sys.argv[2]
### END if len(args)<1:

interactive=0
if interactive == 1:
    print ('INTERACTIVE TEST IS ON')
    dir_wrk='/home/gcambon/HCONFIGS_SPASSO/RESILIENCE/Wrk'
    dir_cruise='/home/gcambon/HCONFIGS_SPASSO/RESILIENCE'

# to load cruise configuration
exec(open(dir_cruise+"/domain_limits.py").read())

######## CHL FILES #####################
#############
### L4 GLOBAL MULTI 4km resolution
#############
filelist=glob.glob(dir_wrk+'/*120000-UKMO-L4_GHRSST-SSTfnd-GMPE-GLOB-v03.0-fv03.0.nc');
filenc=filelist[0]
filemat=filenc[0:-2]+'mat'
print(filenc)

# Load netcdf data
ncdata=xr.open_dataset(filenc,decode_times='False')

# extract for the domain
data_domain = ncdata.sel(lon=slice(Lon[0],Lon[1]), lat=slice(Lat[0],Lat[1]))

# save the data in a mat file
output = {'lon':data_domain['lon'].values,
          'lat':data_domain['lat'].values,
          'time':data_domain['time'].values.astype("datetime64[ns]"),
          'sst':data_domain['analysed_sst'].values-273.15}

sio.matlab.savemat(filemat, output)

#############
### L3 GLOBAL 
#############
filelist=glob.glob(dir_wrk+'/*-IFR-L3C_GHRSST-SSTsubskin-ODYSSEA-GLOB_010_adjusted-v2.0-fv1.0.nc');
filenc=filelist[0]
filemat=filenc[0:-2]+'mat'
print(filenc)

# Load netcdf data
ncdata=xr.open_dataset(filenc,decode_times='False')

# extract for the domain
data_domain = ncdata.sel(lon=slice(Lon[0],Lon[1]), lat=slice(Lat[0],Lat[1]))

# save the data in a mat file
output = {'lon':data_domain['lon'].values,
          'lat':data_domain['lat'].values,
          'time':data_domain['time'].values.astype("datetime64[ns]"),
          'sst':data_domain['sea_surface_temperature'].values-273.15}

sio.matlab.savemat(filemat, output)


# ### JPL Global SST L4
# files=glob.glob(dir_wrk+'/*-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc')
# for nc_name in files:
#     filemat=nc_name[0:-2]+'mat'

#     # Load netcdf data
#     ncdata=xr.open_dataset(nc_name)

#     # extract for the domain
#     data_domain = ncdata.sel(lon=slice(Lon[0],Lon[1]), lat=slice(Lat[0],Lat[1]))

#     # save the data in a mat file
#     output = {'lon':data_domain['lon'].values,
#               'lat':data_domain['lat'].values,
#               'time':data_domain['time'].values,
#               'sst':data_domain['analysed_sst'].values-273.15}
#     sio.matlab.savemat(filemat, output)

