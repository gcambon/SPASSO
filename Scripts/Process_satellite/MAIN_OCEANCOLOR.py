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
#                  MAIN_OCEANCOLOR.py
#
# Python script for convertion of oceancolor CHL data
# Input: CMEMS *_CHL*.nc files
#
# Output: .mat file containing Longitude (-180 to 180) for mediterranean products, Latitude and SST[degC] 
# Output: .mat file containing Longitude (0 to 360) for global products, Latitude and SST[degC] 
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
import sys 
from netCDF4 import Dataset
import cmocean as cm_oc
import colormaps as cmaps
import csv
import matplotlib.rcsetup as rcsetup
import numpy as np
import scipy.io as sio 
import matplotlib.pyplot as plt
import xarray as xr
import matplotlib as mpl 
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
import matplotlib 
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import glob 
from mpl_toolkits.basemap import Basemap 
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
#from cftime import utime
#from netcdftime import utime
import cftime
import time
import datetime 
#############################################
args=sys.argv
### TEST MODE 
if len(args)==1:
    dir_wrk='/home/gcambon/HCONFIGS_SPASSO/DEMO2/Wrk'
    dir_cruise='/home/gcambon/HCONFIGS_SPASSO/DEMO2'
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
filelist=glob.glob(dir_wrk+'/*_d-ACRI-L4-CHL-MULTI_4KM-GLO-NRT.nc');
filenc=filelist[0]
filemat=filenc[0:-2]+'mat'
print(filenc)

#Load netcdf data
ncdata=xr.open_dataset(filenc,decode_times='False')

#  extract for the domain
data_domain = ncdata.sel(lon=slice(Lon[0],Lon[1]), lat=slice(Lat[1],Lat[0]))

#  save the data in a mat file
output = {'lon':data_domain['lon'].values, 
          'lat':np.flipud(data_domain['lat'].values),
          'time':data_domain['time'].values.astype("datetime64[ns]"), 
          'Chl':np.fliplr(data_domain['CHL'].values)}

#  save the data in a mat file
sio.matlab.savemat(filemat, output)

#############
### L3 GLOBAL MULTI 4km resolution
#############
filelist=glob.glob(dir_wrk+'/*_d-ACRI-L3-CHL-MULTI_4KM-GLO-NRT.nc');
filenc=filelist[0]
filemat=filenc[0:-2]+'mat'
print(filenc)
    
#  Load netcdf data
ncdata=xr.open_dataset(filenc,decode_times='False')
    
#  extract for the domain
data_domain = ncdata.sel(lon=slice(Lon[0],Lon[1]), lat=slice(Lat[1],Lat[0]))
    
#  save the data in a mat file
output = {'lon':data_domain['lon'].values,
          'lat':np.flipud(data_domain['lat'].values),
          'time':data_domain['time'].values.astype("datetime64[ns]"),
          'Chl':np.fliplr(data_domain['CHL'].values)}
#  save the data in a mat file
sio.matlab.savemat(filemat, output)
