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
#                  MAIN_ALTI.py
#
# Python script for convertion of altimetric ADT and current data
# Input: *_h*.nc and *_uv*.nc files
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
#from netcdftime import utime
import cftime
#from cftime import utime
import time
import datetime 
#############################################
args=sys.argv
### TEST MODE 
if len(args)==1:
    dir_wrk='/home/gcambon/HCONFIGS_SPASSO/DEMO2/Wrk'
    dir_cruise='/home/gcambon/HCONFIGS_SPASSO/DEMO2/'
    ### OPERATIONAL MODE
else:
    dir_wrk=sys.argv[1]
    dir_cruise=sys.argv[2]
    ### END if len(args)<1:
        
        
##### ALTI FILES #####################
filelist=glob.glob(dir_wrk+'/*allsat_phy*.nc'); 
filenc=filelist[0]
filemat=filenc[0:-2]+'mat'
print(filenc)
# Load netcdf data
ncdata=xr.open_dataset(filenc,decode_times='False')
data = ncdata.sel()

# save the data in a mat file
output = {'lon':data['longitude'].values, 'lat':data['latitude'].values,
          'time':data['time'].values,'adt':np.squeeze(data['adt'].values),
          'u':np.squeeze(data['ugos'].values),'v':np.squeeze(data['vgos'].values)}

sio.matlab.savemat(filemat, output)









