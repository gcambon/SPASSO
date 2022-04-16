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
#                  plot_ALTI.py
# Python script for plotting of all AVISO and CMEMS data
#
# Input: .mat file containing from MAIN_test.py
#        station_coord.txt and waypoints_coord.txt files
# 
# Output: 	.png files with all the variables from .mat files
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
import matplotlib

matplotlib.use('Agg')
##matplotlib.use('Tkagg') # interactive only

import sys 
import cmocean as cm_oc
import colormaps as cmaps
import csv 
import matplotlib.rcsetup as rcsetup
import numpy as np
import xarray as xr
import scipy.io as sio 
import matplotlib.pyplot as plt
import matplotlib as mpl 
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
import matplotlib 
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import glob 
from mpl_toolkits.basemap import Basemap, shiftgrid 
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)


#from netcdftime import utime
#from cftime import utime
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
    #--END if len(args)<1:

## for interactive test
interactive=0
if interactive == 1:
    print ('INTERACTIVE TEST IS ON')
    dir_wrk='/home/gcambon/HCONFIGS_SPASSO/RESILIENCE/Wrk'
    dir_cruise='/home/gcambon/HCONFIGS_SPASSO/RESILIENCE'
    plt.ion()
    plt.close('all')
    print('DIR_WRK:',dir_wrk)
    


#  Load data
filelist=glob.glob(dir_wrk+'/*allsat_phy*.mat')
print(filelist)
filemat=filelist[0]


filefig=[filemat[0:-3]+'png', filemat[0:-4]+'_zoom.png'] 
filefig_d=[dir_wrk+'/oftheday/'+filemat[len(dir_wrk):-22]+'.png', dir_wrk+'/oftheday/'+filemat[len(dir_wrk):-22]+'_zoom.png'] 

reso_meridians = [3,2] 

# to load cruise configuration
filelist_domain = glob.glob(dir_cruise+'/domain_limits*.py')
exec(open(dir_cruise+"/cruise_params.py").read())

# load the SST data and extract lon and lat from the projection
data_aviso=sio.loadmat(filemat)
lon = data_aviso['lon'][0]
lat = data_aviso['lat'][0]

#load the cruise stations
dico_stations=open(dir_cruise+'/station_coord.txt')
lon_stations=[]
lat_stations=[]
for line in dico_stations:
    xx,yy=line.split()
    lat_stations.append(float(xx))
    lon_stations.append(float(yy))

#load the cruise waypoints
dico_waypoint=open(dir_cruise+'/waypoint_coord.txt')
lon_waypoint=[]
lat_waypoint=[]
for line in dico_waypoint:
    xxx,yyy=line.split()
    lon_waypoint.append(float(xxx))
    lat_waypoint.append(float(yyy))
    
# #load the glider trajectory
# dico_glider_lon=open(dir_cruise+'/gliderLon.txt') 
# dico_glider_lat=open(dir_cruise+'/gliderLat.txt')

# lon_glider=[]
# lat_glider=[]
# for line in dico_glider_lon:
#     x_gl= float(line)
#     lon_glider.append(float(x_gl))
# for line in dico_glider_lat:
#     y_gl= float(line)
#     lat_glider.append(float(y_gl))

# #load the ZEE limits
# dico_zee=open(dir_cruise+'/Algeria.txt')
# lon_zee=[]
# lat_zee=[]
# for line in dico_zee:
#     xxxx,yyyy=line.split()
#     lon_zee.append(float(xxxx))
#     lat_zee.append(float(yyyy))

# dico_zee_spain=open(dir_cruise+'/Spain.txt')
# lon_zee_sp=[]
# lat_zee_sp=[]
# for line in dico_zee_spain:
#     xxxx_sp,yyyy_sp=line.split()
#     lon_zee_sp.append(float(xxxx_sp))
#     lat_zee_sp.append(float(yyyy_sp))

# dico_zee_spain=open(dir_cruise+'/Italy_ZEE_allwithoutcoast.txt')
# lon_zee_sp=[]
# lat_zee_sp=[]
# for line in dico_zee_spain:
#     xxxx_sp,yyyy_sp=line.split()
#     lon_zee_sp.append(float(xxxx_sp))
#     lat_zee_sp.append(float(yyyy_sp))

#load the SWOT trajectories  
#dico_extra=open(dir_cruise+'/extra_coord.txt')
#lon_extra=[]
#lat_extra=[]
#for line in dico_extra:
#        xxxxx,yyyyy=line.split()
#       lon_extra.append(float(xxxxx))
#       lat_extra.append(float(yyyyy))

# #load the S3B trajectories
# dico_extra=open(dir_cruise+'/S3B.txt') 
# lon_extra=[]
# lat_extra=[]
# for line in dico_extra:
#     xxxxx,yyyyy=line.split()
#     lon_extra.append(float(xxxxx))
#     lat_extra.append(float(yyyyy))
    
#data 
lon, lat = np.meshgrid(lon, lat)
adt = data_aviso['adt']
longitudes=data_aviso['lon'] [0]
latitudes=data_aviso['lat'] [0]
u=data_aviso['u']
v=data_aviso['v']


# datagrid for lon [-180;180] (example: mediterranean product)
adt_newgrid = adt 
lon_newgrid, lat_newgrid = np.meshgrid(longitudes, latitudes)
adt_newgrid = np.ma.masked_where(np.isnan(adt_newgrid), adt_newgrid) 
u_newgrid = u 				
v_newgrid = v 				
lon_newgrid_2 = longitudes 	


# datagrid for lon [0;360], shift to the following grid  (example: global product)
#adt_newgrid, lon_newgrid = shiftgrid(180.,adt,longitudes,start=False,cyclic=360.0) 
#lon_newgrid, lat_newgrid = np.meshgrid(longitudes, latitudes)
#adt_newgrid = np.ma.masked_where(np.isnan(adt_newgrid), adt_newgrid) 
#u_newgrid,lon_newgrid_2 = shiftgrid(180.,u,longitudes,start=False) 
#v_newgrid,lon_newgrid_2 = shiftgrid(180.,v,longitudes,start=False) 


count_domain = 0
for file_domain in filelist_domain:
    ##execfile(file_domain)
    exec(open(file_domain).read())
    
    #define the figure setup
    fig=plt.figure()
    #fig,ax=plt.subplots()
    
    #define the geographical projection with Basemap
    mymap=Basemap(projection='merc',llcrnrlat=Lat[0],urcrnrlat=Lat[1],llcrnrlon=Lon[0],urcrnrlon=Lon[1],resolution='h')#equivalent to m_proj

    #project the data on the figure axis
    x, y = mymap(lon, lat) #[km]

    #project the stations on the figure axis
    x_stations,y_stations = mymap(lon_stations,lat_stations)

    #project the waypoint on the figure axis
    x_waypoint,y_waypoint = mymap(lon_waypoint,lat_waypoint)

    # #project the glider on the figure axis
    # x_gl,y_gl = mymap(lon_glider,lat_glider)

    # #project the ZEE limits on the figure axis
    lon_zee = zee[:,0]
    lat_zee = zee[:,1]
    x_zee,y_zee = mymap(lon_zee,lat_zee)
    #x_zee_sp,y_zee_sp = mymap(lon_zee_sp,lat_zee_sp)

    # #project the SWOT trajectories on the figure axis
    # x_extra,y_extra = mymap(lon_extra,lat_extra)        

    x_newgrid, y_newgrid = mymap(lon_newgrid, lat_newgrid)
    cax1=mymap.pcolormesh(x_newgrid,y_newgrid,adt_newgrid*100,cmap=cm_oc.cm.ice,zorder=-1)

    clevels_adt=np.arange(adtmin[count_domain],adtmax[count_domain]+adt_int,adt_int)
    cax2=mymap.contour(x_newgrid,y_newgrid,adt_newgrid*100,levels=clevels_adt,colors='grey')
    plt.clabel(cax2,cax2.levels,inline=True, fmt='%1.0f', fontsize=9,colors='black')
    #ax.clabel(CS, inline=1, fontsize=10, manual=manual_locations)
    
    # # plot wind vectors on projection grid.
    # scale_quiv = 2
    # uproj,vproj,xx,yy = mymap.transform_vector(u_newgrid,v_newgrid,lon_newgrid_2,latitudes,48*scale_quiv,22*scale_quiv,returnxy=True,masked=True) 

    # # now plot.
    # Q = mymap.quiver(xx,yy,uproj,vproj,linewidth=0.05,color='r')

    # # make quiver key.
    # qk = plt.quiverkey(Q, 1.25, 1.1, 0.5, '0.5 m/s', labelpos='N', color='r')

    #add the coastline (data from Basemap)
    mymap.drawcoastlines()

    #define the tick and the grid (done for global)
    #from south pole to northen pole with a resolution of 1 degree 
    myparallels=np.arange(-90,90+1,1)
    mymeridians = np.arange(-180,180+1,reso_meridians[count_domain])

    #draw the ticks labels
    mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
    mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)

    #draw the continent (data from Basemap)
    mymap.fillcontinents(color='0.83',lake_color='0.83',zorder=100) 

    #draw the stations
    mymap.plot(x_stations,y_stations,'-',color='w',zorder=1)

    if count_domain > 0:
        mymap.scatter(x_stations,y_stations,s=15,color='w',zorder=1)
        
    #draw the waypoint
    #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
    #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
    #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)

    # #draw the glider
    # mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1) 

    # #draw the ZEE limits 
    mymap.plot(x_zee,y_zee,'--',color='firebrick',lw=2,zorder=1)
    
    # #draw the swot trajectories
    # mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
    # mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
    # mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
    # mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)

    #add the colorbar
    if (len(adtmin) >= count_domain+1 and len(adtmax) >= count_domain+1):
        cax1.set_clim(adtmin[count_domain],adtmax[count_domain]) 
    else:
        cax1.set_clim(adtmin[len(adtmin)-1],adtmax[len(adtmax)-1])


    cbar1=fig.colorbar(cax1, orientation='vertical',shrink=0.9)
    cbar1.ax.set_ylabel('ADT [cm]')
    
    #add the title
    plt.title(filemat[len(dir_wrk)+1:-4]) 
    
    plt.savefig(filefig[count_domain]) 
    plt.savefig(filefig_d[count_domain])
    

    count_domain = count_domain +1










































