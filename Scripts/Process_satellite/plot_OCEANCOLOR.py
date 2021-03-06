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
#                  plot_OCEANCOLOR.py
# Python script for plotting of all AVISO and CMEMS data
#
# Input: .mat file containing from MAIN_test.py
#        station_coord.txt and waypoints_coord.txt files
# 
# Output:       .png files with all the variables from .mat files
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
#matplotlib.use('Tkagg')
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
from cftime import utime
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

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

# for interactive tests
interactive=0
if interactive == 1:
    print ('INTERACTIVE TEST IS ON')
    dir_wrk='/home/gcambon/HCONFIGS_SPASSO/RESILIENCE/Wrk'
    dir_cruise='/home/gcambon/HCONFIGS_SPASSO/RESILIENCE'
    plt.ion()
    plt.close('all')
    print('DIR_WRK:',dir_wrk)

# to load cruise configuration
filelist_domain = glob.glob(dir_cruise+'/domain_limits*.py')
exec(open(dir_cruise+"/cruise_params.py").read())
fileext=['','_zoom']
reso_meridians = [3,2]

#to add SSH to the plot
filelist_ssh=glob.glob(dir_wrk+'/*allsat_phy*.mat')
print(filelist_ssh)
filemat=filelist_ssh[0]
data_aviso=sio.loadmat(filemat)
#  data
adt = data_aviso['adt']
longitudes=data_aviso['lon'][0]
latitudes=data_aviso['lat'][0]
#  datagrid for lon [-180;180] (example: mediterranean product)
adt_newgrid = adt 
lon_newgrid, lat_newgrid = np.meshgrid(longitudes, latitudes)
adt_newgrid = np.ma.masked_where(np.isnan(adt_newgrid), adt_newgrid)

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
    
for i in range(0,2):

    count_domain = 0
    for file_domain in filelist_domain: 
        exec(open(file_domain).read())
        
        #define the geographical projection with Basemap
        mymap=Basemap(projection='merc',llcrnrlat=Lat[0],urcrnrlat=Lat[1],llcrnrlon=Lon[0],urcrnrlon=Lon[1],resolution='h')

        #project the stations on the figure axis
        x_stations,y_stations = mymap(lon_stations,lat_stations)

        #project the waypoint on the figure axis
        x_waypoint,y_waypoint = mymap(lon_waypoint,lat_waypoint)

        # #project the glider on the figure axis
        # x_gl,y_gl = mymap(lon_glider,lat_glider)

        # #project the ZEE limits on the figure axis
        # x_zee,y_zee = mymap(lon_zee,lat_zee)
        # x_zee_sp,y_zee_sp = mymap(lon_zee_sp,lat_zee_sp)
        lon_zee = zee[:,0]
        lat_zee = zee[:,1]
        x_zee,y_zee = mymap(lon_zee,lat_zee)
        
        # #project the S3B trajectories on the figure axis
        # x_extra,y_extra = mymap(lon_extra,lat_extra)


        if (i==0):
        #############
        # L4 CHL DATA
        #############
            print(' ')
            files=glob.glob(dir_wrk+'/*d-ACRI-L4-CHL-MULTI_4KM-GLO-NRT.mat')
            print(files)
            
            if (len(files)>0):
                #define the figure setup
                fig=plt.figure()
                variables=sio.loadmat(files[0])
                lonv=variables['lon']
                latv=variables['lat']
                Chl=np.squeeze(variables['Chl'])
                if (len(lonv)>0):
                    long,latg=np.meshgrid(lonv,latv)
                    (x,y)=mymap(long,latg)
                    
                Chl = np.ma.masked_where(np.isnan(Chl), Chl)
                cax1=mymap.pcolormesh(x,y,Chl,cmap=cmaps.viridis,zorder=-1) 
                
                if (len(chlmin) >= count_domain+1 and len(chlmax) >= count_domain+1):
                    cax1.set_clim(chlmin[count_domain],chlmax[count_domain]) 
                else:
                    cax1.set_clim(chlmin[len(chlmin)-1],chlmax[len(chlmax)-1])
                    
                mymap.drawcoastlines()

                #Add ssh contour of the day
                #  clevels_adt
                clevels_adt=np.arange(adtmin[count_domain],adtmax[count_domain]+adt_int,adt_int)
                #  contour
                x_newgrid, y_newgrid = mymap(lon_newgrid, lat_newgrid)
                cax2=mymap.contour(x_newgrid,y_newgrid,adt_newgrid*100,levels=clevels_adt,colors='grey',zorder=1)
                plt.clabel(cax2,cax2.levels,inline=True, fmt='%1.0f', fontsize=9,colors='black')
                
                #define the tick and the grid (done for global)
                #from south pole to nortehn pole with a resolution of 1 degree
                myparallels=range(-90,90+1,1)
                mymeridians = np.arange(-180,180+1,reso_meridians[count_domain]) 

                #draw the ticks labels
                mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
                mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)

                #draw the continent (data from Basemap)
                mymap.fillcontinents(color='0.83',lake_color='0.83',zorder=100) 

                #draw the stations
                mymap.plot(x_stations,y_stations,'-',color='r',zorder=1)
                if count_domain>0:
                    mymap.scatter(x_stations,y_stations,s=15,color='r',zorder=1) 

                # draw the waypoint
                #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                # draw the glider trajectory
                #mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)

                # # draw the ZEE limits
                # mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1)
                mymap.plot(x_zee,y_zee,'--',color='firebrick',lw=2,zorder=1)
                
                # # draw the S3B trajectories
                # mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1) 
                # mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
                # mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
                # mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)

                # add the colorbar
                cbar1=fig.colorbar(cax1, orientation='vertical',shrink=0.9) 
                cbar1.ax.set_ylabel('Chl [$\mu$g/L]') 

                # add the title
                filefig1=files[0]+fileext[count_domain]+'.png'
                print('filefig is '+ filefig1)
                titlefig1=files[0]
                plt.title(titlefig1[len(dir_wrk)+1:len(dir_wrk)+9]+' - CHL MULTI - L4 - 4KM')
                plt.savefig(filefig1)
                
                #filefig1_d=dir_wrk+'/oftheday/'+titlefig1[len(dir_wrk)+10:-8].replace('.','_')+'_mat'+fileext[count_domain]+'.png'
                filefig1_d=dir_wrk+'/oftheday/'+titlefig1[len(dir_wrk)+10:-8]+'.mat'+fileext[count_domain]+'.png'
                print(filefig1_d)
                plt.savefig(filefig1_d)

        #############
        # L3 CHL DATA
        #############
        if (i==1):
            print(' ')
            files=glob.glob(dir_wrk+'/*d-ACRI-L3-CHL-MULTI_4KM-GLO-NRT.mat')
            print(files)
            
            if (len(files)>0):
                #define the figure setup
                fig=plt.figure()
                variables=sio.loadmat(files[0])
                lonv=variables['lon']
                latv=variables['lat']
                Chl=np.squeeze(variables['Chl'])
                if (len(lonv)>0):
                    long,latg=np.meshgrid(lonv,latv)
                    (x,y)=mymap(long,latg)
                    
                Chl = np.ma.masked_where(np.isnan(Chl), Chl)
                cax1=mymap.pcolormesh(x,y,Chl,cmap=cmaps.viridis,zorder=-1) 
                
                if (len(chlmin) >= count_domain+1 and len(chlmax) >= count_domain+1):
                    cax1.set_clim(chlmin[count_domain],chlmax[count_domain]) 
                else:
                    cax1.set_clim(chlmin[len(chlmin)-1],chlmax[len(chlmax)-1])
                    
                mymap.drawcoastlines()

                #Add ssh contour of the day
                #  clevels_adt
                clevels_adt=np.arange(adtmin[count_domain],adtmax[count_domain]+adt_int,adt_int)
                #  contour
                x_newgrid, y_newgrid = mymap(lon_newgrid, lat_newgrid)
                cax2=mymap.contour(x_newgrid,y_newgrid,adt_newgrid*100,levels=clevels_adt,colors='grey',zorder=1)
                plt.clabel(cax2,cax2.levels,inline=True, fmt='%1.0f', fontsize=9,colors='black')
                #define the tick and the grid (done for global)
                #from south pole to nortehn pole with a resolution of 1 degree
                myparallels=range(-90,90+1,1)
                mymeridians = np.arange(-180,180+1,reso_meridians[count_domain]) 
                
                #draw the ticks labels
                mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
                mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
                
                #draw the continent (data from Basemap)
                mymap.fillcontinents(color='0.83',lake_color='0.83',zorder=100) 
                
                #draw the stations
                mymap.plot(x_stations,y_stations,'-',color='r',zorder=1)
                if count_domain>0:
                    mymap.scatter(x_stations,y_stations,s=15,color='r',zorder=1) 
                    
                # draw the waypoint
                #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                # draw the glider trajectory
                #mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)
                
                # # draw the ZEE limits
                # mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1)
                mymap.plot(x_zee,y_zee,'--',color='firebrick',lw=2,zorder=1)
                
                # # draw the S3B trajectories
                # mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1) 
                # mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
                # mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
                # mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)
                
                # add the colorbar
                cbar1=fig.colorbar(cax1, orientation='vertical',shrink=0.9) 
                cbar1.ax.set_ylabel('Chl [$\mu$g/L]') 
                
                # add the title
                filefig1=files[0]+fileext[count_domain]+'.png'
                print('filefig is '+ filefig1)
                titlefig1=files[0]
                plt.title(titlefig1[len(dir_wrk)+1:len(dir_wrk)+9]+' - CHL MULTI - L3 - 4KM')
                plt.savefig(filefig1)
                
                #filefig1_d=dir_wrk+'/oftheday/'+titlefig1[len(dir_wrk)+10:-8].replace('.','_')+'_mat'+fileext[count_domain]+'.png'
                filefig1_d=dir_wrk+'/oftheday/'+titlefig1[len(dir_wrk)+10:-8]+'.mat'+fileext[count_domain]+'.png'
                print(filefig1_d)
                plt.savefig(filefig1_d)

        ######################
        ######################                   
        count_domain = count_domain+1

























