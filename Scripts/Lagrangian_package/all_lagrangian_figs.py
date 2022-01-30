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
#                  all_lagrangian_figs.py
#
# Python script for plotting Lagrangian diagnostics
#
# Input: .mat file containing output from MAIN_Lagrangian.py
#        station_coord.txt and waypoints_coord.txt files
#
# Output: .png files with maps of (some or all) diagnostics
#	  use 'keys' to decide which ones are plotted 	  
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
import sys 
import cmocean as cm_oc
import csv 
import matplotlib.rcsetup as rcsetup
import numpy as np
import scipy.io as sio 
import matplotlib.pyplot as plt
from matplotlib.dates import DayLocator, HourLocator, DateFormatter, drange
import matplotlib
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
import glob 
from mpl_toolkits.basemap import Basemap, shiftgrid 
from matplotlib.patches import Polygon
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
	print('all_lagrangian_figs operational mode')
#############################################

# Labels for each figure
figure_keys=['NYYYYYN']
figure_labels=['OW_disp','OW_param','FSLE','Lon_adv','Lat_adv','AVISO_vel','Time_From_Bathy [day]'];

# Load data
filename=dir_wrk+'/[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_*d[0-9].mat'
filelist=glob.glob(filename)
filemat=filelist[0]


if len(filelist)<0: 
   print('FILE: '+filename+'.mat NOT FOUND!')

filelist_domain = sorted(glob.glob(dir_cruise+'/domain_limits*.py'))
execfile(dir_cruise+"/cruise_params.py") 
fileext=['','_zoom']
reso_meridians = [2,1] 

# Extracts all the common variables 
variables=sio.loadmat(filemat)
long=variables['long']
latg=variables['latg']
day0=variables['day0'] 
filefigbase="%02d%02d%02d" % (day0[0,0],day0[0,1],day0[0,2]) 

# load the cruise stations
dico_stations=open(dir_cruise+'/station_coord.txt')
lon_stations=[]
lat_stations=[]
for line in dico_stations:
    xx,yy=line.split()
    lat_stations.append(float(xx))
    lon_stations.append(float(yy))

# load the cruise waypoints
dico_waypoint=open(dir_cruise+'/waypoint_coord.txt')
lon_waypoint=[]
lat_waypoint=[]
for line in dico_waypoint:
    xxx,yyy=line.split()
    lon_waypoint.append(float(xxx))
    lat_waypoint.append(float(yyy))

# load the glider trajectory
dico_glider_lon=open('/home/glider/realTimePosition/gliderLon.txt')
dico_glider_lat=open('/home/glider/realTimePosition/gliderLat.txt')
lon_glider=[]
lat_glider=[]
for line in dico_glider_lon:
    x_gl= float(line)
    lon_glider.append(float(x_gl))
for line in dico_glider_lat:
    y_gl= float(line)
    lat_glider.append(float(y_gl))
    
# load the ZEE limits
dico_zee=open(dir_cruise+'/Algeria.txt')
lon_zee=[]
lat_zee=[]
for line in dico_zee:
    xxxx,yyyy=line.split()
    lon_zee.append(float(xxxx))
    lat_zee.append(float(yyyy))

dico_zee_spain=open(dir_cruise+'/Spain.txt')
lon_zee_sp=[]
lat_zee_sp=[]
for line in dico_zee_spain:
    xxxx_sp,yyyy_sp=line.split()
    lon_zee_sp.append(float(xxxx_sp))
    lat_zee_sp.append(float(yyyy_sp))

dico_zee_spain=open(dir_cruise+'/Italy_ZEE_allwithoutcoast.txt') 
lon_zee_sp=[]
lat_zee_sp=[]
for line in dico_zee_spain:
        xxxx_sp,yyyy_sp=line.split()
        lon_zee_sp.append(float(xxxx_sp))
        lat_zee_sp.append(float(yyyy_sp))
    
# load the SWOT trajectories  
#dico_extra=open(dir_cruise+'/extra_coord.txt')
#lon_extra=[]
#lat_extra=[]
#for line in dico_extra:
#        xxxxx,yyyyy=line.split()
#       lon_extra.append(float(xxxxx))
#       lat_extra.append(float(yyyyy))

# load the S3B trajectories
dico_extra=open(dir_cruise+'/S3B.txt')
lon_extra=[]
lat_extra=[]
for line in dico_extra:
        xxxxx,yyyyy=line.split()
        lon_extra.append(float(xxxxx))
        lat_extra.append(float(yyyyy))

# Defines the map
for ii in range(0,7):
        count_domain = 0
	for file_domain in filelist_domain:
                execfile(file_domain)
                mymap=Basemap(projection='merc',llcrnrlat=Lat[0],urcrnrlat=Lat[1],llcrnrlon=Lon[0],urcrnrlon=Lon[1],resolution='h')
                # project the data on the figure axis
                x, y = mymap(long, latg) #[km]
                # project the stations on the figure axis
                x_stations,y_stations = mymap(lon_stations,lat_stations)
                # project the waypoint on the figure axis
                x_waypoint,y_waypoint = mymap(lon_waypoint,lat_waypoint)
                # project the glider on the figure axis
                x_gl,y_gl = mymap(lon_glider,lat_glider)
                # project the ZEE limits on the figure axis
                x_zee,y_zee = mymap(lon_zee,lat_zee)
                x_zee_sp,y_zee_sp = mymap(lon_zee_sp,lat_zee_sp)
                # project the S3B trajectories on the figure axis
                x_extra,y_extra = mymap(lon_extra,lat_extra)
                
	        if figure_keys[0][ii] in 'Y':
		
		        if ii==0:
			        OWdispersion_bi=variables['OWdispersion_bi']
                                OWdispersion_bi = np.ma.masked_where(np.isnan(OWdispersion_bi), OWdispersion_bi)
			        fig1=plt.figure()
			        cax1=mymap.pcolormesh(x,y,OW_dispersion_bi,cmap=cm_oc.cm.deep_r,zorder=-1)
			        # add the coastline (data from Basemap)
			        mymap.drawcoastlines()	
			        # define the tick and the grid (done for global)
                                # from south pole to nortehn pole with a resolution of 1 degree
			        myparallels=np.arange(-90,90+1,1)
                                mymeridians = np.arange(-180,180+1,reso_meridians[count_domain]) 
			        # draw the ticks labels
			        mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
			        mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
			        # draw the continent (data from Basemap)
			        mymap.fillcontinents(color='0.6',lake_color='0.6',zorder=100)
			        # draw the stations
			        mymap.plot(x_stations,y_stations,'-',color='w',zorder=1)
	                        if count_domain>0:
                                        mymap.scatter(x_stations,y_stations,s=15,color='w',zorder=1) 
			        # draw the waypoint
			        #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                                # draw the glider
				mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)
                                # draw the ZEE limits
	                        mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1)
			        # draw the swot trajectories
        			mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
        			mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)
			        # add the colorbar
			        cbar1=fig1.colorbar(cax1, orientation='vertical',shrink=0.9)
			        cax1.set_clim(1,20)
			        cbar1=fig1.colorbar(cax1,ticks=[5,10,15,20], orientation='vertical',shrink=0.9)
			        cbar1.ax.set_ylabel('Eddy retention [d]')
			        # add the title
			        filefig1=dir_wrk+'/'+figure_labels[ii]+'_'+filefigbase+fileext[count_domain]+'.png' 
			        plt.title(figure_labels[ii]+'_'+filefigbase)
			        plt.savefig(filefig1)
			        filefig1_d=dir_wrk+'/oftheday/'+figure_labels[ii]+fileext[count_domain]+'.png'
			        plt.savefig(filefig1_d)

		
		        if ii==1:
			        owmi=variables['owmi']
                                owmi = np.ma.masked_where(np.isnan(owmi), owmi) 
			        fig2=plt.figure()
			        cax2=mymap.pcolormesh(x,y,(60*60*24)**2*owmi,cmap=cm_oc.cm.balance,zorder=-1)
			        # add the coastline (data from Basemap)
			        mymap.drawcoastlines()	
			        # define the tick and the grid (done for global)
                                # from south pole to nortehn pole with a resolution of 1 degree
			        myparallels=np.arange(-90,90+1,1)
			        mymeridians = np.arange(-180,180+1,reso_meridians[count_domain])
			        mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
			        mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
			        # draw the continent (data from Basemap)
			        mymap.fillcontinents(color='0.6',lake_color='0.6',zorder=100)
			        # draw the stations
			        mymap.plot(x_stations,y_stations,'-',color='k',zorder=1)
	                        if count_domain>0:
                                        mymap.scatter(x_stations,y_stations,s=15,color='k',zorder=1)
			        # draw the waypoint
			        #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                                # draw the glider
                                mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)
                                # draw the ZEE limits
	                        mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1) 	
			        # draw the swot trajectories
        			mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
        			mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)
			        # add the colorbar
			        cbar2=fig2.colorbar(cax2, orientation='vertical',shrink=0.9)
			        if (len(owmin) >= count_domain+1 and len(owmax) >= count_domain+1):
                                        cax2.set_clim(owmin[count_domain],owmax[count_domain])
                                else:
                                        cax2.set_clim(owmin[len(owmin)-1],owmax[len(owmax)-1])
   
			        cbar2.ax.set_ylabel('Okubo-Weiss [d$^{-2}$]')
			        # add the title
			        filefig2=dir_wrk+'/'+figure_labels[ii]+'_'+filefigbase+fileext[count_domain]+'.png'
			        plt.title(figure_labels[ii]+'_'+filefigbase)
			        plt.savefig(filefig2)
	                        filefig2_d=dir_wrk+'/oftheday/'+figure_labels[ii]+fileext[count_domain]+'.png'
			        plt.savefig(filefig2_d)
                                
		        if ii==2:
			        fsle=variables['lambda']
			        fig3=plt.figure()
                                fsle = np.ma.masked_where(np.isnan(fsle), fsle) 
			        cax3=mymap.pcolormesh(x,y,-60*60*24*fsle,cmap=cm_oc.cm.thermal,zorder=-1)
			        # add the coastline (data from Basemap)
			        mymap.drawcoastlines()	
			        # define the tick and the grid (done for global)
                                # from south pole to nortehn pole with a resolution of 1 degree
			        myparallels=np.arange(-90,90+1,1)
			        mymeridians = np.arange(-180,180+1,2)
			        mymeridians = np.arange(-180,180+1,reso_meridians[count_domain]) 
			        mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
			        mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
			        # draw the continent (data from Basemap)
			        mymap.fillcontinents(color='0.6',lake_color='0.6',zorder=100)
			        # draw the stations
			        mymap.plot(x_stations,y_stations,'-',color='w',zorder=1)
	                        if count_domain>0:
                                        mymap.scatter(x_stations,y_stations,s=15,color='w',zorder=1) 
			        # draw the waypoint
			        #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                                # draw the glider
				mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)
                                # draw the ZEE limits
	                        mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1) 
			        # draw the swot trajectories
        			mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
        			mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)
			        # add the colorbar
                                if (len(fslemin) >= count_domain+1 and len(fslemax) >= count_domain+1):
                                        cax3.set_clim(fslemin[count_domain],fslemax[count_domain])
                                else:
                                        cax3.set_clim(fslemin[len(fslemin)-1],fslemax[len(fslemax)-1]) 
                                                
			        cbar3=fig3.colorbar(cax3, orientation='vertical',shrink=0.9)
			        cbar3.ax.set_ylabel('FSLE [d$^{-1}$]')
			        # add the title			       
			        filefig3=dir_wrk+'/'+figure_labels[ii]+'_'+filefigbase+fileext[count_domain]+'.png' 
			        plt.title(figure_labels[ii]+'_'+filefigbase)
			        plt.savefig(filefig3)
                                filefig3_d=dir_wrk+'/oftheday/'+figure_labels[ii]+fileext[count_domain]+'.png'
			        plt.savefig(filefig3_d)
                                
		        if ii==3:
			        lonf15=variables['lonf15']
			        fig4=plt.figure()
                                lonf15 = np.ma.masked_where(np.isnan(lonf15), lonf15) 
                                palette = cmap=cm_oc.cm.delta 
                                palette.set_bad('0.75',1.0) 
			        cax4=mymap.pcolormesh(x,y,long-lonf15,cmap=palette,zorder=-1) 
			        # add the coastline (data from Basemap)
			        mymap.drawcoastlines()	
			        # define the tick and the grid (done for global)
                                # from south pole to nortehn pole with a resolution of 1 degree
			        myparallels=np.arange(-90,90+1,1)
			        mymeridians = np.arange(-180,180+1,reso_meridians[count_domain]) 
                                # draw the ticks labels
			        mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
			        mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
			        # draw the continent (data from Basemap)
			        mymap.fillcontinents(color='0.6',lake_color='0.6',zorder=100)
			        # draw the stations
			        mymap.plot(x_stations,y_stations,'-',color='k',zorder=1)
	                        if count_domain>0:
                                        mymap.scatter(x_stations,y_stations,s=15,color='k',zorder=1) 
			        # draw the waypoint
			        #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                                # draw the glider
				mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)     
                                # draw the ZEE limits
	                        mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1)
			        # draw the swot trajectories
        			mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
        			mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
			        # add the colorbar
			        if (len(deltalonmin) >= count_domain+1 and len(deltalonmax) >= count_domain+1):
                                        cax4.set_clim(deltalonmin[count_domain],deltalonmax[count_domain]) 
                                else:
                                        cax4.set_clim(deltalonmin[len(deltalonmin)-1],deltalonmax[len(deltalonmax)-1])
			        cbar4=fig4.colorbar(cax4, orientation='vertical',shrink=0.9)
			        cbar4.ax.set_ylabel('$\Delta$Longitude (15d) [$^\circ$]')  
			        # add the title			        
			        filefig4=dir_wrk+'/'+figure_labels[ii]+'_'+filefigbase+fileext[count_domain]+'.png'
			        plt.title(figure_labels[ii]+'_'+filefigbase)
			        plt.savefig(filefig4)	
                                filefig4_d=dir_wrk+'/oftheday/'+figure_labels[ii]+fileext[count_domain]+'.png'
			        plt.savefig(filefig4_d)
                                
		        if ii==4:
			        latf15=variables['latf15']
			        fig5=plt.figure()
                                latf15 = np.ma.masked_where(np.isnan(latf15), latf15) 
                                palette = cmap=cm_oc.cm.curl 
                                palette.set_bad('0.75',1.0)
			        cax5=mymap.pcolormesh(x,y,latg-latf15,cmap=palette,zorder=-1)
			        # add the coastline (data from Basemap)
			        mymap.drawcoastlines()	
			        # define the tick and the grid (done for global)
                                # from south pole to nortehn pole with a resolution of 1 degree
			        myparallels=np.arange(-90,90+1,1)
			        mymeridians = np.arange(-180,180+1,reso_meridians[count_domain])
			        mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
			        mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
			        # draw the continent (data from Basemap)
			        mymap.fillcontinents(color='0.6',lake_color='0.6',zorder=100)
			        # draw the stations
			        mymap.plot(x_stations,y_stations,'-',color='k',zorder=1)
	                        if count_domain>0:
                                        mymap.scatter(x_stations,y_stations,s=15,color='k',zorder=1) 
			        # draw the waypoint
			        #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                                # draw the glider
				mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)
                                # draw the ZEE limits
	                        mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1)
			        # draw the swot trajectories
        			mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
        			mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)
			        # add the colorbar
                                if (len(deltalatmin) >= count_domain+1 and len(deltalatmax) >= count_domain+1):
			                cax5.set_clim(deltalatmin[count_domain],deltalatmax[count_domain])
                                else:
                                        cax5.set_clim(deltalatmin[len(deltalatmin)-1],deltalatmax[len(deltalatmax)-1]) 
			        cbar5=fig5.colorbar(cax5, orientation='vertical',shrink=0.9)
			        cbar5.ax.set_ylabel('$\Delta$Latitude (15d) [$^\circ$]') 
			        # add the title			      
			        filefig5=dir_wrk+'/'+figure_labels[ii]+'_'+filefigbase+fileext[count_domain]+'.png'
			        plt.title(figure_labels[ii]+'_'+filefigbase)
			        plt.savefig(filefig5)
	                        filefig5_d=dir_wrk+'/oftheday/'+figure_labels[ii]+fileext[count_domain]+'.png'
			        plt.savefig(filefig5_d)

		        if ii==5:
			        Ucms=variables['Ucms']
			        Vcms=variables['Vcms']
                                E= Ucms**2+Vcms**2 
                                E = np.ma.masked_where(np.isnan(E), E) 
			        fig6=plt.figure()
			        cax6=mymap.pcolormesh(x,y,E,cmap=cm_oc.cm.thermal,zorder=-1)
			        # add the coastline (data from Basemap)
			        mymap.drawcoastlines()	
			        # define the tick and the grid (done for global)
                                # from south pole to nortehn pole with a resolution of 1 degree
			        myparallels=np.arange(-90,90+1,1)
			        mymeridians = np.arange(-180,180+1,reso_meridians[count_domain]) 
			        # draw the ticks labels
			        mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
			        mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
			        # draw the continent (data from Basemap)
			        mymap.fillcontinents(color='0.6',lake_color='0.6',zorder=100)
			        # draw the stations
			        mymap.plot(x_stations,y_stations,'-',color='w',zorder=1)
	                        if count_domain>0:
                                        mymap.scatter(x_stations,y_stations,s=15,color='w',zorder=1)
			        # draw the waypoint
			        #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                                # draw the glider
				mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)
                                # draw the ZEE limits
	                        mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1) 	
			        # draw the swot trajectories
        			mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
        			mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)
			        # add the colorbar
                                if (len(kinemin) >= count_domain+1 and len(kinemax) >= count_domain+1):
			                cax6.set_clim(kinemin[count_domain],kinemax[count_domain])
                                else:
                                        cax6.set_clim(kinemin[len(kinemin)-1],kinemax[len(kinemax)-1])                                           
                                cbar6=fig6.colorbar(cax6, orientation='vertical',shrink=0.9)
			        cbar6.ax.set_ylabel('Kinetic Energy [$cm^2/s^2$]')
			        # add the title			
			        filefig6=dir_wrk+'/'+figure_labels[ii]+'_'+filefigbase+fileext[count_domain]+'.png'
			        plt.title(figure_labels[ii]+'_'+filefigbase)
			        plt.savefig(filefig6)
			        filefig6_d=dir_wrk+'/oftheday/'+figure_labels[ii]+fileext[count_domain]+'.png'
			        plt.savefig(filefig6_d)

		        if ii==6:
			        touched=variables['touched']
			        fig7=plt.figure()
                                touched = np.ma.masked_where(np.isnan(touched), touched) 
			        cax7=mymap.pcolormesh(x,y,touched,cmap=cm_oc.cm.deep,zorder=-1)
			        # add the coastline (data from Basemap)
			        mymap.drawcoastlines()	
			        # define the tick and the grid (done for global)
                                # from south pole to nortehn pole with a resolution of 1 degree
			        myparallels=np.arange(-90,90+1,1)
			        mymeridians = np.arange(-180,180+1,reso_meridians[count_domain])
			        mymap.drawparallels(myparallels,labels=[1,0,0,0],fontsize=10)
			        mymap.drawmeridians(mymeridians,labels=[0,0,0,1],fontsize=10)
			        # draw the continent (data from Basemap)
			        mymap.fillcontinents(color='0.6',lake_color='0.6',zorder=100)
			        # draw the stations
			        mymap.plot(x_stations,y_stations,'-',color='w',zorder=1)
	                        if count_domain>0:
                                        mymap.scatter(x_stations,y_stations,s=15,color='w',zorder=1) 
			        # draw the waypoint
			        #mymap.plot(x_waypoint,y_waypoint,color='r',zorder=1)
                                #mymap.plot(x_waypoint[0:6],y_waypoint[0:6],'-',color='#66ffff',zorder=1)
                                #mymap.plot(x_waypoint[7:13],y_waypoint[7:13],'-',color='#66ccff',zorder=1)
                                #mymap.plot(x_waypoint[13:18],y_waypoint[13:18],'-',color='#0000ff',zorder=1)
                                # draw the glider
				mymap.plot(x_gl,y_gl,'*',color='r',markersize=0.5,zorder=1)
                                # draw the ZEE limits
	                        mymap.plot(x_zee_sp,y_zee_sp,color='w',lw=0.5,zorder=1) 
			        # draw the swot trajectories
        			mymap.plot(x_extra[0:32],y_extra[0:32],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[33:62],y_extra[33:62],'-',color=(1, 0.6, 0.6), zorder=1) 
        			mymap.plot(x_extra[63:92],y_extra[63:92],'-',color=(1, 0.6, 0.6), zorder=1)
        			mymap.plot(x_extra[93:126],y_extra[93:126],'-',color=(1, 0.6, 0.6), zorder=1)
			        # add the colorbar
			        cax7.set_clim(0,15)
			        cbar7=fig7.colorbar(cax7, orientation='vertical',shrink=0.9)
			        cbar7.ax.set_ylabel('Time from bathymetry [d]')
			        # add the title			        
			        filefig7=dir_wrk+'/'+figure_labels[ii]+'_'+filefigbase+fileext[count_domain]+'.png'
			        plt.title(figure_labels[ii]+'_'+filefigbase)
			        plt.savefig(filefig7)		
   		                filefig7_d=dir_wrk+'/oftheday/'+figure_labels[ii]+fileext[count_domain]+'.png'
			        plt.savefig(filefig7_d)

                count_domain = count_domain +1

