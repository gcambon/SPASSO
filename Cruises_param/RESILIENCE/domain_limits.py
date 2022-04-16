#############################################################
##                   SPASSO                  vs 0.0.0beta  ###
##############################################################
##  This  file  is  free software ;  it  is distributed  in  #
##  the hope that  it  will  be  useful,  but  without  any  #
##  warranty.  You  can  redistribute  it and/or modify  it  #
##  under  the  terms  of  the  GNU  General Public License  #
##  as  published  by  the   Free  Software  Foundation  at  #
##       http://www.gnu.org/copyleft/gpl.html                #
##############################################################
#                  domain_limits.py
#
# Python script to define the maps domain
#
# 05/07/2018: SB: changed domain lon to [-150;-130] for MOANA_MATY]
# 16/02/2018: SB: changed domain lon to [-2:+9]
# 15/02/2018: SB: changed domain lat to [35.5; 44.5], removed colorbar limits
# 02/05/2017: LR: add parameters for various colorbars
# 30/09/2016: ADP&AD: creation 
# 
# Credits: alice.dellapenna@mio.osupytheas.fr
#          andrea.doglioli@univ-amu.fr
###############################################################
# Set zoom on the domain
#MOVE TO Cruise DIRECTORY

# FUMSECK campaign within:
#Lon=[5, 11]; # AR changed  for FUMSECK 19/03/2017
#Lat=[42, 44.7]; # AR changed  for FUMSECK 19/03/2017

# RESILIENCE campaign
#Lon=[  30,  48];
#Lat=[ -32, -21];  #=> #1

Lon=[  34,  45];
Lat=[ -24, -18];  #=> #2


# Colorbar limits #SB removed 15/02/2018
#sstmin = 15
#sstmax = 20
#chl_ticks=[np.log(0.05),np.log(0.10),np.log(0.15),np.log(0.20),np.log(0.25),np.log(0.4)]
#chl_ticklabels=['0.05','0.10','0.15','0.20','0.25','0.40']  # vertically oriented colorbar
#fslemin = 0
#fslemax = 0.4

