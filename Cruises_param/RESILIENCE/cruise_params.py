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
#                  cruise_params.py
#
# Python script to define the plots limits
#
# 28/08/2018: SB changed sstmax to sstmax = [29,28.5]
# 26/07/2018: SB changed fslemax to fslemax = [0.3,0.01]
# 23/07/2018: SB: preparation for MOANA_MATY multi_campaign --> all variables parameters have a size 2 adtmin = [0.76,0.92] adtmax = [1.12,1.04]sstmin = [25,26] sstmax = [29,29] chlmin = [0.1,0.1] chlmax = [0.5,0.5] owmin = [-0.4,-0.3] owmax = [0.4,0.3] fslemin = [0,0] fslemax = [0.4,0.1]deltalonmin = [-5,-3] deltalonmax = [5,3] deltalatmin = [-3,-2] deltalatmax = [3,2] kinemin = [0,0]  kinemax = [2000,800]
# 11/07/2018: SB: preparation for MOANA_MATY: changed sst to [25;29], ssh to [0.76;1.12]
# 06/07/2018: SB: preparation for MOANA_MATY: changed sst to [20;30]
# 05/07/2018: SB: begin preparation for MOANA_MATY: changed adt to [0.5;1.3]
# 07/05/2018: SB: changed sstmin to 14.5 and sstmax to 18.5
# 07/05/2018: SB: changed chlmin to 0.1 and chlmax to 0.5
# 06/04/2018: SB: added deltalonmin=-5 deltalonmax=5 deltalatmin=-3 deltalatmax=3 
# 15/08/2018: SB: changed chlmin to 0.15 and chlmax to 1.1
# 21/02/2018: SB: changed kinemax to 2000
# 20/02/2018: SB: changed kinemax to 1500
# 16/02/2018: SB: changed kinemax to 1000, adtmin to -0.36
# 15/02/2018: SB: changed kinemax to 700
# 14/02/2018: SB: changed adt to [-0.3; 0.16]
# 13/02/2018: SB: changed chlmax to 0.4, kinemax to 500,fslemax to 0.5, ow to [-0.4;0.4] 
# 12/02/2018: SB: added kinemin (0) & kinemax (400) 
# 12/02/2018: SB: changed fslemin (0) & fslemax (0.4) 
# 12/02/2018: SB: added owmin (-0.3) & owmax (0.3) 
# 12/02/2018: SB: added chlmin (0.01) & chlmax (0.3) 
# 12/02/2018: SB: changed sstmin (12) & sstmax(16.5)
# 12/02/2018: SB: added adtmin (-0.4) & adtmax (0.2)
# 12/02/2018: SB: creation 
#
# Credits: alice.dellapenna@mio.osupytheas.fr
#          andrea.doglioli@univ-amu.fr
###############################################################
#
# Colorbar limits
# Values for [ domain_limits.py, domain_limits_zoom.py ]
#adtmin = -0.36
#adtmax = 0.16
adtmin = [ 0.5, 0.5]   #AR added 21/03/2019
adtmax = [ 1.8, 1.8]     #AR added 21/03/2019

#sstmin = 14.5#14
#sstmax = 18.5 #18
sstmin = [23,28]#14      #AR added 21/03/2019
sstmax = [30,30] #18     #AR added 21/03/2019

#chlmin = [0.1,0.1] #0.15 # ! not 0 this will be taken in log
#chlmax = [0.65,0.5] #1.1
chlmin = [0.1,0.1] #0.15 # ! not 0 this will be taken in log
chlmax = [0.45,0.45] #1.1
chl_ticks=[np.log(0.05),np.log(0.10),np.log(0.15),np.log(0.20),np.log(0.25),np.log(0.4)]
chl_ticklabels=['0.05','0.10','0.15','0.20','0.25','0.40']  # vertically oriented colorbar

owmin = [-0.4,-0.3]
owmax = [0.3,0.3]

#fslemax = [0.5,0.1]
#fslemax = [0.4,0.1]
fslemin = [0,0]
fslemax = [0.2,0.01]

deltalonmin = [-5,-3]
deltalonmax = [ 5, 3]

deltalatmin = [-3,-2]
deltalatmax = [ 3, 2]

kinemin = [0,0] 
kinemax = [1700,800]



