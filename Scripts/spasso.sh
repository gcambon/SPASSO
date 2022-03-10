#!/bin/bash
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%                 SPASSO                  vs 0.0.0beta  %%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%  This  file  is  free software ;  it  is distributed  in  %
#%%  the hope that  it  will  be  useful,  but  without  any  %
#%%  warranty.  You  can  redistribute  it and/or modify  it  %
#%%  under  the  terms  of  the  GNU  General Public License  %
#%%  as  published  by  the   Free  Software  Foundation  at  %
#%%       http://www.gnu.org/copyleft/gpl.html                %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#F.d'Ovidio, F.Nencioli, L.Rousselet, A.Doglioli, Alice Della Penna, Stephanie Barrillon, Anaïs Ricout.

###############################################################
echo
echo
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '%%                        SPASSO                         %%'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '%%                                                       %%'
echo '%%                                                       %%'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

echo
echo
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '%%                    Python ENV activation              %%'
source ~/.bash_profile
conda activate mypy37_base
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

#set -e
#set -x

############# 1st step: DEFINE THE PATHS ############# 
# Home
export main_path=/home/gcambon/HCONFIGS_SPASSO #echo $main_path to display it

# Current cruise
export cruise_path=$main_path/RESILIENCE 

# Data
export dir_data=$cruise_path/Data

#=============================
# Dataset 

##### ALT #####
#ALTI_product='nrt.cmems-du.eu/Core/SEALEVEL_EUR_PHY_L4_NRT_OBSERVATIONS_008_060/dataset-duacs-nrt-europe-merged-allsat-phy-l4'
ALTI_product='nrt.cmems-du.eu/Core/SEALEVEL_GLO_PHY_L4_NRT_OBSERVATIONS_008_046/dataset-duacs-nrt-global-merged-allsat-phy-l4'

#### CHL #####
#-- CHL_L4
CHL_L4_product='nrt.cmems-du.eu/Core/OCEANCOLOUR_GLO_CHL_L4_NRT_OBSERVATIONS_009_033/dataset-oc-glo-bio-multi-l4-chl_interpolated_4km_daily-rt'

#-- CHL_L3
CHL_L3_MULTI_med_product='nrt.cmems-du.eu/Core/OCEANCOLOUR_MED_CHL_L3_NRT_OBSERVATIONS_009_040/dataset-oc-med-chl-multi-l3-chl_1km_daily-rt-v02'   
CHL_L3_OLCI_a_med_product='nrt.cmems-du.eu/Core/OCEANCOLOUR_MED_CHL_L3_NRT_OBSERVATIONS_009_040/dataset-oc-med-chl-olci_a-l3-chl_1km_daily-rt-v02'   

#CHL_L3_MULTI_glob_product='nrt.cmems-du.eu/Core/OCEANCOLOUR_GLO_CHL_L3_NRT_OBSERVATIONS_009_032/dataset-oc-glo-bio-multi-l3-pft_4km_daily-rt'
CHL_L3_MULTI_glob_product='nrt.cmems-du.eu/Core/OCEANCOLOUR_GLO_CHL_L3_NRT_OBSERVATIONS_009_032/dataset-oc-glo-bio-multi-l3-chl_4km_daily-rt'
CHL_L3_OLCI_a_glob_product='nrt.cmems-du.eu/Core/OCEANCOLOUR_GLO_CHL_L3_NRT_OBSERVATIONS_009_032/dataset-oc-glo-chl-olci_a-l3-av_4km_daily-rt-v02'

#### SST ####
SST_L4_product='nrt.cmems-du.eu/Core/SST_GLO_SST_L4_NRT_OBSERVATIONS_010_005/METOFFICE-GLO-SST-L4-NRT-OBS-GMPE-V3' 
SST_L3_product='nrt.cmems-du.eu/Core/SST_GLO_SST_L3S_NRT_OBSERVATIONS_010_010/IFREMER-GLOB-SST-L3-NRT-OBS_FULL_TIME_SERIE'
#SST_L3_product='nrt.cmems-du.eu/Core/SST_MED_SST_L3S_NRT_OBSERVATIONS_010_012/SST_MED_SST_L3S_NRT_OBSERVATIONS_010_012_a' 

SST_L4_JPL_product='podaac-ftp.jpl.nasa.gov/OceanTemperature/ghrsst/data/L4/GLOB/JPL_OUROCEAN/G1SST'
#=============================

#=============================
# Type of product for lagrangian diagnostics
lagrangian_type='nrt_cmems' 

# To differentiate the paths between NRT and DT products
product_type='nrt' 

dir_ALTI=$dir_data/ALTI/$ALTI_product/ 
#
dir_CHL_L4=$dir_data/CHL/$CHL_L4_product/
dir_CHL_L3_MULTI_glob=$dir_data/CHL/$CHL_L3_MULTI_glob_product/ 
dir_CHL_L3_OLCI_a_glob=$dir_data/CHL/$CHL_L3_OLCI_a_glob_product/
#
dir_CHL_L3_MULTI_med=$dir_data/CHL/$CHL_L3_MULTI_med_product/ 
dir_CHL_L3_OLCI_a_med=$dir_data/CHL/$CHL_L3_OLCI_a_med_product/ 
# 
dir_SST_L4=$dir_data/SST/$SST_L4_product/ 
dir_SST_L3=$dir_data/SST/$SST_L3_product/ 
dir_SST_L4_JPL=$dir_data/SST/$SST_L4_JPL_product/ 

# Work
dir_wrk=$cruise_path/Wrk

# Scripts
dir_scripts=$cruise_path/Process_satellite
dir_lagrang=$cruise_path/Lagrangian_package
dir_lamta=$cruise_path/Lagrangian_package/lamta.dev/

# Processed data
dir_PROC=$cruise_path/Processed

# bulletin
dir_BUL=$cruise_path/Bulletin

# Figures
dir_FIG=$cruise_path/Figures

# Logs dir
logs_dir=$cruise_path/Logs

# Web 
web_path=$cruise_path/Web
dir_FIG_web=$web_path/Figures_web
dir_FIG_web_oftheday=$web_path/Figures_web_oftheday
dir_PROC_web=$web_path/Processed_web
dir_BUL_web=$web_path/Bulletin_web
#
web_path_remote=yy@xx.fr:XX/public_html/SPASSO/Web
dir_FIG_web_remote=$web_path_remote/Figures_web
dir_FIG_web_oftheday_remote=$web_path_remote/Figures_web_oftheday
dir_PROC_web_remote=$web_path_remote/Processed_web
dir_BUL_web_remote=$web_path_remote/Bulletin_web

############# 2d step: CHOICE OF DATA [Y/N] #############
ALTI_YESNO=Y          
user_alti=XXX 
pwd_alti=XXX 
ALTI_data=( phy ) 
#===================================================
CHL_L4_YESNO=Y
user_chl_l4=XXX
pwd_chl_l4=XXX 
CHL_data=( CHL ) 

CHL_L3_glob_YESNO=Y        
user_chl_l3=XXX 
pwd_chl_l3=XXX
CHL_data=( CHL )

CHL_L3_med_YESNO=N
user_chl_l3=XXX 
pwd_chl_l3=XXX
CHL_data=( CHL ) 
#===================================================
SST_L4_YESNO=Y
user_sst_l4=XXX 
pwd_sst_l4=XXX 
SST_data=( sst ) 

SST_L3_YESNO=Y         
user_sst_l3=XXX 
pwd_sst_l3=XXX 
SST_data=( sst ) 

SST_L4_JPL_YESNO=N   
SST_data=( sst ) 
#===================================================
# Extra Options [Y/N]
LAGRANGIAN_YESNO=N
WEB_YESNO=Y
MAILING_YESNO=N

DOWNLOAD_DATA_YESNO=Y
PROCESS_DATA_YESNO=Y

COMPRESS_DATANC_YESNO=N

COMPRESS_DATAMAT_YESNO=N
# N is for plot debuging
# Y is for prod

############# 3rd step: DATES OF FILES NAME OF DATASET #############
# orig ALTI 0day ago arrive des données 0daysago ALT CMEMS L4 00h00 
date_alti=`date +%Y%m%d --date='0day ago'`   
year_alti=`date +%Y --date='0day ago'`       
month_alti=`date +%m --date='0day ago'`  

# orig SST_L4 arrive des données -1daysago SST CMEMS L4& L3S
# a midi-2days ago OK
date_SST_L4=`date +%Y%m%d --date='2day ago'`  
year_SST_L4=`date +%Y --date='2day ago'`	 
month_SST_L4=`date +%m --date='2day ago'`

# orig SST_L3 arrive des données -1daysago SST CMEMS L4& L3S
# a midi -2days ago OK
date_SST_L3=`date +%Y%m%d --date='2day ago'`  
year_SST_L3=`date +%Y --date='2day ago'`	 
month_SST_L3=`date +%m --date='2day ago'`

# orig SST_JPL_L4
year_SST_L4_JPL=`date +%Y --date='2day ago'`
day_SST_L4_JPL=`date +%j --date='2day ago'`  
date_SST_L4_JPL=`date +%Y%m%d --date='2day ago'`

# orig CHL_L4 1day ago
# a midi -2days ago OK
date_CHL_L4=`date +%Y%m%d --date='2day ago'`  
year_CHL_L4=`date +%Y --date='2day ago'`     
month_CHL_L4=`date +%m --date='2day ago'`

# orig CHL_L3 1day ago
# a midi -2days ago OK
date_CHL_L3=`date +%Y%m%d --date='2day ago'`  
year_CHL_L3=`date +%Y --date='2day ago'`     
month_CHL_L3=`date +%m --date='2day ago'`    

date_day=`date +%D --date='1day ago'` 

if [ "$DOWNLOAD_DATA_YESNO" == Y ]; then
    ############# CLEAN WRK DIRECTORY #############
    cd $dir_wrk
    rm -rf $dir_wrk/*
    mkdir $dir_wrk/oftheday 
    cd -
    
    ############# 4th step: START DOWNLOADING DATA #############
    cd  $dir_data
    ################################################
    # Download NRT (Near-real-time) SSH and velocity data
    ################################################
    if [ "$ALTI_YESNO" == Y ]; then \
	echo
	echo '------------------ DOWNLOAD ALTI FILES in Data/ and make a copy in Wrk/-------------------'
	for data in ${ALTI_data[@]}  
	do
	    wget -r --mirror -nd --directory-prefix=$dir_data/ALTI/$ALTI_product/  -nv --no-proxy --no-proxy --user=$user_alti --password=$pwd_alti ftp://$ALTI_product/$year_alti/$month_alti/*_allsat_${data}_*_${date_alti}_*.nc  	
	    cp $dir_data/ALTI/$ALTI_product/*_allsat_${data}_*_${date_alti}_*.nc $dir_wrk/ 
	done
    fi  
    ################################################
    # Download NRT (Near-real-time) Oceancolor data
    ################################################
    if [ "$CHL_L4_YESNO" == Y ]; then 
	echo '------------------ DOWNLOAD CHL L4 FILES in in Data/ and make a copy in Wrk/-------------------'
	echo
	wget -r --mirror -nd --directory-prefix=$dir_data/CHL/$CHL_L4_product/  -nv --no-proxy --user=$user_chl_l4 --password=$pwd_chl_l4  ftp://$CHL_L4_product/${year_CHL_L4}/${month_CHL_L4}/${date_CHL_L4}_d-ACRI-L4-CHL-MULTI_4KM-GLO-NRT.nc
	cp $dir_data/CHL/$CHL_L4_product/${date_CHL_L4}_d-ACRI-L4-CHL-MULTI_4KM-GLO-NRT.nc $dir_wrk/
	
    fi
    if [ "$CHL_L3_glob_YESNO" == Y ]; then 
	echo '------------------ DOWNLOAD CHL L3 glob FILES in Data/ and make a copy in Wrk/-------------------'
	echo
	wget -r --mirror -nd --directory-prefix=$dir_data/CHL/$CHL_L3_MULTI_glob_product/  -nv --no-proxy --user=$user_chl_l3 --password=$pwd_chl_l3  ftp://$CHL_L3_MULTI_glob_product/${year_CHL_L3}/${month_CHL_L3}/${date_CHL_L3}_d-ACRI-L3-CHL-MULTI_4KM-GLO-NRT.nc 
	cp $dir_data/CHL/$CHL_L3_MULTI_glob_product/${date_CHL_L3}_d-ACRI-L3-CHL-MULTI_4KM-GLO-NRT.nc $dir_wrk/
    fi
    if [ "$CHL_L3_med_YESNO" == Y ]; then 
	echo '------------------ DOWNLOAD CHL L3 med FILES in Data/ and make a copy in Wrk/-------------------'
	echo
	wget -r --mirror -nd --directory-prefix=$dir_data/CHL/$CHL_L3_MULTI_med_product/  -nv --no-proxy --user=$user_chl_l3 --password=$pwd_chl_l3  ftp://$CHL_L3_MULTI_med_product/${year_CHL_L3}/${month_CHL_L3}/${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-NRT-v02.nc
	cp $dir_data/CHL/$CHL_L3_MULTI_med_product/${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-NRT-v02.nc $dir_wrk/ 
    fi

    ################################################
    # Download NRT (Near-real-time) Sea_Surface_Temperature data
    ################################################
    if [ "$SST_L4_YESNO" == Y ]; then 
	echo '------------------ DOWNLOAD SST L4 FILES in Data/ and make a copy in Wrk/-------------------'
	echo 
	wget -r --mirror -nd --directory-prefix=$dir_data/SST/$SST_L4_product/  -nv --no-proxy --user=$user_sst_l4 --password=$pwd_sst_l4 ftp://$SST_L4_product/${year_SST_L4}/${month_SST_L4}/${date_SST_L4}120000-UKMO-L4_GHRSST-SSTfnd-GMPE-GLOB-v03.0-fv03.0.nc 
	cp $dir_data/SST/$SST_L4_product/${date_SST_L4}120000-UKMO-L4_GHRSST-SSTfnd-GMPE-GLOB-v03.0-fv03.0.nc $dir_wrk/ 
    fi
    
    if [ "$SST_L3_YESNO" == Y ]; then 
	echo '------------------ DOWNLOAD SST L3 FILES in Data/ and make a copy in Wrk/-------------------'
	echo
	wget -r --mirror -nd --directory-prefix=$dir_data/SST/$SST_L3_product/  -nv --no-proxy --user=$user_sst_l3 --password=$pwd_sst_l3 ftp://$SST_L3_product/${year_SST_L3}/${month_SST_L4}/${date_SST_L3}-IFR-L3C_GHRSST-SSTsubskin-ODYSSEA-GLOB_010_adjusted-v2.0-fv1.0.nc
	cp $dir_data/SST/$SST_L3_product/${date_SST_L3}-IFR-L3C_GHRSST-SSTsubskin-ODYSSEA-GLOB_010_adjusted-v2.0-fv1.0.nc $dir_wrk/ 
    fi

    if [ "$SST_L4_JPL_YESNO" == Y ]; then 
	echo '------------------ DOWNLOAD SST L4 JPL FILES in Data/ and make a copy in Wrk/-------------------'
	echo
	wget -r --mirror -nd --directory-prefix=$dir_data/SST/$SST_L4_JPL_product/  -nv ftp://$SST_L4_JPL_product/${year_SST_L4_JPL}/${day_SST_L4_JPL}/${date_SST_L4_JPL}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc.bz2
	bunzip2 -k $dir_data/SST/$SST_L4_JPL_product/${date_SST_L4_JPL}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc.bz2
	cp $dir_data/SST/$SST_L4_JPL_product/${date_SST_L4_JPL}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc $dir_wrk/ 
    fi
    cd -
fi  # DOWNLOAD_DATA_YESNO

if [ "$PROCESS_DATA_YESNO" == Y ]; then
    ############# 5th step: Launch Python and Octave scripts #############
    echo
    echo '-------------------SCRIPTS EXECUTION IN PYTHON/OCTAVE--------------------'
    cd $dir_scripts
    # Process ALTI DATA 
    if [ "$ALTI_YESNO" == Y ]; then     
	for data in ${ALTI_data[@]}         
	do
	    echo ' Process ALTI '${data}   
	    python -W ignore $dir_scripts/MAIN_ALTI.py $dir_wrk $cruise_path  
	    python -W ignore $dir_scripts/plot_ALTI.py $dir_wrk $cruise_path   
	done
    fi
    
    # Process CHL DATA 
    if [ "$CHL_L4_YESNO" == Y ] || [ "$CHL_L3_med_YESNO" == Y ]; then
	for data in ${CHL_data[@]}			
	do 						
	    echo ' Process CHL '${data}  
	    python -W ignore $dir_scripts/MAIN_OCEANCOLOR.py $dir_wrk $cruise_path  
	    python -W ignore $dir_scripts/plot_OCEANCOLOR.py $dir_wrk $cruise_path		  
	done
    fi
    
    # Process SST DATA 		 			
    if [ "$SST_L3_YESNO" == Y ] || [ "$SST_L4_YESNO" == Y ] ; then			
	for data in ${SST_data[@]}			
	do					
	    echo ' Process SST '${data}   				
	    python -W ignore $dir_scripts/MAIN_SST.py $dir_wrk $cruise_path 		 
	    python -W ignore $dir_scripts/plot_SST.py $dir_wrk $cruise_path		
	done
    fi
    echo '-------------------END SCRIPT/EXIT PYTHON/OCTAVE--------------------'
    cd -
fi  #end PROCESS_DATA_YESNO


if [ "$WEB_YESNO" == Y ]; then
    # Duplicate plot copy to the web site to have the first figures sooner 
    echo
    echo '---------------PUTTING FIGURES ON A WEBSITE-----------------'
    
    rsync -auv --include '*.png' --exclude '*' $dir_wrk/ $dir_FIG_web/
    /bin/sh $cruise_path/Process_web/indicepagweb.sh $dir_FIG_web/

    rsync -auv --include '*.png' --exclude '*' $dir_wrk/oftheday/ $dir_FIG_web_oftheday/
    /bin/sh $cruise_path/Process_web/indicepagweb_Figures_oftheday.sh $dir_FIG_web_oftheday/
    
    rsync -auv $dir_PROC/ $dir_PROC_web/
    
    echo '---------------END PUTTING ON A WEBSITE-----------------'
fi


if [ "$LAGRANGIAN_YESNO" == Y ]; then
    # start Lagrangian analysis
    echo
    echo '----------------- START LAGRANGIAN ANALYSIS ---------------------------'
    cd $dir_ALTI 
    octave -q $dir_lagrang/MAIN_Lagrangian.m $dir_lagrang $dir_ALTI $dir_lamta $cruise_path $date_alti $lagrangian_type	
    python -W ignore $dir_lagrang/all_lagrangian_figs.py $dir_wrk $cruise_path
    echo '----------------- END LAGRANGIAN ANALYSIS ----------------------------'
    cd -
fi

if [ "$COMPRESS_DATANC_YESNO" == Y ]; then
    ############# 6th step: Compressing and coying phases #############
    echo
    echo '--------COMPRESS ORIGINAL DATA (SST and CHL ) FILES in DATA/ ------------------'
    cd $dir_data   

    # SST data 
    if [ "$SST_L4_YESNO" == Y ]; then 	  
	cd $dir_SST_L4 
	/bin/gzip -f ${date_SST_L4}120000-UKMO-L4_GHRSST-SSTfnd-GMPE-GLOB-v03.0-fv03.0.nc 
    fi
    
    if [ "$SST_L3_YESNO" == Y ]; then 	
	cd $dir_SST_L3
	/bin/gzip -f ${date_SST_L3}-IFR-L3C_GHRSST-SSTsubskin-ODYSSEA-GLOB_010_adjusted-v2.0-fv1.0.nc
    fi

    if [ "$SST_L4_JPL_YESNO" == Y ]; then 
	cd $dir_SST_L4_JPL 
	/bin/gzip -f ${date_SST_L4_JPL}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.nc 
    fi

    # CHL data
    if [ "$CHL_L4_YESNO" == Y ]; then 	
	cd $dir_CHL_L4
	/bin/gzip -f -v ${date_CHL_L4}_d-ACRI-L4-CHL-MULTI_4KM-GLO-NRT.nc
    fi
    if [ "$CHL_L3_glob_YESNO" == Y ]; then 	
	cd $dir_CHL_L3_MULTI_glob 
	/bin/gzip  -f ${date_CHL_L3}_d-ACRI-L3-CHL-MULTI_4KM-GLO-NRT-v02.nc 
	
	#cd $dir_CHL_L3_OLCI_a_glob  
	#/bin/gzip  -f ${date_CHL_L3}_d-ACRI-L3-CHL-AV_Oa_4KM-GLO-NRT-v02.nc
    fi

    if [ "$CHL_L3_med_YESNO" == Y ]; then 
	cd $dir_CHL_L3_MULTI_med 
	/bin/gzip  -f ${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-NRT-v02.nc 	
	#	cd $dir_CHL_L3_OLCI_a_med 
	#	/bin/gzip  -f ${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4AD4_Oa_1KM-MED-NRT-v02.nc 
    fi
    cd -
fi

if [ "$COMPRESS_DATAMAT_YESNO" == Y ]; then
    echo
    echo '--------COMPRESS AND COPY .MAT DATA FILES FROM WRK/ to PROCESSED/ ------------------'
    cd $dir_wrk
    #### ALTI DATA
    if [ "$ALTI_YESNO" == Y ]; then 
	/bin/gzip *_allsat_phy_*.mat
	cp -v *_allsat_phy_*.mat.gz $dir_PROC
	
	if [ "$LAGRANGIAN_YESNO" == Y ]; then 
	    gzip ${date_alti}_*_d0.mat 					
	    cp -v ${date_alti}_*_d0.mat.gz $dir_PROC 	
	    
	    gzip ${date_alti}_*_d0_lambda_only.mat 					
	    cp -v ${date_alti}_*_d0_lambda_only.mat.gz $dir_PROC
	fi
    fi
    #
    #### SST DATA
    if [ "$SST_L4_YESNO" == Y ]; then 
	/bin/gzip ${date_SST_L4}120000-UKMO-L4_GHRSST-SSTfnd-GMPE-GLOB-v03.0-fv03.0.mat  
	cp -v ${date_SST_L4}120000-UKMO-L4_GHRSST-SSTfnd-GMPE-GLOB-v03.0-fv03.0.mat.gz $dir_PROC 
    fi
    #
    if [ "$SST_L3_YESNO" == Y ]; then 
	/bin/gzip ${date_SST_L3}-IFR-L3C_GHRSST-SSTsubskin-ODYSSEA-GLOB_010_adjusted-v2.0-fv1.0.mat
	cp -v ${date_SST_L3}-IFR-L3C_GHRSST-SSTsubskin-ODYSSEA-GLOB_010_adjusted-v2.0-fv1.0.mat.gz $dir_PROC 
    fi
    #
    if [ "$SST_L4_JPL_YESNO" == Y ]; then 
	/bin/gzip  ${date_SST_L4_JPL}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.mat  
	cp -v ${date_SST_L4_JPL}-JPL_OUROCEAN-L4UHfnd-GLOB-v01-fv01_0-G1SST.mat.gz $dir_PROC
    fi
    
    #### CHL DATA
    if [ "$CHL_L4_YESNO" == Y ]; then     
	/bin/gzip  ${date_CHL_L4}_d-ACRI-L4-CHL-MULTI_4KM-GLO-NRT.mat
	cp -v ${date_CHL_L4}_d-ACRI-L4-CHL-MULTI_4KM-GLO-NRT.mat.gz $dir_PROC 
    fi
    if [ "$CHL_L3_glob_YESNO" == Y ]; then 
	/bin/gzip  ${date_CHL_L3}_d-ACRI-L3-CHL-MULTI_4KM-GLO-NRT.mat
	cp -v ${date_CHL_L3}_d-ACRI-L3-CHL-MULTI_4KM-GLO-NRT.mat.gz $dir_PROC 
	
	#/bin/gzip  ${date_CHL_L3}_d-ACRI-L3-CHL-AV_Oa_4KM-GLO-NRT-v02.mat 
	#cp -v ${date_CHL_L3}_d-ACRI-L3-CHL-AV_Oa_4KM-GLO-NRT-v02.mat.gz $dir_PROC 
    fi
    
    if [ "$CHL_L3_med_YESNO" == Y ]; then 
	/bin/gzip  ${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-NRT-v02.mat 
	cp -v ${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4AD4_MULTI_1KM-MED-NRT-v02.mat.gz $dir_PROC
    fi
    if [ "$CHL_L3_med_YESNO" == Y ]; then
	/bin/gzip  ${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4Ad4_Oa_1KM-MED-NRT-v02.mat 
	cp -v ${date_CHL_L3}_d-OC_CNR-L3-CHL-MedOC4Ad4_Oa_1KM-MED-NRT-v02.mat.gz $dir_PROC 
    fi
    cd -
fi  #end COMPRESS_DATAMAT_YESNO

echo
echo '------------------ COPY FILES FROM WRK/ to FIGURES/ ------------------'
cd $dir_wrk
cp *.png $dir_FIG   
tar -zcvf spasso_figures_${date_alti}.tar.gz *.png
cp spasso_figures_${date_alti}.tar.gz  $dir_FIG_web
cd -

# 7th step: Putting figures and data on the website #############
if [ "$WEB_YESNO" == Y ]; then
# #ORIG
#     echo
#     echo '---------------PUTTING FIGURES ON A WEBSITE-----------------'

#     rsync -auv --include '*.png' --exclude '*' $dir_wrk/ $dir_FIG_web/
#     /bin/sh $cruise_path/Process_web/indicepagweb.sh $dir_FIG_web/
#     #=> of the day
#     rsync -auv --include '*.png' --exclude '*' $dir_wrk/oftheday/ $dir_FIG_web_oftheday/
#     /bin/sh $cruise_path/Process_web/indicepagweb_Figures_oftheday.sh $dir_FIG_web_oftheday/

#     # syncrhonization for processed data on the web
#     rsync -auv $dir_PROC/ $dir_PROC_web/
#     /bin/sh $cruise_path/Process_web/indicepagweb_Processed.sh $dir_PROC_web/
    
#     # syncrhonization for bulletion on the web 
#     rsync -auv $dir_BUL/ $dir_BUL_web/
#     /bin/sh $cruise_path/Process_web/indicepagweb_BULLETIN.sh $dir_BUL_web/
    
#     # # syncrhonize gliderMap on the web
#     # rsync -auv /home/glider/realTimePosition/gliderMap.html $web_path/Glider_web/gliderMap.html
#     # /bin/sh $cruise_path/Process_web/indicepagweb_Glider.sh $web_path/Glider_web/
#     echo '---------------END PUTTING ON A WEBSITE-----------------'
# fi
#ON MY WEBSITE
    echo
    echo '---------------PUTTING FIGURES ON A WEBSITE-----------------'
    # local rsync
    RSYNC_CMD="rsync -auv"
    ${RSYNC_CMD} --include '*.png' --exclude '*' $dir_wrk/ $dir_FIG_web/
    /bin/sh $cruise_path/Process_web/indicepagweb.sh $dir_FIG_web/
    #=> of the day
    ${RSYNC_CMD} --include '*.png' --exclude '*' $dir_wrk/oftheday/ $dir_FIG_web_oftheday/
    /bin/sh $cruise_path/Process_web/indicepagweb_Figures_oftheday.sh $dir_FIG_web_oftheday/

    # syncrhonization for processed data on the web
    ${RSYNC_CMD} $dir_PROC/ $dir_PROC_web/
    /bin/sh $cruise_path/Process_web/indicepagweb_Processed.sh $dir_PROC_web/
    
    # syncrhonization for bulletion on the web 
    ${RSYNC_CMD} $dir_BUL/ $dir_BUL_web/
    /bin/sh $cruise_path/Process_web/indicepagweb_BULLETIN.sh $dir_BUL_web/
    
    # # syncrhonize gliderMap on the web
    # ${RSYNC_CMD} /home/glider/realTimePosition/gliderMap.html $web_path/Glider_web/gliderMap.html
    # /bin/sh $cruise_path/Process_web/indicepagweb_Glider.sh $web_path/Glider_web/
    echo '---------------END LOCAL RSYNC        -----------------'

    
    RSYNC_REMOTE_CMD="rsync -e ssh -av"
    ${RSYNC_REMOTE_CMD} $dir_FIG_web/ $dir_FIG_web_remote

    ${RSYNC_REMOTE_CMD} $dir_FIG_web_oftheday/ $dir_FIG_web_oftheday_remote

    ${RSYNC_REMOTE_CMD} $dir_PROC_web/ $dir_PROC_web_remote

    ${RSYNC_REMOTE_CMD} $dir_BUL_web/ $dir_BUL_web_remote
    
    echo '---------------END PUTTING ON A WEBSITE-----------------'
fi
    
############# 8th step:EMAILING #############
if [ "$MAILING_YESNO" == Y ]; then
    echo
    echo '-------------------EMAILING--------------------'
    cd $dir_wrk

    mailing_list="xx@xx,yy@yy"
    
    #tar -zcvf spasso_DEMO_Lagr_figures_${date_ALTI}.tar.gz *adv*png OW*png FSLE*.png *vel*.png
    tar -zcvf spasso_RESILIENCE_figures_${date_ALTI}.tar.gz nrt_*.png *SST*.png *CHL*.png 
    
    #echo "Please, find attached the SPASSO Lagrangian figures of the day for the DEMO cruise." | mutt -s "[DEMO]: SPASSO Lagrangian Figures" $mailing_list -a spasso_DEMO_Lagr_figures_${date_ALTI}.tar.gz
    
    echo "Please, find attached the SPASSO figures of the day for the RESILIENCE cruise." | mutt -s "[RESILIENCE]: SPASSO Figures" $mailing_list -a spasso_RESILIENCE_figures_${date_ALTI}.tar.gz
    
    #echo "Please, find attached the SEA003 Glider positions until today for the DEMO cruise." | mutt -s "[DEMO]: SEA003 Glider positions" $mailing_list -a /home/glider/realTimePosition/gliderLon.txt /home/glider/realTimePosition/gliderLat.txt

    echo '---------------END EMAILING     -----------------'
    cd -
fi

echo
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '%%                    END SPASSO                         %%'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '%%                                                       %%'
echo '%%                                                       %%'
echo '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo

############# END #############
