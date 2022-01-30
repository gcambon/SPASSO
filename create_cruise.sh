#!/bin/bash

set -x
set -e

#===========================================================
#user definition
SPASSO_SOURCES=/home/gcambon/MYGITHUB/SPASSO/

HOMEDIR_SPASSO=/home/gcambon/HCONFIGS_SPASSO/
OUTDIR_SPASSO=/local/tmp/3/gcambon/CONFIGS_SPASSO/
CRUISE_NAME=DEMO2
#===========================================================

# Create HOME for specific cruises
mkdir -p ${HOMEDIR_SPASSO}/$CRUISE_NAME
cp -Rf create_cruise.sh ${HOMEDIR_SPASSO}/$CRUISE_NAME/create_cruise.sh.back

# Create OUTDIR for specific cruise
CRUISE_OUT_DIR=${OUTDIR_SPASSO}/${CRUISE_NAME}
mkdir -p $CRUISE_OUT_DIR

# Create the subdirectory
# P
mkdir -p $CRUISE_OUT_DIR/Wrk
mkdir -p $CRUISE_OUT_DIR/Data
mkdir -p $CRUISE_OUT_DIR/Processed
mkdir -p $CRUISE_OUT_DIR/Bulletin
mkdir -p $CRUISE_OUT_DIR/Figures
mkdir -p $CRUISE_OUT_DIR/Logs
#  Web
mkdir -p $CRUISE_OUT_DIR/Web/Processed_web
mkdir -p $CRUISE_OUT_DIR/Web/Bulletin_web
mkdir -p $CRUISE_OUT_DIR/Web/Glider_web
mkdir -p $CRUISE_OUT_DIR/Web/Figures_web
mkdir -p $CRUISE_OUT_DIR/Web/Figures_web_oftheday

# Create links
ln -sf $CRUISE_OUT_DIR/Wrk       ${HOMEDIR_SPASSO}/$CRUISE_NAME
ln -sf $CRUISE_OUT_DIR/Data      ${HOMEDIR_SPASSO}/$CRUISE_NAME
ln -sf $CRUISE_OUT_DIR/Processed ${HOMEDIR_SPASSO}/$CRUISE_NAME
ln -sf $CRUISE_OUT_DIR/Bulletin  ${HOMEDIR_SPASSO}/$CRUISE_NAME
ln -sf $CRUISE_OUT_DIR/Figures   ${HOMEDIR_SPASSO}/$CRUISE_NAME
ln -sf $CRUISE_OUT_DIR/Logs      ${HOMEDIR_SPASSO}/$CRUISE_NAME

# Create links for web
ln -sf $CRUISE_OUT_DIR/Web       ${HOMEDIR_SPASSO}/$CRUISE_NAME

#======================================================================
# Copy the scripts file
cp -Rf $SPASSO_SOURCES/Scripts/spasso.sh          ${HOMEDIR_SPASSO}/$CRUISE_NAME
# - for download and process of satelitte data
cp -Rf $SPASSO_SOURCES/Scripts/Process_satellite  ${HOMEDIR_SPASSO}/$CRUISE_NAME
cp -Rf $SPASSO_SOURCES/Scripts/Lagrangian_package ${HOMEDIR_SPASSO}/$CRUISE_NAME

# - for the web processing
cp -Rf $SPASSO_SOURCES/Scripts/Process_web        ${HOMEDIR_SPASSO}/$CRUISE_NAME

# - for cruise parameter and domaine
cp -Rf $SPASSO_SOURCES/Cruises_param              ${HOMEDIR_SPASSO}/$CRUISE_NAME

#======================================================================
# Create the dir for land textfile
cp -Rf $SPASSO_SOURCES/Cruises_param              ${HOMEDIR_SPASSO}/$CRUISE_NAME
