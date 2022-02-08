#!/bin/sh

#subdir=/home/SPASSO/Web/PREBIOSWOT/Figures_web_oftheday/
subdir=$1
dir=$subdir
chmod a+r $dir/*.*

#cp /home/SPASSO.beta/Scripts/generate_index_Figures_oftheday.sh $dir/
cp ${cruise_path}/Process_web/generate_index_Figures_oftheday.sh $dir/
chmod a+x $dir/generate_index_Figures_oftheday.sh

cd $dir
./generate_index_Figures_oftheday.sh

/bin/rm $dir/generate_index_Figures_oftheday.sh

chmod a+rx $dir
cd -
#--------------------------

