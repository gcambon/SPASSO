#!/bin/sh

#subdir=/home/SPASSO/Web/PREBIOSWOT/Figures_web/
subdir=$1
echo $subdir
dir=$subdir
chmod a+r $dir/*.*

cp /home/SPASSO.beta/Scripts/generate_index.sh $dir/
chmod a+x $dir/generate_index.sh

cd $dir
./generate_index.sh
#/bin/sh $dir/generate_index.sh

/bin/rm $dir/generate_index.sh

chmod a+rx $dir
cd -

#--------------------------

