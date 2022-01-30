#!/bin/sh

#subdir=/home/SPASSO/Web/PREBIOSWOT/Bulletin_web/
subdir=$1
dir=$subdir
chmod a+r $dir/*.*

cp /home/SPASSO.beta/Scripts/generate_index.sh $dir
chmod a+x $dir/generate_index.sh

cd $dir
./generate_index.sh

/bin/rm $dir/generate_index.sh

chmod a+rx $dir

#--------------------------

