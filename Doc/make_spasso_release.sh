#! /bin/sh
echo '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '#%%%%                 SPASSO                  vs 0.0.0beta  %%%'
echo '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
#%%  This  file  is  free software ;  it  is distributed  in  %
#%%  the hope that  it  will  be  useful,  but  without  any  %
#%%  warranty.  You  can  redistribute  it and/or modify  it  %
#%%  under  the  terms  of  the  GNU  General Public License  %
#%%  as  published  by  the   Free  Software  Foundation  at  %
#%%       http://www.gnu.org/copyleft/gpl.html                %
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo '# this script build a .tar.gz fils containing a standard SPASSO'
echo '# directory tree and source files except data'
echo '#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
echo ' '
echo ' '
#echo '!!! WARNING : this script removes password from spasso.sh file !!!'
#echo 'press Enter to continue'
#read pipo

cruise_dir=$1
backup_date=`date +%Y%m%d`

maindir=SPASSO.beta_$backup_date

mkdir $maindir/
#mkdir $maindir/Data
#cp -rf Data/aviso_upd_global_merged_madt_uv_2010  $maindir/Data/
#cp -rf Data/aviso_upd_global_merged_madt_uv_2014  $maindir/Data/
#cp -rf Data/bathymetry  $maindir/Data/

mkdir $maindir/Cruises
if [ -n "$cruise_dir" ]; then
  #  cp -rf Cruises/$cruise_dir $maindir/Cruises/
    rsync -avz --exclude 'Processed'  /home/SPASSO.beta/Cruises/$cruise_dir $maindir/Cruises/
else
  cp -rf Cruises/DEMO $maindir/Cruises/
fi

mkdir $maindir/Scripts
cp -rf Scripts/* $maindir/Scripts/

mkdir $maindir/Doc
cp -rf Doc/* $maindir/Doc/

mkdir $maindir/Web
cp -rf Web/*.* $maindir/Web/

# copy all in maindir except Data/ and maindir
#rsync -avz --exclude 'Data' --exclude $maindir /home/SPASSO/ $maindir/

#### cp spasso.sh with password
#cruisepwd=/home/SPASSO/Cruises/*
#for f in $cruisepwd
#do
#cp $f/spasso.sh $f/spasso.sh.pwd
#done
#
#### Remove usernames and passwords form spasso.sh file
#cruises=/home/SPASSO/$maindir/Cruises/*
#for c in $cruises
#do
#sed -i -e 's/user_cls=lob_doglioli/user_cls=???/g' $c/spasso*.sh
#sed -i -e 's/pwd_cls=coral15/pwd_cls=???/g' $c/spasso*.sh
#sed -i -e 's/user_aviso=lob_doglioli/user_aviso=???/g' $c/spasso*.sh
#sed -i -e 's/pwd_aviso=coral15/pwd_aviso=???/g' $c/spasso*.sh
#sed -i -e 's/user_myocean=adoglioli/user_myocean=???/g' $c/spasso*.sh
#sed -i -e 's/pwd_myocean=Andrea01/pwd_myocean=???/g' $c/spasso*.sh
#sed -i -e 's/user_aviso=adoglioli/user_aviso=???/g' $c/spasso*.sh
#sed -i -e 's/pwd_aviso=Andrea01/pwd_aviso=???/g' $c/spasso*.sh
##    grep -rl 'user_aviso=mio_rousselet' $c/spasso*.sh | xargs sed -i 's/user_aviso=mio_rousselet/user_aviso=???/g'
##    grep -rl 'pwd_aviso=iez11r6' $c/spasso*.sh | xargs sed -i 's/pwd_aviso=iez11r6/pwd_aviso=???/g'
##    grep -rl 'user_aviso=lob_doglioli' $c/spasso*.sh | xargs sed -i 's/user_aviso=lob_doglioli/user_aviso=???/g'
##    grep -rl 'pwd_aviso=coral15' $c/spasso*.sh | xargs sed -i 's/pwd_aviso=coral15/pwd_aviso=???/g'
##    grep -rl 'user_myocean=adoglioli' $c/spasso*.sh | xargs sed -i 's/user_myocean=adoglioli/user_myocean=???/g'
##    grep -rl 'pwd_myocean=Andrea01' $c/spasso*.sh | xargs sed -i 's/pwd_myocean=Andrea01/pwd_myocean=???/g'
##    grep -rl 'user_cls=lob_doglioli' $c/spasso*.sh | xargs sed -i 's/user_cls=lob_doglioli/user_cls=???/g'
##    grep -l 'pwd_cls=coral15' $c/spasso*.sh | xargs sed -i 's/pwd_cls=coral15/pwd_cls=???/g'
#done
#
##### replace spasso.sh with password in /home/SPASSO/Cruises/
#for f in $cruisepwd
#do
#mv $f/spasso.sh.pwd $f/spasso.sh
#done


#cp *.sh $maindir/
#cp *README*.txt $maindir/

tar cvfz $maindir.tar.gz $maindir/

rm -rf $maindir

#mv $maindir.tar.gz /mnt/OSCAHR/Figures_web/SPASSO_releases/

echo '%%% The End %%%'

