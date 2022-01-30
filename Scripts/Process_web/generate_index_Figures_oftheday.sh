#!/bin/sh
#
# generate_index_html
#
# Generates index.html files with a link to each file in the directory, and
# recursively for its sub-directories.
#
# Note: This script recursively calls itself, so it had better be in its
# own path!
#

INDEX=index.html

    #
    # These are the standard locations for apache's icons
    #
img_txt="<img src=/icons/text.gif>"
img_dir="<img src=/icons/folder.gif>"
img_bak="<img src=/icons/back.gif>"

echo "<html> <h3>This web directory contains the figures created each 3 hours by SPASSO (Software Package for A Satellite-based Strategy for Oceanographic cruises) treating and analyzing the near-real time data from AVISO, NASA-Oceancolor, MyOceans." > $INDEX
#echo "<H2>Directory listing</H2>" > $INDEX
echo "<hr>" >> $INDEX
echo "$img_bak <a href=../>Parent Directory</a><br>" >> $INDEX
echo "<table width=90% border=0>" >> $INDEX

#echo "<tr><td><b> filename </td><td><b> size </td><td><b> creation date (UTC) </td><td><b> rights </td></tr>" >> $INDEX

for file in `ls -t`
do
    if [ $file != index.html ] && [ $file != README ]  && [ $file != generate_index_Figures_oftheday.sh ]
    then
        if [ -d $file ]
        then
            #AD    (cd $file && generate_index_html)
            #echo "$img_dir <a href=$file/>$file/</a><br>" >> $INDEX
	    echo " "
#<tr><td>$img_dir <a href=$file/>$file/</a></td></tr>" >> $INDEX
	else
            size=`ls -lrt $file | awk '{ print $5; }'`
            sizekb=`expr $size / 1024`
            diritti=`ls -lrt $file | awk '{ print $1; }'`
            day=`ls -lrt $file | awk '{ print $7; }'`
            month=`ls -lrt $file | awk '{ print $6; }'`
            year=`ls -lrt $file | awk '{ print $8; }'`
            #echo "<tr><td>$img_txt <a href=$file>$file</a></td><td>($sizekb kilobytes)</td><td>$day $month $year</td><td>$diritti</td></tr>" >> $INDEX
           # echo "<a href=$file><img src=$file height=10%></a> $file ($sizekb kilobytes)<br>" >> $INDEX
            echo "<a href=$file><img src=$file height=10%></a> $file ($sizekb kilobytes) <i> $day $month $year </i> <br>" >> $INDEX
        fi
    fi
done

#echo "<hr>" >> $INDEX
echo "</table>" >> $INDEX
echo "<hr>" >> $INDEX

if [ -f README ]
then
    echo "<pre>" >> $INDEX
    cat README >> $INDEX
    echo "</pre>" >> $INDEX
fi

chmod a+r,g+w $INDEX
#AD# chgrp mercury $INDEX

