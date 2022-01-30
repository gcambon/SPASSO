#!/bin/sh
# Script to extract glider position and update a webpage
# The aim is to regularly execute this script with a crontab.

# Download emails + extract lat lon (creates gliderLat.txt and gliderLon.txt)
python /home/glider/realTimePosition/getPosition.py

# Build html googleMap file (gliderMap.html)
python /home/glider/realTimePosition/map.py

### copy gliderMap to the web #AR added 25/047/2019
rsync -auv /home/glider/realTimePosition/gliderMap.html /home/SPASSO.beta/Web/DEMO/Glider_web/gliderMap.html
/bin/sh /home/SPASSO.beta/Cruises/DEMO/indicepagweb_Glider.sh /home/SPASSO.beta/Web/DEMO/Glider_web/

# Goodbye message
#echo 'Glider webpage now up-to-date!'





