#!/bin/sh

cd /home/SPASSO.beta/Data/GLIDER_SEA003/
`/usr/bin/smbclient -d0 //139.124.16.155/SEA003 -U glideradm%MIO2017 << EOF
recurse ON
prompt OFF
mget M*
exit EOF`

