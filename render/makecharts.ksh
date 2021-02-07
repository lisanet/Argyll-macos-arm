#!/bin/sh

# create some gamut mapping test charts
./timage cube.tif
./timage -x cube16.tif
./timage -g 20 cubeG20.tif
./timage -g 50 cubeG50.tif
./timage -4 cmykcube.tif
./timage -t -s surface.tif
./timage -t -s -x surface16.tif
./timage -t -s -g 20 surfaceG20.tif
./timage -t -s -g 50 surfaceG50.tif
