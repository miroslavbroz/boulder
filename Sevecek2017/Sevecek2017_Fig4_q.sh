#!/bin/sh

FILE=Sevecek2017_Fig4_q

map2dat.awk $FILE.map | awk '(NF>1){ print 10**$1,$2; }' > $FILE.dat


