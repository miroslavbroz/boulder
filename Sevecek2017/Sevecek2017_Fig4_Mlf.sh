#!/bin/sh

FILE=Sevecek2017_Fig4_Mlf

map2dat.awk $FILE.map | awk '(NF>0){ print 10**$1,10**$2; }' > $FILE.dat


