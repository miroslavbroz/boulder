#!/usr/bin/gnuplot

set xl "log D [km]"
set yl "log N(>D)"

p "Harris_etal_2015.dat" u 1:2 w lp

pa -1

