#!/usr/bin/gnuplot

set xl "log D [km]"
set yl "log N(>D)"

p "Bottke_etal_2005.dat" u 1:2 w lp
pa -1

