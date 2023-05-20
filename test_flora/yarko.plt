#!/usr/bin/gnuplot

set term x11

set xl "{/Helvetica-Oblique D} [km]"
set yl "{/Symbol t} [My]"

set xr [0.8*1.0e-3:1000]
set yr [1:1e6]
set logscale xy
set xtics (\
  "1000"    1e3  0,\
  "100"     1e2  0,\
  "10"      1e1  0,\
  "1"       1e0  0,\
  "0.1"     1e-1 0,\
  "0.01"    1e-2 0,\
  "10^{-3}" 1e-3 0,\
  "10^{-4}" 1e-4 0,\
  "10^{-5}" 1e-5 0,\
  "10^{-6}" 1e-6 0,\
  "10^{-7}" 1e-7 0,\
  )
set ytics (\
  "1"    1e0 0,\
  ""     2e0 1,\
  ""     3e0 1,\
  ""     4e0 1,\
  ""     5e0 1,\
  ""     6e0 1,\
  ""     7e0 1,\
  ""     8e0 1,\
  ""     9e0 1,\
  "10"   1e1 0,\
  ""     2e1 1,\
  ""     3e1 1,\
  ""     4e1 1,\
  ""     5e1 1,\
  ""     6e1 1,\
  ""     7e1 1,\
  ""     8e1 1,\
  ""     9e1 1,\
  "100"  1e2 0,\
  ""     2e2 1,\
  ""     3e2 1,\
  ""     4e2 1,\
  ""     5e2 1,\
  ""     6e2 1,\
  ""     7e2 1,\
  ""     8e2 1,\
  ""     9e2 1,\
  "10^3" 1e3 0,\
  ""     2e3 1,\
  ""     3e3 1,\
  ""     4e3 1,\
  ""     5e3 1,\
  ""     6e3 1,\
  ""     7e3 1,\
  ""     8e3 1,\
  ""     9e3 1,\
  "10^4" 1e4 0,\
  ""     2e4 1,\
  ""     3e4 1,\
  ""     4e4 1,\
  ""     5e4 1,\
  ""     6e4 1,\
  ""     7e4 1,\
  ""     8e4 1,\
  ""     9e4 1,\
  "10^5" 1e5 0,\
  ""     2e5 1,\
  ""     3e5 1,\
  ""     4e5 1,\
  ""     5e5 1,\
  ""     6e5 1,\
  ""     7e5 1,\
  ""     8e5 1,\
  ""     9e5 1,\
  "10^6" 1e6 0,\
  ""     2e6 1,\
  ""     3e6 1,\
  ""     4e6 1,\
  ""     5e6 1,\
  ""     6e6 1,\
  ""     7e6 1,\
  ""     8e6 1,\
  ""     9e6 1,\
  "10^7" 1e7 0,\
  )
set key left width -0.0 samplen 1.5 spacing 1.05
set tics front

set arrow from 1,graph 0 rto 0,graph 1 nohead lt 0
set arrow from 0.001,graph 0 rto 0,graph 1 nohead lt 0

set lmargin 7.0
set rmargin 2.0
set bmargin 3.1
set tmargin 0.8

p \
  "yarko_dadt_EDITED.out" u 1:2 w l lw 7 dt 1 lc '#0066ff' t "MB",\
  "yarko_dadt_EDITED.out" u 1:2 w l lw 5 dt 2 lc '#ff6600' t "Flora",\
  "yarko_dadt_inner_belt_tau.out"  u 1:2 w l dt 2 lc 'gray'    t "inner belt",\
  "yarko_dadt_middle_belt_tau.out" u 1:2 w l dt 2 lc '#999999' t "middle belt",\
  "yarko_dadt_outer_belt_tau.out"  u 1:2 w l dt 2 lc 'black'   t "outer belt",\
  "yarko_bottke_tau.dat" u 1:2 w l t "{/=14 Bottke etal. (2005)}",\

pa -1

set term post eps enh color solid "Helvetica,18"
set out "yarko.eps"
set size 0.675,0.9
rep

system("eps2png yarko.eps")

q


