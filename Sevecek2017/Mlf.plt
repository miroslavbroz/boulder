#!/usr/bin/gnuplot

Mtot = 1.0

#Mlf(Q_Qstar) = 0.008e0*((Q_Qstar)**1.4e0*exp(-(Q_Qstar/3.0e0)**1.6e0)) * Mtot
Mlf(Q_Qstar) = 0.022e0*((Q_Qstar)**1.0e0*exp(-(Q_Qstar/5.0e0)**1.2e0)) * Mtot

set xl "Q/Qstar"
set yl "Mlr [Mtot]"

set xr [1e-3:1e3]
set yr [1e-4:]

set samples 1000
set logscale x
set logscale y
set zeroaxis
set key bottom

tmp=1.0; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0

p \
  "Sevecek_etal_2017_Fig4_Mlf.dat" u 1:2 w p,\
  Mlf(x) w l,\
  1.0 w l lt 0

pa -1

q


