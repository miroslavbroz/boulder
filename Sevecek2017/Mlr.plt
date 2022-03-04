#!/usr/bin/gnuplot

Mtot = 1.0

#Mlr(Q_Qstar) = Q_Qstar < 1.0 ? (-0.5e0*(Q_Qstar-1.e0)+0.5e0)*Mtot : (-0.35e0*(Q_Qstar-1.e0)+0.5e0)*Mtot
Mlr(Q_Qstar) = Q_Qstar < 1.0 ? (-0.5e0*(Q_Qstar-1.e0)+0.5e0)*Mtot : (-0.20e0*(Q_Qstar-1.e0)+0.5e0)*Mtot

set xl "Q/Qstar"
set yl "Mlr [Mtot]"

set xr [1e-3:1e3]
set yr [1e-3:]

set samples 1000
set logscale x
set logscale y
set zeroaxis
set key bottom

tmp=1.0; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0

p \
  "Sevecek_etal_2017_Fig4_Mlr.dat" u 1:2 w p,\
  Mlr(x) w l,\
  1.0 w l lt 0

pa -1

q


