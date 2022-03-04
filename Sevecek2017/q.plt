#!/usr/bin/gnuplot

Mtot = 1.0

#qf(Q_Qstar) = -10.e0+7.e0*(Q_Qstar)**0.4e0*exp(-Q_Qstar/7.e0)
qf(Q_Qstar) = -10.e0+7.3e0*(Q_Qstar)**0.05e0*exp(-Q_Qstar/100.e0)

set xl "Q/Qstar"
set yl "Mlr [Mtot]"

set xr [1e-3:1e3]

set samples 1000
set logscale x
set zeroaxis
set key bottom

tmp=1.0; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0

p \
  "Sevecek_etal_2017_Fig4_q.dat" u 1:2 w p,\
  qf(x) w l

pa -1

q


