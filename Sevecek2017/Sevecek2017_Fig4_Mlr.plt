#!/usr/bin/gnuplot

vkms = 5.0
Mtot = 1.0

Mlr1(Q_Qstar) = Mtot/(1.e0+(0.6e0+56.e0*exp(-1.e0*vkms)) * Q_Qstar**(0.8e0+8.e0*exp(-0.7e0*vkms)))
Mlr2(Q_Qstar) = Q_Qstar < 1.0 ? (-0.5e0*(Q_Qstar-1.e0)+0.5e0)*Mtot : (-0.35e0*(Q_Qstar-1.e0)+0.5e0)*Mtot

set xl "Q_Qstar"
set yl "Mlf [Mtot]"

set xr [1e-2:20]
set yr [0.2:1.0]
#set xr [1e-3:1e3]
#set yr [1e-3:1.0]

set logscale xy
set samples 1000

p "Sevecek2017_Fig4_Mlr.dat" w p,\
  Mlr1(x),\
  Mlr2(x)

pa -1



