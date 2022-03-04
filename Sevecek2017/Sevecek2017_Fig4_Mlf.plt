#!/usr/bin/gnuplot

vkms = 5.0
Mtot = 1.0

Mlf1(Q_Qstar) = Mtot/(0.24e0*vkms**3 * (Q_Qstar)**(-0.6e0-2.e0*exp(-0.3e0*vkms)) + exp(-0.3e0*vkms)*Q_Qstar + 11.e0 + 2.e0*vkms)
Mlf2(Q_Qstar) = Mtot*0.008e0*((Q_Qstar)**1.4e0*exp(-(Q_Qstar/3.0e0)**1.6e0))

set xl "Q_Qstar"
set yl "Mlf [Mtot]"

set xr [1e-2:20]
set yr [1e-4:1e-1]

set logscale xy

p "Sevecek2017_Fig4_Mlf.dat" w p,\
  Mlf1(x),\
  Mlf2(x)

pa -1



