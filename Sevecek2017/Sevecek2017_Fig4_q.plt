#!/usr/bin/gnuplot

vkms = 5.0
Mtot = 1.0

qf1(Q_Qstar) = -12.3e0+0.75e0*vkms+(11.5e0-1.e0*vkms)*exp(-5.e-3*Q_Qstar) / (1.e0+0.1*Q_Qstar**(-0.4e0))
qf2(Q_Qstar) = -10.e0+7.e0*(Q_Qstar)**0.4e0*exp(-Q_Qstar/7.e0)

set xl "Q_Qstar"
set yl "q []"

set xr [1e-2:20]
set yr [-5.5:-1.5]

set logscale x
set samples 1000

p "Sevecek2017_Fig4_q.dat" w p,\
  qf1(x),\
  qf2(x)

pa -1



