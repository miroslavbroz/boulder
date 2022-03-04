#!/usr/bin/gnuplot

vkms = 5.0
vkms = 2.0
Mtot = 1.0

Mlr1(Q_Qstar) = Mtot/(1.e0+(0.6e0+56.e0*exp(-1.e0*vkms)) * Q_Qstar**(0.8e0+8.e0*exp(-0.7e0*vkms)))
Mlr2(Q_Qstar) = Q_Qstar < 1.0 ? (-0.5e0*(Q_Qstar-1.e0)+0.5e0)*Mtot : (-0.35e0*(Q_Qstar-1.e0)+0.5e0)*Mtot

Mlf1(Q_Qstar) = Mtot/(0.24e0*vkms**3 * (Q_Qstar)**(-0.6e0-2.e0*exp(-0.3e0*vkms)) + exp(-0.3e0*vkms)*Q_Qstar + 11.e0 + 2.e0*vkms)
Mlf2(Q_Qstar) = Mtot*0.008e0*((Q_Qstar)**1.4e0*exp(-(Q_Qstar/3.0e0)**1.6e0))

set xl "Q_Qstar"
set yl "Mlr+Mlf [Mtot]"

set xr [1e-8:20]
#set yr [0.2:1.0]
#set xr [1e-3:1e3]
tmp=2.e-3
set yr [-tmp:tmp]

set logscale x
set zeroaxis
set samples 1000

p \
  Mlr1(x)+Mlf1(x)-Mtot,\
  Mlr2(x)+Mlf2(x)-Mtot

pa -1



