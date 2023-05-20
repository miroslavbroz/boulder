#!/usr/bin/gnuplot

load "Q.lab"

# the scaling law is given in CGS units

Q_star_D(r1eff,rho,Q_0,a_benz,BB,b_benz,rho_bas,q_fact) = (Q_0*r1eff**a_benz + BB*rho_bas*r1eff**b_benz)/q_fact

########################################################################

km=1.e3; m =1.e0;  cm=1.e-2; mu=1.e-6
km=1.e0; m =1.e-3; cm=1.e-5; mu=1.e-9

set term x11

set xl "{/Helvetica-Oblique D} [km]"
set yl "{/Helvetica-Oblique Q}^{*}_D [erg/g]"

set xr [100.*mu:1000.*km]
set yr [1e5:1.e11]
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
  "10^{5}"  1e5  0,\
  ""        2e5  1,\
  ""        3e5  1,\
  ""        4e5  1,\
  ""        5e5  1,\
  ""        6e5  1,\
  ""        7e5  1,\
  ""        8e5  1,\
  ""        9e5  1,\
  "10^{6}"  1e6  0,\
  ""        2e6  1,\
  ""        3e6  1,\
  ""        4e6  1,\
  ""        5e6  1,\
  ""        6e6  1,\
  ""        7e6  1,\
  ""        8e6  1,\
  ""        9e6  1,\
  "10^{7}"  1e7  0,\
  ""        2e7  1,\
  ""        3e7  1,\
  ""        4e7  1,\
  ""        5e7  1,\
  ""        6e7  1,\
  ""        7e7  1,\
  ""        8e7  1,\
  ""        9e7  1,\
  "10^{8}"  1e8  0,\
  ""        2e8  1,\
  ""        3e8  1,\
  ""        4e8  1,\
  ""        5e8  1,\
  ""        6e8  1,\
  ""        7e8  1,\
  ""        8e8  1,\
  ""        9e8  1,\
  "10^{9}"  1e9  0,\
  ""        2e9  1,\
  ""        3e9  1,\
  ""        4e9  1,\
  ""        5e9  1,\
  ""        6e9  1,\
  ""        7e9  1,\
  ""        8e9  1,\
  ""        9e9  1,\
  "10^{10}" 1e10 0,\
  ""        2e10 1,\
  ""        3e10 1,\
  ""        4e10 1,\
  ""        5e10 1,\
  ""        6e10 1,\
  ""        7e10 1,\
  ""        8e10 1,\
  ""        9e10 1,\
  "10^{11}" 1e11 0 \
  )
set key left width -10 samplen 1.5 spacing 1.05

set arrow from 1,graph 0 rto 0,graph 1 nohead lt 0
set arrow from 0.001,graph 0 rto 0,graph 1 nohead lt 0

set lmargin 7.0
set rmargin 2.0
set bmargin 3.1
set tmargin 0.8

x1=1.e2*cm*2
x2=2.e7*cm*2
y1=1.e4
y2=2.e9

p Q_star_D(x/2./cm,rho_1,Q_0_1,a_benz_1,BB_1,b_benz_1,rho_bas_1,q_fact_1) lw 7 dt 1 lc '#0066ff' t "MB",\
  Q_star_D(x/2./cm,rho_2,Q_0_2,a_benz_2,BB_2,b_benz_2,rho_bas_2,q_fact_2) lw 5 dt 2 lc '#aaaaaa' t "NEO",\
  Q_star_D(x/2./cm,rho_3,Q_0_3,a_benz_3,BB_3,b_benz_3,rho_bas_3,q_fact_3) lw 3 dt 1 lc '#000099' t "Veritas",\
  Q_star_D(x/2./cm,rho_4,Q_0_4,a_benz_4,BB_4,b_benz_4,rho_bas_4,q_fact_4) lw 2 dt 1 lc '#ffaaff' t "Karin",\
  Q_star_D(x/2./cm, 3.0, 9.0e7, -0.36, 0.5 , 1.36, 3.0, 1.0) t "{/=14 basalt, 5 km/s, 3 g/cm^3}" lw 1 dt 2 lc 'gray',\
  Q_star_D(x/2./cm, 3.0, 7.0e7, -0.45, 2.1 , 1.19, 3.0, 1.0) t "{/=14 ice, 0.5 km/s, 3 g/cm^3}"  lw 1 dt 3 lc 'black'

pa -1

set term post eps enh color solid "Helvetica,18"
set out "Q.eps"
set size 0.675,0.9
rep

system("eps2png Q.eps")

q

