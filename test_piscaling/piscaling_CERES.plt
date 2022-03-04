#!/usr/bin/gnuplot

deg = pi/180.
km = 1.e3
G = 6.67e-11

# Zajacek (2012):

# pi_D = C_D*pi_2**(-beta)
# pi_D = D*(rho_t/m)**(1./3.)
# pi_2 = 1.61*g*L/V_imp**2

# Schmidt & Housen (1987); competent rock or saturated soil
C_D = 1.6
beta = 0.22

# target; (1) Ceres
rho_t = 2161.
D_t = 946.*km
M_t = 4./3.*pi*(D_t/2.)**3 * rho_t
g = G*M_t/(D_t/2.)**2

# projectile
V_imp = 4.27*km
rho_p = 3000.
m_p(d_p) = 4./3.*pi*(d_p/2.)**3*rho_p

# impact angle
phi = 45.*deg

# crater diameter
D_c(d_p) = C_D*(1.61*g*d_p/V_imp**2)**(-beta) * (m_p(d_p)/rho_t)**(1./3.) * (sin(phi))**(1./3.)

print "g = ", g, " m/s^2"

set table "piscaling.dat"
set xr [1e1:1e5]
p D_c(x) t "Ceres"
unset table

########################################################################

set term x11

set xl "d_p [m]"
set yl "D_c [m]"

set xr [1e1:1e5]
set yr [1e2:1e7]
set logscale x
set logscale y

tmp=1.71*km; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0
tmp=2.87*km; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0
tmp=4.15*km; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0
tmp=13.0*km; set arrow from tmp,graph 0 to tmp,graph 1 nohead lt 0

p D_c(x) t "Ceres",\
  "../MB_Pallas/piscaling.dat" w l dt 2 t "Pallas",\
  13.0*x w l t "Marchi etal. (2016) upper limit",\
   7.6*x w l t "lower",\
  20*km w l lt 0,\
  30*km w l lt 0,\
  40*km w l lt 0,\
  100*km w l lt 0

pa -1

set term png small
set out "piscaling_CERES.png"
rep

q
  11.5*x dt 2,\

