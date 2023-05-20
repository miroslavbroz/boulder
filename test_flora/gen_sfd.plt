#!/usr/bin/gnuplot

set colors classic
set term x11

set xl "D [km]"
set yl "N dN"

set xr [0.8*1.e-3:]
set yr [0.8:]
set xtics 10
set mxtics 10
set ytics 10
set mytics 10
set logscale xy

set style line 1 lt 1
set style line 2 lt 3
set style line 3 lt 5
set style line 4 lt 7

set arrow from 1.0,graph 0 rto 0,graph 1 lt 0
set arrow from 0.001,graph 0 rto 0,graph 1 lt 0
set arrow from 100.e-9,graph 0 rto 0,graph 1 lt 0

p \
  "<./hc_sfd.awk sfd1.dat_EDITED" w lp lc '#0066ff' t "initial MB",\
  "<./hc_sfd.awk sfd2.dat_EDITED" w lp lc '#ff6600' t "initial Flora",\
  "<./hc_sfd.awk sfd1.dat" w l lc '#0066ff' not,\
  "<./hc_sfd.awk sfd2.dat" w l lc '#ff6600' not,\
  "<./hc_sfd.awk sfd3.dat" w l lc '#0000cc' not,\
  "<./hc_sfd.awk sfd4.dat" w l lc '#ffaaff' not,\
  "../../dat/Bottke_etal_2015.dat" u (10.**$1):(10.**$2) w p pt 7 lc 'gray' t "Bottke et al. (2015)",\
  "MB_2018.dat_hc"        w l lc 'gray' t "observed MB",\
  "8_Flora_family.dat_hc" w l lc 'gray' t "observed Flora",\
  "<echo 100.e-9 1.2e25" w p pt 1 ps 3 lw 3 lc '#0066ff' t "zodiacal cloud",\
  "<echo 100.e-9 2.7e24" w p pt 1 ps 3 lw 3 lc '#0000cc' t "IRAS 8{/Symbol\260} band",\
  "<echo 100.e-9 1.3e24" w p pt 1 ps 3 lw 3 lc '#cc00ff' t "IRAS 2.1{/Symbol\260} band"

pa -1

set term png small
set out "gen_sfd.png"
rep


q
  "<awk '(FNR==1){ print \"\"; } ($2==1)' ../../MB/[0-9]*/tmp/4000.out | ./hc.awk" u 1:2 t "evolved MB" w l,\

  "../Cibulkova_etal_2014_MB.dat_hc"                 w l,\
  "../Gladman_etal_2009.dat" u (10.**$1):(10.**$2)   w l dt 2

