#!/bin/sh

#make gen_sfd

#           rho [g/cm^3] d1 [km]  d2    d3    d4    qa    qb    qc    qd    qe    dmin        dmax   d_norm  n_norm  boulder_mfactor
./gen_sfd5  3.0          100.     15.   3.6   0.1   -4.55 -2.6  -3.6  -2.7  -2.7  0.000987654 4000.  100.    250.    1.5  > sfd1.dat	# MB-like, less populous to get fit @ t = 1.4 Gyr, adjusted scaling law
./gen_sfd5  3.04         3.0      0.5   0.2   0.1   -4.2  -4.2  -3.5  -2.7  -2.7  0.000987654 4000.  3.0     3000.   1.5  > sfd2.dat	# initial Flora, 2x observed

./realsfd.awk  < sfd1.dat > sfd1.dat_EDITED
./realflora.awk < sfd2.dat > sfd2.dat_EDITED

./gen_ic < gen_ic.in
./gen_sfd.plt

exit




