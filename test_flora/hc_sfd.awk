#!/usr/bin/awk -f

# hc_sfd.awk
# Construct a cumulative histogram from a differential one (gen_sfd output).
# Miroslav Broz (miroslav.broz@email.cz), Feb 24th 2011

(NR>1){
  i++;
  D[i] = 2.*$2;	# km
  dN[i] = $3;	# -
}
END{
  n=i;
  N_gt_D = 0.;
  for (i=n; i>=1; i=i-1) {
    N_gt_D = N_gt_D + dN[i];
    print D[i], N_gt_D;
  }
}

