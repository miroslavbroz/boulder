#!/usr/bin/awk -f

# hc.awk
# Construct a cumulative histogram from a differential one (ob.data format).
# Miroslav Broz (miroslav.broz@email.cz), Feb 24th 2011

BEGIN{
  i=0;
}
(NF>0){
  i++;
  D[i] = 2.*$4/1.e5;	# cm -> km
  dN[i] = $5/$3;	# total mass/one-body mass
}
(NF==0){
  PRN();
}
END{
  PRN();
}

func PRN(){
  n=i;
  N_gt_D = 0.;
  for (i=n; i>=1; i=i-1) {
    N_gt_D = N_gt_D + dN[i];
    print D[i], N_gt_D;
  }
  i=0;
  print "";
}

