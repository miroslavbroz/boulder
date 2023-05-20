#!/usr/bin/awk -f

BEGIN{
  i=0;
}
!/^\/\*/{
  i++;
  prn(i);
}
END{
  for (j=i; j<=4; j++) {
    prn(j);
  }
}

func prn(i){
  gsub("d","e");
  printf("rho_%d    = %f\n", i, $1+0.0);
  printf("Q_0_%d    = %f\n", i, $2+0.0);
  printf("a_benz_%d = %f\n", i, $3+0.0);
  printf("BB_%d     = %f\n", i, $4+0.0);
  printf("b_benz_%d = %f\n", i, $5+0.0);
  printf("rho_bas_%d= %f\n", i, $6+0.0);
  printf("q_fact_%d = %f\n", i, $7+0.0);
}

