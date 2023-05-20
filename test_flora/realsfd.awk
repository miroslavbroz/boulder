#!/usr/bin/awk -f

# realsfd.awk
# Modify the SFD produced by gen_sfd to resemble the observed one.
# Miroslav Broz (miroslav.broz@email.cz), Sep 13th 2011

BEGIN{
  eps = 10.;
  fmt = "%12d %20.14f %25.16f     \n";
}
{
  if      (abs($2-479.) < eps) { printf(fmt, $1, $2, $3+1); }
  else if (abs($2-279.) < eps) { printf(fmt, $1, $2, $3+1); }
  else if (abs($2-243.) < eps) { printf(fmt, $1, $2, $3+1); }  # incl. Vesta
  else if (abs($2-213.) < eps) { printf(fmt, $1, $2, $3+0); }
  else if (abs($2-186.) < eps) { printf(fmt, $1, $2, $3+0); }
  else if (abs($2-162.) < eps) { printf(fmt, $1, $2, $3+0); }
  else if (abs($2-142.) < eps) { printf(fmt, $1, $2, $3+0); }
  else { print; }
}


func abs(x) {
  if (x > 0) {
    return x;
  } else {
    return -x;
  }
}


