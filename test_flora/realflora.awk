#!/usr/bin/awk -f

# realeos.awk
# Modify the SFD produced by gen_sfd to resemble the observed one.
# Miroslav Broz (miroslav.broz@email.cz), Sep 13th 2011

BEGIN{
  eps =  0.1;
  fmt = "%12d %20.14f %25.16f     \n";
}
{
  if      (abs($2-100.0) < eps) { printf(fmt, $1, $2, $3+0); }
  else if (abs($2-  6.3) < eps) { printf(fmt, $1, $2, $3-0); }
  else if (abs($2-  7.2) < eps) { printf(fmt, $1, $2, $3-0); }
  else if (abs($2-  8.3) < eps) { printf(fmt, $1, $2, $3-3); }
  else if (abs($2-  9.5) < eps) { printf(fmt, $1, $2, $3-1); }
  else if (abs($2- 10.8) < eps) { printf(fmt, $1, $2, $3-2); }
  else if (abs($2- 12.4) < eps) { printf(fmt, $1, $2, $3-1); }
  else if (abs($2- 14.2) < eps) { printf(fmt, $1, $2, $3+1); }
  else if (abs($2- 28.0) < eps) { printf(fmt, $1, $2, $3+1); }
  else if (abs($2- 72.2) < eps) { printf(fmt, $1, $2, $3+1); }
  else { print; }
}


func abs(x) {
  if (x > 0) {
    return x;
  } else {
    return -x;
  }
}


