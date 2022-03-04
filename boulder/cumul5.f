
c---------------------------------------------------------------

      real*8 function cumul5(rmax,r1,r2,r3,r4,qa,qb,qc,qd,qe,
     :  r_norm,n_norm,rad)

      implicit none
      real*8 rmax,r1,r2,r3,r4,qa,qb,qc,qd,qe,r_norm,n_norm,rad

      if(rad.gt.rmax) then
        cumul5 = 0.d0
      end if

      if(rad.le.rmax.and.rad.gt.r1) then
        cumul5 = n_norm * (rad/r_norm)**qa
      end if

      if(rad.le.r1.and.rad.gt.r2) then
        cumul5 = n_norm * (r1/r_norm)**qa * (rad/r1)**qb
      end if

      if(rad.le.r2.and.rad.gt.r3) then
        cumul5 = n_norm * (r1/r_norm)**qa *(r2/r1)**qb *(rad/r2)**qc
      end if

      if (rad.le.r3.and.rad.gt.r4) then
        cumul5 = n_norm * (r1/r_norm)**qa *(r2/r1)**qb *(r3/r2)**qc
     :    *(rad/r3)**qd
      end if

      if (rad.le.r4) then
        cumul5 = n_norm * (r1/r_norm)**qa *(r2/r1)**qb *(r3/r2)**qc
     :    *(r4/r3)**qd * (rad/r4)**qe
      end if

      return

      end 

