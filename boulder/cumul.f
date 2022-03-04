
c---------------------------------------------------------------

      real*8 function cumul(rmax,r1,r2,qa,qb,qc,r_norm,n_norm,rad)
      real*8 rmax,r1,r2,qa,qb,qc,r_norm,n_norm,rad

      if(rad.gt.rmax) then
         cumul = 0.d0
      end if

      if(rad.le.rmax.and.rad.gt.r1) then
         cumul = n_norm * (rad/r_norm)**qa
      end if

      if(rad.le.r1.and.rad.gt.r2) then
         cumul = n_norm * (r1/r_norm)**qa * (rad/r1)**qb
      end if

      if(rad.le.r2) then
         cumul = n_norm*(r1/r_norm)**qa *(r2/r1)**qb *(rad/r2)**qc
      end if

      return

      end 

