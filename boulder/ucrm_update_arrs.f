
c****************************************************************************
c                                  UCRM_UPDATE_ARRS
c****************************************************************************
c     This procedure updates marrs and sarrs using new marr and mpop
c     (the moving-bin algorithm) 
c****************************************************************************

      subroutine ucrm_update_arrs(i1st,j,nbinneg,nbins,marr,sarr,mpop
     &     ,npop,axe,da,ecc,inc)

      include 'ucrm3.4.inc'

c...  Inputs and Outputs
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX)
      real*8 mpop(Manuli,0:BINMAX),npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 ecc(Manuli,BINNEG:BINMAX),inc(Manuli,BINNEG:BINMAX)
      real*8 axe(Manuli),da(Manuli)
c...  Inputs
      integer j,nbins,nbinneg,i1st

c...  Internals
      integer k
      integer jj,kk,nint,nearest,nearest_lower,nearest_higher,difference
      real*8 onethird,ftp,temp
      integer nbin_max,ibin_plus,nbins_new
      real*8 Msmallest,r_low,r_high
      real*8 pop_ref,msum,nsum,slope,nsumnew,msumnew
      real*8 rl3,rh3,rl1,rh1,rc3,rc1,dndr3,dndr1,sl,slmin,ncol,vcol
      real*8 Mminus,eminus,iminus
      save ftp,onethird

c...  Init saved Variables
      if (i1st.eq.0) then
       ftp      = 4.d0/3.d0*acos(-1.d0)
       onethird = 1.d0/3.d0
      endif
     
c...  For real bins: new mass bin size = mass/pop 
      do jj = 1,nbins  
       if(npop(j,jj,1).lt.1.d4)then ! too big for an integer otherwise
        nint=int(npop(j,jj,1)+0.5d0) 
        npop(j,jj,1)=nint
       endif
       if(npop(j,jj,1).ge.1.d0)then
        marr(j,jj)=mpop(j,jj)/npop(j,jj,1)
       else
        mpop(j,jj)=0.d0
       endif
      enddo

c...  (2) Check for out of order marr bins and re-sort
 993  continue
      do jj = 2,nbins 
       if (marr(j,jj).lt.marr(j,jj-1)) then              
        if (dbg) then
         write(*,*) '# Re-sorting Bins at j,jj = ',j,jj
        endif

c (a) SFD bins
        temp = marr(j,jj)
        marr(j,jj) = marr(j,jj-1)
        marr(j,jj-1) = temp              

        temp = mpop(j,jj)
        mpop(j,jj) = mpop(j,jj-1)
        mpop(j,jj-1) = temp
              
        temp = npop(j,jj,1)
        npop(j,jj,1) = npop(j,jj-1,1)
        npop(j,jj-1,1) = temp
                          
        temp = ecc(j,jj)
        ecc(j,jj) = ecc(j,jj-1)
        ecc(j,jj-1) = temp
              
        temp = inc(j,jj)
        inc(j,jj) = inc(j,jj-1)
        inc(j,jj-1) = temp

        goto 993            ! start over until we are error free  
       endif
      enddo
     
c... chasing for bins confused with trash bin
 666  continue
      do jj=1,nbins
       if(marr(j,jj).le.marr(j,0))then 
        if (dbg) then
          write(*,*)'# trashing:',marr(j,jj),marr(j,0)
        endif
        if(jj.eq.1.and.mpop(j,2).lt.mpop(j,1)/mfactor**2)then 
         if (dbg) then
           write(*,*)'# WARNING: new bin 1 of annulus ',j,
     :     ' has small mass:',mpop(j,2),marr(j,2)
         endif
        endif
        mpop(j,0)=mpop(j,0)+mpop(j,jj)
        do k=jj+1,nbins
         mpop(j,k-1)=mpop(j,k)
         npop(j,k-1,1)=npop(j,k,1)
         marr(j,k-1)=marr(j,k)
         ecc(j,k-1)=ecc(j,k)
         inc(j,k-1)=inc(j,k)
        enddo
        nbins=nbins-1
        goto 666
       endif
      enddo
      
c... checking if we need to modify nbin by adding new empty bins up to a mass 
c... twice of that of the last filled bin, or removed empty bins 
c... if they are too many
      do jj=nbins,1,-1
       if(mpop(j,jj).gt.0.d0)then
        nbin_max=jj
        goto 881
       endif
      enddo
 881  continue
      ibin_plus=log(2.d0)/log(mfactor)+1.d0
      nbins_new=nbin_max+ibin_plus
      if(nbins_new.gt.nbins)then
       if(nbins_new.gt.BINMAX)then
        write(*,*) '# not enough space for positive bins ',
     :    nbins_new, BINMAX, nbin_max, ibin_plus
        stop
       else
        do jj=nbins+1,nbins_new
         mpop(j,jj)=0.d0
         npop(j,jj,1)=0.d0
         marr(j,jj)=marr(j,jj-1)*mfactor
         ecc(j,jj)=0.d0
         inc(j,jj)=0.d0
        enddo
       endif
      endif
      nbins=nbins_new

c...  (3) Detect bins too far apart and create new bins between
 995  continue
      do jj = 1,nbins-1 
         
c...  If mass ratio is larger that the critical mass ratio, create new bins
       if(marr(j,jj)/marr(j,jj-1).gt.mfactor**1.5d0) then 
c...  First pull all bins forward by one to create a space for new bin 
c...  There has to be space at nbins! 
        nbins=nbins+1
        if (dbg) then
          write(*,*) '# generating bin ',jj,marr(j,jj),marr(j,jj-1)
        endif
        if(nbins.gt.BINMAX)then
         write(*,*) '# not enough space for positive bins'
         stop
        else
         do kk = nbins,jj+1,-1 
          marr(j,kk) = marr(j,kk-1)
          mpop(j,kk) = mpop(j,kk-1)
          npop(j,kk,1) = npop(j,kk-1,1)
          ecc(j,kk)  = ecc(j,kk-1)
          inc(j,kk)  = inc(j,kk-1)
         enddo
            
c...  Derive the mass of objects in the newly-created bin
c   we don't want that bins close to marr(0) are created empty.
c for instance if bin(2) is created empty and then bin(1) is trashed
c the new bin(1) is empty. Thsu the negative tail (see below) will drop empty 
c too. On the contrary, we don't want to fill by interpolation the bins at 
c  Otherwise when a runaway body grows in mass we artificially introduce bodies
c between itself and the rest of the distribution. When bins are fulled or not
c is governed by the parameter Nmul
         marr(j,jj) = sqrt(marr(j,jj-1)*marr(j,jj+1)) ! geometric mean
         if(marr(j,jj).le.Nmul*marr(j,0).and.
     %      npop(j,jj+1,1).gt.Nmul)then ! fill the array by interp
          if (dbg) then
            write(*,*) '# interp ',marr(j,0),marr(j,jj)
          endif
          if(jj.gt.1)then
           Mminus=mpop(j,jj-1)
           eminus=ecc(j,jj-1)
           iminus=inc(j,jj-1)
          else
           Mminus=0.d0
           eminus=ecc(j,jj+1)
           iminus=inc(j,jj+1)
          endif
          msum=mpop(j,jj+1)+Mminus
c     eccentricities and inclinations are filled with the mass weighted mean
          ecc(j,jj)=(mpop(j,jj+1)*ecc(j,jj+1)+Mminus*eminus)/msum
          inc(j,jj)= (mpop(j,jj+1)*inc(j,jj+1)+Mminus*iminus)/msum
c     the mass ig filled by linear interpolation
          slope=(mpop(j,jj+1)-Mminus)/(marr(j,jj+1)-marr(j,jj-1))
          mpop(j,jj)=Mminus+slope*(marr(j,jj)-marr(j,jj-1))
          msumnew=mpop(j,jj+1)+mpop(j,jj)+Mminus !renormalize for mass conservation
          mpop(j,jj+1)=mpop(j,jj+1)/msumnew*msum
          mpop(j,jj)=mpop(j,jj)/msumnew*msum
          Mminus=Mminus/msumnew*msum
          if(jj.gt.1)then
           mpop(j,jj-1)=Mminus
c now compute npop, round npop to an integer, and re-adjust marr
           npop(j,jj-1,1)=mpop(j,jj-1)/marr(j,jj-1)
           if(npop(j,jj-1,1).lt.1.d4)then 
            nint=int(npop(j,jj-1,1)+0.5d0)
            npop(j,jj-1,1)=nint
           endif
           if(npop(j,jj-1,1).ge.1.d0)then
            marr(j,jj-1)=mpop(j,jj-1)/npop(j,jj-1,1)
           else
            mpop(j,jj-1)=0.d0
           endif
          endif
          npop(j,jj,1)=mpop(j,jj)/marr(j,jj)
          if(npop(j,jj,1).lt.1.d4)then 
           nint=int(npop(j,jj,1)+0.5d0) 
           npop(j,jj,1)=nint
          endif
          if(npop(j,jj,1).ge.1)then
           marr(j,jj)=mpop(j,jj)/npop(j,jj,1)
          else
           mpop(j,jj)=0.d0
          endif
          npop(j,jj+1,1)=mpop(j,jj+1)/marr(j,jj+1)
          if(npop(j,jj+1,1).lt.1.d4)then 
           nint=int(npop(j,jj+1,1)+0.5d0) 
           npop(j,jj+1,1)=nint
          endif
          if(npop(j,jj+1,1).ge.1)then
           marr(j,jj+1)=mpop(j,jj+1)/npop(j,jj+1,1)
          else
           mpop(j,jj+1)=0.d0
          endif
         else  ! fill with zeros 
          if (dbg) then
            write(*,*) '# zeroing ',marr(j,0),marr(j,jj)
          endif
          ecc(j,jj)=0.d0
          inc(j,jj)=0.d0
          mpop(j,jj)=0.d0
          npop(j,jj,1)=0.d0
         endif
        endif
        goto 995            ! restart checks
       endif
      enddo
          
c...  (4) Check for the need to join bins 
 996   continue
       do jj = 2,nbins 
c...  Check if bins too close. The test on npops is because one does not want
c...  to create an artificial runaway body in largest bins
        if(npop(j,jj,1).ge.1.d0 .and. npop(j,jj-1,1).ge.1.d0 .and.
     &     marr(j,jj)/marr(j,jj-1) .lt. mfactor**0.5d0) then
         if (dbg) then
           write(*,*) '# Bins too close, joining them, jj = ',jj
         endif
         ecc(j,jj-1) = (ecc(j,jj-1)*mpop(j,jj-1)+
     $        ecc(j,jj)*mpop(j,jj))/(mpop(j,jj-1)+mpop(j,jj))
         inc(j,jj-1) = (inc(j,jj-1)*mpop(j,jj-1)+
     $        inc(j,jj)*mpop(j,jj))/(mpop(j,jj-1)+mpop(j,jj))
         mpop(j,jj-1) = mpop(j,jj-1) + mpop(j,jj)            
         npop(j,jj-1,1) = npop(j,jj-1,1)+npop(j,jj,1)
         marr(j,jj-1) = mpop(j,jj-1)/npop(j,jj-1,1)
c...  Pull the bins beginning at jj+1 each back 1
         do kk  = jj+1,nbins 
          marr(j,kk-1) = marr(j,kk)
          mpop(j,kk-1) = mpop(j,kk)
          npop(j,kk-1,1) = npop(j,kk,1)
          ecc(j,kk-1)=ecc(j,kk)
          inc(j,kk-1)=inc(j,kk)
         enddo
         nbins=nbins-1
         goto 996          ! start over until we are error free 
        endif
       enddo

c... now create the negative tail bins, imposing a slope that is either the extrapolation of the slope in the vicinity of marr(1) or a Donhany slope if the extrapolation cannot be done. This is important to avoid artificial waves in the SFD that would occur in case of a truncated distribution.
      if(mpop(j,1).gt.0.d0)then
       call coll_rate(j,j,
     &      ncol,vcol)   ! compute the collision probability to know how small can be a projectile to kill a target in bin 1.
       Msmallest=.05d0*marr(j,1)*1.d5/(vcol)**2 !compute the mass of the smallest impactor capable to break an object in bin 1, with *D=1.d5 at the current typical collision velocity. Take a factor of 40 for safety (hence .05 instead of 2 in the formula
       nbinneg=log(marr(j,0)/Msmallest)/log(2.d0) ! compute nbinneg using a factor 2. in mass for negative binning.
       if(nbinneg.le.0)nbinneg=-1
      else
       nbinneg=-1
      endif
      if(-nbinneg.lt.BINNEG)then
       write(*,*) '# not enough space for negative bins'
       stop
      elseif(nbinneg.gt.0)then
c here we try an extrapolation of the SFD between bin 1 and bin 23. The value
c 23 is arbitrary, but for massfactor=1.4 it corresponds to approximately 3 
c orders of magnitude in mass and 1 in size. If the extrapolation fails 
c (for instance bin 23 is empty) then a Donhany slope is assumed
       if(npop(j,23,1).gt.0.d0.and.npop(j,1,1).gt.0.d0)then
        rl3=(sqrt(marr(j,23)*marr(j,22))/ftp/rho(j))**onethird
        rh3=(sqrt(marr(j,23)*marr(j,24))/ftp/rho(j))**onethird
        rc3=(marr(j,23)/ftp/rho(j))**onethird
        rl1=(sqrt(marr(j,1)*marr(j,0))/ftp/rho(j))**onethird
        rh1=(sqrt(marr(j,1)*marr(j,2))/ftp/rho(j))**onethird
        rc1=(marr(j,1)/ftp/rho(j))**onethird
        dndr3=npop(j,23,1)/(rh3-rl3)  ! differential number
        dndr1=npop(j,1,1)/(rh1-rl1)
        sl=(log10(dndr3)-log10(dndr1))/(log10(rc3)-log10(rc1))
       else
        sl=-3.5d0
       endif
       sl=sl+1.d0   ! pass from differential to cumulative
       if (dbg) then
         write(*,*) '# tail slope',sl
       endif
c now build the tail distribution. Now that we have the slope we need
c to calibrate the absoulte number of objects (i.e. the vertical scale)
c This is done so that the negative tail extrapolated to bin 1 
c would have the smae number of objects as in npop(1)
       r_high=(sqrt(marr(j,1)*marr(j,2))/ftp/rho(j))**onethird
       r_low=(sqrt(marr(j,1)*marr(j,0))/ftp/rho(j))**onethird
c pop_ref is the normalization coefficient
       pop_ref=-sl*npop(j,1,1)/(r_low**(sl)-r_high**(sl))
       do jj=0,-nbinneg,-1 ! now fill the bins
        if(jj.lt.0)then
         marr(j,jj)=marr(j,jj+1)/2.d0 ! Notice: we force a factor 2 spacing here instead of mfactor
        endif
        r_low=(marr(j,jj)/sqrt(2.d0)/ftp/rho(j))**onethird
        r_high=(sqrt(marr(j,jj)*marr(j,jj+1))/ftp/rho(j))**onethird
        npop(j,jj,1)=pop_ref/sl*(r_high**(sl)-r_low**(sl))
        if(npop(j,jj,1).lt.1.d4)then ! too big for an integer otherwise
           npop(j,jj,1)=int(npop(j,jj,1)+0.5d0)
        endif
        ecc(j,jj)=ecc(j,1)  ! put e and i of bin 1 (which is non-empty otherwise we would not be here
        inc(j,jj)=inc(j,1)
       enddo
      endif

c...  Compute sarr's directly from marr's
      do jj = -nbinneg,nbins
       sarr(j,jj)  = (marr(j,jj)/(ftp*rho(j)))**onethird
      end do

      return
      end                       ! ucrm_checks

