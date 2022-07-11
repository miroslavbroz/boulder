cccccccccccccccccccccccccccccccccccccccc
      subroutine pop_change(i,j,nbins,mpop,npop,
     %     mpop_change,npop_change,mpop_change_tot,npop_change_tot,
     %     fe,fi,fe_tot,fi_tot,dtcol,dtstern,ifail)
c     bufferizes the population changes in number and mass. 
c     Checks that the population  does not decrease negative and computes 
c     the appropriate timestep for the population not to decrease more than 
c     tollS

      include 'ucrm3.4.inc'
c     inputs
      integer i,j,nbins(Manuli)
      real*8 mpop(Manuli,0:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 npop_change(2,0:BINMAX,Ndata), mpop_change(2,0:BINMAX)
      real*8 dtcol
      real*8 fe(2,0:BINMAX),fi(2,0:BINMAX)
c     inputs/Outputs
      real*8 npop_change_tot(Manuli,0:BINMAX,Ndata)
      real*8 mpop_change_tot(Manuli,0:BINMAX)
      real*8 fe_tot(Manuli,0:BINMAX),fi_tot(Manuli,0:BINMAX)
      real*8 dtstern
      integer ifail
c     internals
      integer jj
      real*8 dtstern_tmp

      ifail=0
      do jj = 1,nbins(j)   ! first check the target bins. Notice, check is done only if population decrease. increasing the population > tollS is fine. For instance if a bin has one object the number and mass can double!
         if(mpop_change(1,jj).lt.0.d0.and.mpop(j,jj).gt.0.d0)then
            dtstern_tmp=dtcol*tollS/abs(mpop_change(1,jj)/mpop(j,jj))
            dtstern=min(dtstern_tmp,dtstern)
            if(dtstern_tmp.lt.dtcol*tollS)then
               dtcol=dtstern
               if (dbg) then
                 write(*,*)'# dtstern = ',dtstern,' y too small'
               endif
               ifail=1
               return
            endif
         endif
         if(npop_change(1,jj,1).lt.0.d0.and.npop(j,jj,1).gt.0.d0)then
            dtstern_tmp=dtcol*tollS/abs(npop_change(1,jj,1)/
     %           npop(j,jj,1))
            dtstern=min(dtstern_tmp,dtstern)
            if(dtstern_tmp.lt.dtcol*tollS)then
               dtcol=dtstern
               if (dbg) then
                 write(*,*)'# dtstern = ',dtstern,' y too small'
               endif
               ifail=1
               return
            endif
         endif
      end do
      do jj = 1,nbins(i)   ! then check the projectile bins
         if(mpop_change(2,jj).lt.0.d0.and.mpop(i,jj).gt.0.d0)then
            dtstern_tmp=dtcol*tollS/abs(mpop_change(2,jj)/mpop(i,jj))
            dtstern=min(dtstern_tmp,dtstern)
            if(dtstern_tmp.lt.dtcol*tollS)then
               dtcol=dtstern
               if (dbg) then
                 write(*,*)'# dtstern = ',dtstern,' y too small'
               endif
               ifail=1
               return
            endif
         endif
         if(npop_change(2,jj,1).lt.0.d0.and.npop(i,jj,1).gt.0.d0)then
            dtstern_tmp=dtcol*tollS/abs(npop_change(2,jj,1)
     %           /npop(i,jj,1))
            dtstern=min(dtstern_tmp,dtstern)
            if(dtstern_tmp.lt.dtcol*tollS)then
               dtcol=dtstern
               if (dbg) then
                 write(*,*)'# dtstern = ',dtstern,' y too small'
               endif
               ifail=1
               return
            endif
         endif
      end do
      
c...  The step was successfull for this pair  (no return hit), update global change variables
      do jj = 0,nbins(j)
         mpop_change_tot(j,jj) = mpop_change_tot(j,jj) 
     &        + mpop_change(1,jj)
         npop_change_tot(j,jj,1) = npop_change_tot(j,jj,1) 
     &        + npop_change(1,jj,1)
         fe_tot(j,jj)=fe_tot(j,jj)+fe(1,jj)
         fi_tot(j,jj)=fi_tot(j,jj)+fi(1,jj)
      end do
      do jj = 0,nbins(i)
         mpop_change_tot(i,jj) = mpop_change_tot(i,jj) 
     &        + mpop_change(2,jj)
         npop_change_tot(i,jj,1) = npop_change_tot(i,jj,1) 
     &        + npop_change(2,jj,1)
         fe_tot(i,jj)=fe_tot(i,jj)+fe(2,jj)
         fi_tot(i,jj)=fi_tot(i,jj)+fi(2,jj)
      end do
      
      return
      end

