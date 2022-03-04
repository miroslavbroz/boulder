cccccccccccccccccccccccccccccccccccc
      subroutine update_pops(Nanuli,nbins,mpop,npop,ecc,inc,
     %     mpop_change_tot,npop_change_tot,fe_tot,fi_tot)
c     update masses,number, e and i of all populations
      include 'ucrm3.4.inc'
c     inputs
      integer Nanuli,nbins(Manuli)
      real*8 npop_change_tot(Manuli,0:BINMAX,Ndata)
      real*8 mpop_change_tot(Manuli,0:BINMAX)
      real*8 fe_tot(Manuli,0:BINMAX),fi_tot(Manuli,0:BINMAX)
c     input/outputs
      real*8 ecc(Manuli,BINNEG:BINMAX),inc(Manuli,BINNEG:BINMAX)
      real*8 mpop(Manuli,0:BINMAX),npop(Manuli,BINNEG:BINMAX,Ndata)
c     internals
      integer j,jj
c ...............................................................
      do j = 1,Nanuli  
       do jj = 0,nbins(j) 
        ecc(j,jj) = mpop(j,jj)*ecc(j,jj)+fe_tot(j,jj)
        inc(j,jj) = mpop(j,jj)*inc(j,jj)+fi_tot(j,jj)
        npop(j,jj,1) = npop(j,jj,1) + npop_change_tot(j,jj,1) 
        mpop(j,jj) = mpop(j,jj) + mpop_change_tot(j,jj)
        if(mpop(j,jj).gt.0.d0)then
         ecc(j,jj) = ecc(j,jj)/mpop(j,jj)
         inc(j,jj) = inc(j,jj)/mpop(j,jj)
        else
         ecc(j,jj) = 0.d0
         inc(j,jj) = 0.d0
        endif
       enddo
      enddo
c
      return
      end
