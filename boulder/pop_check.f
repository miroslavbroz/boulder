cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine pop_check(Nanuli,nbins,mpop,npop,
     $     mpop_change_tot,npop_change_tot,dtcol,dtstern,ifail)
c...  Check that no population or mass bin is negative. 
c     Return appropriate timestep for this not to be true next time
      include 'ucrm3.4.inc'

c     Inputs
      integer Nanuli,nbins(Manuli)
      real*8 mpop(Manuli,0:BINMAX),npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 npop_change_tot(Manuli,0:BINMAX,Ndata)
      real*8 mpop_change_tot(Manuli,0:BINMAX)
      real*8 dtcol
c     input/output 
      real*8 dtstern
      integer ifail
c     internals
      integer j,jj
      real*8 dtstern_tmp

      ifail=0
      do j = 1,Nanuli    ! for all anuli 
         do jj = 1,nbins(j)  ! for all bins
            if(mpop_change_tot(j,jj).lt.0.d0.and.mpop(j,jj).gt.0.d0)then
               dtstern_tmp=dtcol*tollS/abs(mpop_change_tot(j,jj)
     %           /mpop(j,jj))
               dtstern=min(dtstern_tmp,dtstern)
               if(dtstern_tmp.lt.dtcol*tollS)then
                  dtcol=dtstern
                  if (dbg) then
                    write(*,*)'# dtstern too small:',dtstern
                  endif
                  ifail=1
                  return
               endif            
            endif
            if(npop_change_tot(j,jj,1).lt.0.d0.and.
     %        npop(j,jj,1).gt.0.d0) then
               dtstern_tmp=dtcol*tollS/abs(npop_change_tot(j,jj,1)/
     %              npop(j,jj,1))
               dtstern=min(dtstern_tmp,dtstern)
               if(dtstern_tmp.lt.dtcol*tollS)then
                  dtcol=dtstern
                  if (dbg) then
                    write(*,*)'# dtstern too small:',dtstern
                  endif
                  ifail=1
                  return
               endif            
            endif
         end do
      end do            
      
      return
      end

