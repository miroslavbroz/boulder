c pop_decay_interp.f

      subroutine pop_decay_interp(Nanuli,tpart,npart,n,t,npart_init)

      include 'ucrm3.4.inc'

      integer Nanuli,n
      real*8 t
      real*8 tpart(npart_max),npart(npart_max,Manuli)
      real*8 npart_init(Manuli,BINNEG:BINMAX)

      integer j,jj
      real*8 tmp
      logical extra
      real*8 linterp	! functions

      do j=1,Nanuli
        tmp = linterp(tpart,npart(1,j),n,t,extra)
        if (extra) then
          write(*,*) '# pop_decay_interp: extrapolation is NOT allowed!'
          stop
        endif

        do jj=BINNEG,BINMAX
          npart_init(j,jj)=tmp
        enddo
      enddo

      return
      end


