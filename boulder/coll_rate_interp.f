c coll_rate_interp.f
c Interpolate collisional probabilities and velocities according to the given time.
c Miroslav Broz (miroslav.broz@email.cz), Mar 3rd 2011

      subroutine coll_rate_interp(Nanuli,Pint_time,Pint_arr,vrel_arr,
     :  n,t)

      include 'ucrm3.4.inc'

      integer Nanuli,n
      real*8 t
      real*8 Pint_time(Pint_max)
      real*8 Pint_arr(Pint_max,Manuli,Manuli)
      real*8 vrel_arr(Pint_max,Manuli,Manuli)

      integer i,j
      logical extra
      real*8 linterp	! functions

      do i=1,Nanuli
        do j=1,Nanuli
          Pint(i,j) = linterp(Pint_time,Pint_arr(1,i,j),n,t,extra)

          if (extra) then
            write(*,*) '# coll_rate_interp: extrapolation is NOT',
     :        ' allowed!'
            stop
          endif

          vrel(i,j) = linterp(Pint_time,vrel_arr(1,i,j),n,t,extra)
        enddo
      enddo

      return
      end


