c linterp.f
c A linear interpolation of the function y(x) given as two arrays.
c Miroslav Broz (miroslav.broz@email.cz), Mar 2nd 2011

      real*8 function linterp(x,y,n,x0,extra)

      implicit none
      integer n
      real*8 x(n),y(n),x0
      logical extra

      integer i
      real*8 y0

      if (x0.lt.x(1)) then
        extra = .true.
        i=2
      elseif (x0.gt.x(n)) then
        extra = .true.
        i=n
      else
        extra = .false.
        i=2
        do while (x(i).lt.x0)
          i=i+1
        enddo
      endif
      y0 = y(i-1) + (y(i)-y(i-1)) * (x0-x(i-1))/(x(i)-x(i-1))

      linterp = y0
      return
      end

