c sort2.f
c Minimal sort of two arrays.
c Miroslav Broz (miroslav.broz@email.cz), Nov 15th 2009

      subroutine sort2(n,x,y)
      integer n
      real*8 x(n),y(n)

      integer j,k,min
      real*8 tmp

      do j = 1,n
        min=j
        do k = j+1,n
          if (x(min).gt.x(k)) then
            min=k
          endif
        enddo
        tmp=x(j)
        x(j)=x(min)
        x(min)=tmp

        tmp=y(j)
        y(j)=y(min)
        y(min)=tmp
      enddo
      
      return
      end

