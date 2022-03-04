c read_ascii.f
c Read ASCII file in two arrays.
c Miroslav Broz (miroslav.broz@email.cz), Nov 14th 2009

      subroutine read_ascii(filename,n,x,y)

      implicit none

      character*80 filename
      integer n
      real*8 x(n),y(n)
      integer i,length
      character*80 str

      i = 0
      open(unit=10,file=filename,status="old",form="formatted")
5     continue
        read(10,10,err=20,end=20) str
10      format(a)
        if ((str(1:1).ne."#").and.(length(str).gt.0)) then
          i = i+1
          if (i.gt.n) then
            write(*,*) "# Error in read_ascii(): i.gt.n"
          endif
          read(str,*,err=20,end=20) x(i),y(i)
         endif
      goto 5
20    continue
      close(10)

      n = i
      return
      end

