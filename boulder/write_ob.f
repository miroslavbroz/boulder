c write_ob.f

      subroutine write_ob(outname,Nanuli,t,nbins,marr,sarr,mpop)

      include 'ucrm3.4.inc'

      character*(*) outname
      integer Nanuli,nbins(Manuli)
      real*8 t
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX) 
      real*8 mpop(Manuli,0:BINMAX)

      integer i,jj

c (a) output SFDs
      open(unit=2,file=outname,position='append')
ccc      write(2,*) t0,total_mass_new
ccc      write(2,*) Nanuli
      do i=1,Nanuli
ccc       write(2,*) axe(i),delta_a(i),nbins(i),marr(i,0),mpop(i,0)
       write(2,*) t,i,nbins(i)
       do jj=1,nbins(i)
        write(2,*) marr(i,jj),sarr(i,jj),mpop(i,jj)
ccc        write(2,*) ecc(i,jj),inc(i,jj)
       enddo
      enddo
      close(2)

      return
      end

