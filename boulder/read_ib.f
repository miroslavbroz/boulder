c read_ib.f

      subroutine read_ib(boulder_name, Nanuli,axe,delta_a,nbins,
     :  marr,mpop,sarr,ecc,inc,npop)

      include 'ucrm3.4.inc'

      character*(*) boulder_name
      integer Nanuli,nbins(Manuli)
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX) 
      real*8 mpop(Manuli,0:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 ecc(Manuli,BINNEG:BINMAX),inc(Manuli,BINNEG:BINMAX)
      real*8 axe(Manuli),delta_a(Manuli)

      integer i,jj,ios

c (a) initial populations
      open(unit=1,file=boulder_name,status='old', iostat = ios)
      if (ios .ne. 0) then
       write (6, *) ' No input SFD file, error ', ios, ' Try again '
       stop
      endif
c      call skip(1)
      read(1,*)Nanuli
      do i=1,Nanuli
       read(1,*)axe(i),delta_a(i),nbins(i),marr(i,0),mpop(i,0)
       do jj=1,nbins(i)
        read(1,*)marr(i,jj),sarr(i,jj),mpop(i,jj)
        read(1,*)ecc(i,jj),inc(i,jj)
        npop(i,jj,1)=mpop(i,jj)/marr(i,jj)
       enddo
       ecc(i,0)=ecc(i,1)
       inc(i,0)=inc(i,1)
      enddo
      close(1)

      return
      end


