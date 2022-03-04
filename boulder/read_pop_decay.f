c read_pop_decay.f

      subroutine read_pop_decay(Nanuli,n,tpart,npart)

      include 'ucrm3.4.inc'

      integer Nanuli,n
      real*8 tpart(npart_max),npart(npart_max,Manuli)

      integer i,j,ierr
      character*80 filename

c (e) dynamical decay data
      filename = 'pop_decay.dat'
      open(unit=1,file=filename,status='old')
      call skip(1)

      i=0
5     continue
        i=i+1
        if (i.gt.npart_max) then
          write(*,*) '# read_pop_decay: index i = ', i,
     :      ' .gt. npart_max = ', npart_max
          stop
        endif

        read(1,*,iostat=ierr,end=15) tpart(i),(npart(i,j),j=1,Nanuli)

        if (ierr.ne.0) then
          write(*,*) '# read_pop_decay: error reading ', filename
          stop
        endif
      goto 5

15    continue
      close(1)
      n=i-1

      if (dbg) then
        write(*,*) '# read_pop_decay: n = ', n
      endif

      return
      end


