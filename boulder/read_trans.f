c read_trans.f
c Read transport matrix.
c Miroslav Broz (miroslav.broz@email.cz), May 28th 2019

      subroutine read_trans(Nanuli,trans_m)

      include 'ucrm3.4.inc'
      include 'trans.inc'

      integer Nanuli
      integer i,j,ierr

      open(unit=1,file='transport.dat',status='old',iostat=ierr)
      if (ierr.ne.0) then
        write(*,*) '# read_trans: error opening transport.dat'
        stop
      endif
      call skip(1)

      do i = 1,Nanuli
        do j = 1, Nanuli
          read(1,*,iostat=ierr) trans_m(i,j)
          if (ierr.ne.0) then
            write(*,*) '# read_trans: error reading transport.dat'
            stop
          endif
        enddo
      enddo

      return
      end


