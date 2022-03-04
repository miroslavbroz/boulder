c read_yarko.f
c Read yarko.dat file with Nanuli filenames (containing tau(D) data).
c Miroslav Broz (miroslav.broz@email.cz), Jun 20th 2013

      subroutine read_yarko(Nanuli,yarko_n,yarko_r,yarko_tau)

      include 'ucrm3.4.inc'
      include 'yarko.inc'

      integer Nanuli

      integer j,ierr
      character*80 filename

      open(unit=1,file='yarko.dat',status='old')
      call skip (1)

      do j=1,Nanuli

        read(1,*,iostat=ierr) filename
        if (ierr.ne.0) then
          write(*,*) '# read_yarko: error reading yarko.dat'
          stop
        endif

        call read_yarko_tau(j,filename,yarko_n,yarko_r,yarko_tau)

      enddo

      close(1)

      return
      end


