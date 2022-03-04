c read_yarko_tau.f
c Read Yarkovsky timescale tau vs size D dependence.
c Miroslav Broz (miroslav.broz@email.cz), Jun 20th 2013

      subroutine read_yarko_tau(j,filename,yarko_n,yarko_r,yarko_tau)

      include 'ucrm3.4.inc'
      include 'yarko.inc'

      integer j
      character*80 filename

      integer k, ierr

      open(unit=2,file=filename,status='old')
      call skip(2)

      k = 0
5     continue
        k = k + 1
        if (k.gt.yarko_max) then
          write(*,*) '# read_yarko_tau: index k = ', k,
     :      ' .gt. yarko_max = ', yarko_max
          stop
        endif

        read(2,*,iostat=ierr,end=15) yarko_r(k,j), yarko_tau(k,j) ! units: km, Myr

        yarko_r(k,j) = yarko_r(k,j)*1.d5/2.d0 ! km -> cm, diameters -> radii
        yarko_tau(k,j) = yarko_tau(k,j)*1.d6 ! Myr -> yr

        if (ierr.ne.0) then
          write(*,*) '# read_yarko_tau: error reading ', filename
          stop
        endif
      goto 5

15    continue
      yarko_n(j) = k - 1
      close(2)

      write(*,*) '# read_yarko_tau: yarko_n(', j, ') = ', yarko_n(j)

      return
      end


