c read_param.f

      subroutine read_param(dtcol,dtmax,t0,tend,dtout)

      real*8  dtcol,dtmax,t0,tend,dtout

c (b) init time, final time, timestep etc.
      open(unit=1,file='param.dat',status='old')
      call skip(1)
      read(1,*) dtcol,dtmax,t0,tend,dtout
      close(1)

      write(*,*) '# dtcol = ', dtcol, ' y = ', dtcol/1.d6, ' My'
      write(*,*) '# dtmax = ', dtmax, ' y = ', dtmax/1.d6, ' My'
      write(*,*) '# t0    = ', t0   , ' y = ', t0/1.d6   , ' My'
      write(*,*) '# tend  = ', tend , ' y = ', tend/1.d6 , ' My'
      write(*,*) '# dtout = ', dtout, ' y = ', dtout/1.d6, ' My'

      return
      end

