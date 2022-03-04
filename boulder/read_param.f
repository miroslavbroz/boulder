c read_param.f

      subroutine read_param(dtcol,dtmax,t0,tend,dtout)

      real*8  dtcol,dtmax,t0,tend,dtout

c (b) init time, final time, timestep etc.
      open(unit=1,file='param.dat',status='old')
      call skip(1)
      read(1,*) dtcol,dtmax,t0,tend,dtout
      close(1)

      write(*,*) 'dtcol = ', dtcol, ' yr'
      write(*,*) 'dtmax = ', dtmax, ' yr'
      write(*,*) 't0 = ', t0, ' yr'
      write(*,*) 'tend = ', tend, ' yr'
      write(*,*) 'dtout = ', dtout, ' yr = ', dtout/1.d6, ' Myr'

      return
      end

