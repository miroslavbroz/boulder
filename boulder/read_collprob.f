c read_collprob.f
c Output is in the common block /impacts/.

      subroutine read_collprob(Nanuli)

      include 'ucrm3.4.inc'
      integer Nanuli

      integer i,j

c (d) collision probabilities and relative velocities
      open(unit=1,file='collprob.dat',status='old')
      call skip(1)
      do i=1,Nanuli
       do j=1,Nanuli
        read(1,*)Pint(i,j),vrel(i,j)
        write(*,*) '# Pint(', i, ',', j, ') = ', Pint(i,j),' km^-2 y^-1'
        write(*,*) '# vrel(', i, ',', j, ') = ', vrel(i,j),' km/s'

        Pint(i,j)=Pint(i,j)/(1.d5)**2       ! notice, we convert to cm^2 here (cm/s; cm^-2 y^-1)
        vrel(i,j)=vrel(i,j)*1.d5
       enddo
      enddo
      close(1)

      return
      end

