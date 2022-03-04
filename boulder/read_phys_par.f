c read_phys_par.f
c Output is in the common block /phys_par/.

      subroutine read_phys_par(Nanuli)

      include 'ucrm3.4.inc'

      integer Nanuli

      integer j

c (c) physical parameters for disruptions
      open(unit=1,file='phys_par.dat',status='old')
      call skip (1)
      do j=1,Nanuli
       read(1,*) rho(j),Q0(j),a_benz(j),BB(j),b_benz(j),
     :   rho_bas(j),q_fact(j)
      enddo
      close(1)

      return
      end

