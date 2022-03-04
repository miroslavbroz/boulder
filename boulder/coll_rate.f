ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   collision probability calculation

      subroutine coll_rate(i,j,pcoll,vcoll)

c  this subroutine uses pre-computed collision probabilities and impact
c  velocities; projectile and target in anuli *i* and *j*

      include 'ucrm3.4.inc'
c     Inputs
      integer i,j
c     Outputs
      real*8 vcoll,pcoll
c .....................................................................
      pcoll=Pint(i,j)
      vcoll=vrel(i,j)

c the unit conversion was moved to read_collprob.f
c      vcoll=vcoll*1.d5       ! notice, we convert to cm^2 here (cm/s; cm^-2 y^-1)
c      pcoll=pcoll/(1.d5)**2
      return
      end

