
c----------------------------------------------------------------
c module on gravitational focussing. Using Greeinzweigh and lissauer 1990
c------------------------------------------------------------------------

      subroutine focusing_factor(m1,r1,m2,r2,a1,vcoll,fg)
      include 'ucrm3.4.inc'
c     input
      real*8 m1,r1,m2,r2,a1,vcoll
c     output
      real*8 fg  ! gravitational focussing
c     internals
      real*8 vesc,rHill,d,v2k,F2B,sin2i,ri
      real*8 q,rlF2B,rlFG,rlF2B_liss,F2B_liss
      real*8 LogF25,LogF50,d4,d2,ratio,Vhill,beta
      external LogF25,LogF50

c     compute the mutual escape velocity in cm/s
      vesc=sqrt(2.d0*G*(m1+m2)/(r1+r2))
     
c     compute mutual Hill radius (in cm)
      rHill=a1*AU*((m1+m2)/msun/3.d0)**(1.d0/3.d0)

      v2k=(G*msun/a1/AU) ! square of kepler velocity in cm/s

c     stuff to correct for gaussian velocity speed. (see WS93)
      vhill=vcoll/sqrt(v2k)/(rHill/a1/AU)
      if(vhill.ge.2.d0)then
       beta=2.7d0
      elseif(vhill.lt.1.d0)then
       beta=1.0d0
      else
       beta=1.d0+(vhill-1.d0)*1.7d0
      endif

c     compute normalized sum or radii (in Rhill)
      d=(r1+r2)/rHill

c     compute the two body gravitational focussing
      F2B=(1.d0+beta*vesc**2/vcoll**2)

c     compute the effective relative inclination for Lissauer's formulae
c     i_H in greenzweig and Lissauer 1990, which stands for the relative
c     velocity assuming e_H=2i_H

      sin2i=(vcoll**2/v2k)*4.d0/7.d0 ! using formula A1 in Kenyon and Luu 98
      ri=a1*AU*sqrt(sin2i)/rHill  ! i_H in greenzweig and Lissauer 1990

      q=2.2d0*d/ri**2
      rlF2B=log10(F2B)
      F2B_liss=1.d0+1.80178d0/d/ri**2
      rlF2B_liss=log10(F2B_liss)
      if(ri.le.0.25d0)then
       rlFg=LogF25(ri,d)/rlF2B_liss*rlF2B
      elseif(ri.le.0.5d0)then
       d4=d*(4.d0*ri)**(-2)
       d2=d*(2.d0*ri)**(-2)
       rlFg=(2.d0-4.d0*ri)*LogF25(0.25d0,d4)+(4.d0*ri-1.d0)
     :   *LogF50(0.5d0,d2)
       rlFg=rlFg/rlF2B_liss*rlF2B
      elseif(ri.le.2.5d0)then
       rlFg=LogF50(ri,d)/rlF2B_liss*rlF2B
      else
       rlFg=rlF2B
      endif
 
      fg=10.d0**rlFg

      return
      end

c-----------------------------------------------------------------------
c     auxiliary functions for gravitational focussing calculation
      real*8 function LogF25(ri,d)
      implicit none
      real*8 ri,d,q,rlFg
         
      q=2.2d0*d/ri**2
      if(ri.le.0.25d0.and.q.le.0.4d0)then
         rlFg=0.725d0-log10(d)-log10(ri)
      elseif(ri.le.0.25d0.and.q.le.1.d0)then
         rlFg=0.96d0-0.685d0*log10(d)-1.63d0*log10(ri)
      elseif(ri.le.0.25d0.and.q.le.6.3d0)then
         rlFg=1.195d0-(1.375d0+0.219d0*ri)*log10(q)-3.d0*log10(ri)
      elseif(ri.le.0.25d0.and.q.gt.6.3d0)then
         rlFg=0.78d0-0.175d0*ri-1.5d0*log10(d)
      else
         write(*,*)'i too big'
         stop
      endif

      LogF25=rlFg
      return
      end

c-----------------------------------------------------------------------
      real*8 function LogF50(ri,d)
      implicit none
      real*8 ri,d,q,f2B,rlF2B,rlFg
         
      q=2.2d0*d/ri**2
      F2B=1.d0+1.80178d0/d/ri**2
      rlF2B=log10(F2B)
      if(ri.le.2.5d0)then
         rlFg=rlF2B+(0.0971d0+0.0794d0*log10(ri))*(ri-2.5d0)**2
         if(log10(d).gt.ri-3.d0)then
            rlFg=rlFg-0.0216d0*(3.d0-ri+log10(d))*(ri-2.5d0)**2
         endif
      else
         write(*,*)' i is too big'
         stop
      endif
      LogF50=rlFg
      return
      end

