
ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   computes ejected mass, slope of fragments and size of the largest
c   fragment according to SPH experiments
cccccccccccccccccccccc

      subroutine ucrm_ejecta_mass(j,i,m1,m2,vcoll,Mlr,Mej,Mlf,qf,
     :  Q_Qstar)
      include 'ucrm3.4.inc'
c     input
      real*8 m1,m2,vcoll
      integer j,i
c     Output
      real*8 Mlr,Mej,Mlf,qf,Q_Qstar
c     internals
      real*8 Qstar,Q,r1eff,d1,vkms,Mtot
      real*8 Mlr1,Mlr2,Mlf1,Mlf2,qf1,qf2
      integer k,l
      real*8 rho_m,Q0_m,a_benz_m,BB_m,rho_bas_m,b_benz_m,q_fact_m
c functions
      real*8 interp

      if(ifrag.eq.0)then
        mej = 0.d0
        mlf = 0.d0
        qf = -2.5d0
        Q_Qstar = 0.d0
        return
      endif

      Mtot = m1+m2
      rho_m = (m1*rho(j)+m2*rho(i))/Mtot
      Q0_m = (m1*Q0(j)+m2*Q0(i))/Mtot
      a_benz_m = (m1*a_benz(j)+m2*a_benz(i))/Mtot
      BB_m = (m1*BB(j)+m2*BB(i))/Mtot
      rho_bas_m = (m1*rho_bas(j)+m2*rho_bas(i))/Mtot
      b_benz_m = (m1*b_benz(j)+m2*b_benz(i))/Mtot
      q_fact_m = (m1*q_fact(j)+m2*q_fact(i))/Mtot
      r1eff = (Mtot/(4.d0/3.d0*pi*rho_m))**(1.d0/3.d0)
!      vkms = vcoll*1.d-5
      vkms = 5.d0  ! use nominal speed to prevent Mlr+Mlf>Mtot!

c scaling law
      Qstar = Q0_m*r1eff**a_benz_m+BB_m*rho_bas_m*r1eff**b_benz_m
      Qstar = Qstar/q_fact_m
      Q = 0.5d0*m2*vcoll**2/m1
      Q_Qstar = Q/Qstar

c largest remnant
      Mlr1 = Mtot/(1.d0+(0.6d0+56.d0*exp(-1.d0*vkms))
     :    * Q_Qstar**(0.8d0+8.d0*exp(-0.7d0*vkms)))     ! Sevecek etal. (2017)
      if (Q.lt.Qstar) then
        Mlr2 = (-0.5d0*(Q_Qstar-1.d0)+0.5d0)*Mtot ! Benz & Asphaug (1999)
      else
        Mlr2 = (-0.35d0*(Q_Qstar-1.d0)+0.5d0)*Mtot ! Benz & Asphaug (1999)
      endif
      d1 = 2.d0*(m1/(4.d0/3.d0*pi*rho(j)))**(1.d0/3.d0)
      Mlr = interp(d1,10.d5,100.d5,Mlr1,Mlr2)
      if (Mlr.lt.0.d0) Mlr=0.d0
      Mej = (m1+m2)-Mlr
         
c largest fragment
      Mlf1 = Mtot/(0.24d0*vkms**3
     :    * (Q_Qstar)**(-0.6d0-2.d0*exp(-0.3d0*vkms))
     :    + exp(-0.3d0*vkms)*Q_Qstar + 11.d0 + 2.d0*vkms) ! Sevecek etal. (2017)
      Mlf2 = 0.008d0*((Q_Qstar)**1.4d0*exp(-(Q_Qstar/3.0d0)**1.6d0)) ! Vernazza etal. (2018), App. A
     :    * Mtot
      Mlf = interp(d1,10.d5,100.d5,Mlf1,Mlf2)
      if (Mlf.gt.Mej) Mlf=0.d0
      if (Mlf.ge.Mlr) then
        Mlr = 0.d0
        Mej = (m1+m2)
      endif

c slope
      qf1 = -12.3d0+0.75d0*vkms+(11.5d0-1.d0*vkms)*exp(-5.d-3*Q_Qstar)
     :   / (1.d0+0.1*Q_Qstar**(-0.4d0))                               ! Sevecek etal. (2017)
      qf2 = -10.d0+7.d0*(Q_Qstar)**0.4d0*exp(-Q_Qstar/7.d0) ! Durda etal. (2007)
      qf = interp(d1,10.d5,100.d5,qf1,qf2)

c      write(*,*) 'vkms = ', vkms, ' km/s'  ! dbg
c      write(*,*) 'd1 = ', d1/1.d5, ' km'  ! dbg
c      write(*,*) 'Q_Qstar = ', Q_Qstar  ! dbg
c      write(*,*) 'Mlr/Mtot: ', Mlr1/Mtot, Mlr2/Mtot, Mlr/Mtot  ! dbg
c      write(*,*) 'Mlf/Mtot: ', Mlf1/Mtot, Mlf2/Mtot, Mlf/Mtot  ! dbg
c      write(*,*) 'qf      : ', qf1, qf2, qf  ! dbg

      return
      end


      real*8 function interp(x,x1,x2,y1,y2)
      implicit none
      real*8 x,x1,x2,y1,y2
      if (x.le.x1) then
        interp = y1
      else if (x.ge.x2) then
        interp = y2
      else
        interp = y1 + (y2-y1)*(x-x1)/(x2-x1)
      endif
      return
      end


