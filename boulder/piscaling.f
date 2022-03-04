c piscaling.f
c Pi-scaling of projectile -> crater diameter (Melosh 1989).
c Miroslav Broz 

      real*8 function piscaling(D_t,d_p,rho_t,rho_p,v_imp,phi)
      include 'ucrm3.4.inc'

      real*8 D_t,d_p,rho_t,rho_p,v_imp,phi

      real*8 C_D,beta
      real*8 M_t,m_p,localg,D_c

      C_D = 1.6d0
      beta = 0.22d0

      M_t = 4.d0/3.d0*pi*(D_t/2.d0)**3 * rho_t
      m_p = 4.d0/3.d0*pi*(d_p/2.d0)**3 * rho_p
      localg = G*M_t/(D_t/2.d0)**2

      D_c = C_D*(1.61d0*localg*d_p/v_imp**2)**(-beta)
     :  * (sin(phi)*m_p/rho_t)**(1.d0/3.d0)

      piscaling = D_c
      return
      end


