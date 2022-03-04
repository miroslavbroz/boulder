c write_craters.f
c Write individual craters for selected bodies.
c Miroslav Broz (miroslav.broz@email.cz), May 30th 2019

c j ... target pop.
c i ... projectile
c k ... bin
c l ... another bin

      subroutine write_craters(j,time,sarr_jk,sarr_il,ncoll_kl,
     :  rho_j,rho_i,vcoll,npop_jk)

      include 'ucrm3.4.inc'
      integer j
      real*8 time,sarr_jk,sarr_il,ncoll_kl,rho_j,rho_i,vcoll,npop_jk

c locals
      integer i1st,icra,ierr
      real*8 D_c,ncoll,phi
      save i1st,icra
      data icra /0/
      data i1st /0/

c functions
      real*8 piscaling

      if ((ncoll_kl.gt.0.d0).and.(2.d0*sarr_jk.gt.crater_dpb_min)) then

        phi = pi/4.d0
        D_c = piscaling(2.d0*sarr_jk,2.d0*sarr_il,rho_j,rho_i,
     :    vcoll,phi)
        ncoll = ncoll_kl/npop_jk

        if (D_c.gt.crater_dc_min) then
          icra = icra+1

          open(unit=15,file="craters.out",access="append",
     :      status="unknown", iostat=ierr)

          if (ierr.ne.0) then
            write(*,*) 'write_craters.f: Error opening craters.out.'
            stop
          endif

          if (i1st.eq.0) then
            write(15,*) "# crater_number & ",
     :        "pop_id & ",
     :        "time [Myr] & ",
     :        "D_target [km] & ",
     :        "D_projectile [km] & ",
     :        "D_crater [km] & ",
     :        "number_of_craters_per_target & ",
     :        "v_coll [km/s] & "
            i1st = 1
          endif

          write(15,20) icra, j, time/1.d6,
     :      2.d0*sarr_jk/1.d5, 2.d0*sarr_il/1.d5, D_c/1.d5,
     :      ncoll, vcoll/1.d5
20        format(i8.8,1x,i4,1x,f12.4,
     :      3(1x,f12.3),1x,e16.8,1x,f12.4)
         
          close(15)

        endif
      endif

      return
      end


