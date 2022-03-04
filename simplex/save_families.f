c save_families.f
c Save information about large disruptions to an array.
c Miroslav Broz (miroslav.broz@email.cz), Aug 24th 2011

      subroutine save_families(j,time,sarr_kl,sarr_il,ncoll_kl,vcoll,
     :     mtot,mlr,mej,mlf,qf,rho_target)

      include '../boulder/ucrm3.4.inc'
      include 'cb_fam.inc'

      integer j
      real*8 time,sarr_kl,sarr_il,ncoll_kl,vcoll,mtot,mlr,mej,mlf,qf,
     :  rho_target

      real*8 m_2nd,dlf,dlr

c choose either largest remnant or largest fragment, whichever is larger
      if (mlr.gt.mlf) then
        m_2nd = mlr
      else
        m_2nd = mlf
      endif
      dlr = 2.d0 * (mlr/(4.d0/3.d0*pi*rho_target))**(1.d0/3.d0)
      dlf = 2.d0 * (mlf/(4.d0/3.d0*pi*rho_target))**(1.d0/3.d0)

      if ((2.d0*sarr_kl.ge.family_dpb_treshold)
     :  .and.(m_2nd/mtot.gt.family_lfpb_min)
     :  .and.(m_2nd/mtot.lt.family_lfpb_max)
     :  .and.(dlf.ge.family_dlf_treshold)) then

        FAM_n = FAM_n + 1
        if (FAM_n.gt.FAM_MAX) then
          write(*,*) '# Warning: FAM_n > FAM_MAX = ', FAM_MAX
          FAM_n = FAM_MAX       ! do NOT stop here, due to simplex...
        endif
        FAM_pop(FAM_n) = j
        FAM_time(FAM_n) = time/1.d6                     ! yr -> Myr
        FAM_dtarget(FAM_n) = 2.d0*sarr_kl/1.d5          ! cm -> km
        FAM_dproject(FAM_n) = 2.d0*sarr_il/1.d5         ! cm -> km
        FAM_ncoll_kl(FAM_n) = ncoll_kl
        FAM_vcoll(FAM_n) = vcoll/1.d5                   ! cm/s -> km/s
        FAM_dlr(FAM_n) = dlr/1.d5                       ! cm -> km
        FAM_dlf(FAM_n) = dlf/1.d5                       ! cm -> km
        FAM_m2nd_over_mtot(FAM_n) = m_2nd/mtot

      endif

      return
      end


