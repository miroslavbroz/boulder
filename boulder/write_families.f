c write_families.f
c Write informations about large disruptions to a file.
c Miroslav Broz (miroslav.broz@email.cz), Aug 24th 2011

      subroutine write_families(j,time,sarr_jk,sarr_il,ncoll_kl,vcoll,
     :     mtot,mlr,mej,mlf,qf,Q_Qstar,rho_target,nbins,deb_spect_m,
     :     deb_spect_n,marr)

      include 'ucrm3.4.inc'
      integer j,nbins
      real*8 time,sarr_jk,sarr_il,ncoll_kl,vcoll,mtot,mlr,mej,mlf,qf,
     :  rho_target,Q_Qstar
      real*8 deb_spect_m(0:BINMAX),deb_spect_n(0:BINMAX) ! debris mass and num.
      real*8 marr(Manuli,BINNEG:BINMAX)

      real*8 m_2nd,dlf,dlr,dlrlf,marr_,sarr_,mpop_,nbins_
      integer i1st,ifam,ierr,jj
      character*80 str
      save i1st, ifam
      data i1st /0/
      data ifam /0/

c choose either largest remnant or largest fragment, whichever is larger
      if (mlr.gt.mlf) then
        m_2nd = mlr
      else
        m_2nd = mlf
      endif
      dlr = 2.d0 * (mlr/(4.d0/3.d0*pi*rho_target))**(1.d0/3.d0)
      dlf = 2.d0 * (mlf/(4.d0/3.d0*pi*rho_target))**(1.d0/3.d0)

c this is a 'fake' parent-body size (essentially a lowest limit)
      dlrlf = 2.d0 * ((mlr+mlf)/(4.d0/3.d0*pi*rho_target))**(1.d0/3.d0)

      if ((2.d0*sarr_jk.ge.family_dpb_treshold)
     :  .and.(m_2nd/mtot.gt.family_lfpb_min)
     :  .and.(m_2nd/mtot.lt.family_lfpb_max)
     :  .and.(dlf.ge.family_dlf_treshold)) then

        ifam = ifam+1

        open(unit=15,file="families.out",access="append",
     :    status="unknown", iostat=ierr)

        if (ierr.ne.0) then
          write(*,*) 'write_families.f: Error opening families.out.'
          stop
        endif

        if (i1st.eq.0) then
          write(15,*) "# family_number & ",
     :      "pop_id & ",
     :      "time [Myr] & ",
     :      "D_target [km] & ",
     :      "D_projectile [km] & ",
     :      "number_of_collisions & ",
     :      "v_coll [km/s] & ",
     :      "D_largest_remnant [km] & ",
     :      "D_largest_fragment [km] & ",
     :      "m_largest_remnant_or_fragment/m_total & ",
     :      "m_ejected/m_total & ",
     :      "qf (cumulative slope of the SFD) & ",
     :      "D_largest_remnant_plus_fragment [km] & ",
     :      "Q/Q_star"
          i1st = 1
        endif

        write(15,20) ifam, j, time/1.d6,
     :    2.d0*sarr_jk/1.d5, 2.d0*sarr_il/1.d5,
     :    ncoll_kl, vcoll/1.d5, dlr/1.d5, dlf/1.d5,
     :    m_2nd/mtot, mej/mtot, qf, dlrlf/1.d5, Q_Qstar
20      format(i8.8,1x,i4,1x,f12.4,2(1x,f12.3),1x,e16.8,1x,f12.4,
     :    8x,2(1x,f12.4),1x,f12.6,1x,f12.6,1x,f12.4,8x,f12.4,
     :    1x,f16.8)

        close(15)

c output also the SFD of the family to a separate file

        if (.FALSE.) then  ! dbg

        open(unit=15,file="families_ob.out",access="append",
     :    status="unknown", iostat=ierr)

        if (ierr.ne.0) then
          write(*,*) 'write_families.f: Error opening families_ob.out'
          stop
        endif

        write(str,"(i8.8)") ifam
        write(15,*) time, 1, nbins, "				ifam = ", str(1:8),
     :    "	ncoll_kl = ", ncoll_kl
        do jj = 1, nbins
          if (deb_spect_n(jj).gt.0.d0) then
            marr_ = deb_spect_m(jj)/deb_spect_n(jj)
            sarr_ = (marr_/(4.d0/3.d0*pi*rho_target))**(1.d0/3.d0)
            mpop_ = deb_spect_m(jj)
          else
            marr_ = marr(j,jj)
            sarr_ = (marr_/(4.d0/3.d0*pi*rho_target))**(1.d0/3.d0)
            mpop_ = 0.d0
          endif
          write(15,*) marr_, sarr_, mpop_
        enddo

        close(15)

        endif  ! dbg

      endif

      return
      end


