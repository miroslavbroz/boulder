c chi2_func.f
c Run boulder and calculate chi^2 for a given starting vector x(ndim)
c Miroslav Broz (miroslav.broz@email.cz), Mar 10th 2011

      real*8 function chi2_func(x)

      include '../boulder/ucrm3.4.inc'
      include 'simplex.inc'
      include 'dependent.inc'
      real*8 x(ndim)
      include 'cb_sfd.inc'
      include 'cb_fam.inc'

c internals
      integer OBS_numb(OBS_MAX)
      real*8 OBS_diam(OBS_DATA_MAX, OBS_MAX)
      real*8 OBS_data(OBS_DATA_MAX, OBS_MAX)

      integer i,j,k,i1st,iu,k1,k2,bestk
      real*8 chi2,sigma2,tend,x_syn,y_syn,y_obs,eps,n_fams
      real*8 chi2_arr(SFD_MAX)
      integer n_SFD_arr(SFD_MAX)
      logical extra

c functions
      real*8 loginterp

c data
      data i1st /0/
      save i1st
      save OBS_numb, OBS_diam, OBS_data

c-----------------------------------------------------------------------

c get some parameters from the x() array
      j=0
      do i=1,nparam
        if (variable(i)) then
          j=j+1
          x_param(i) = x(j)
        endif
      enddo

      tend       = x_param(20)

c print and save debugging info
      if (debug) then
        iu=10
        open(unit=iu, file="chi2_func.tmp", access="append")

        write(*,30) "# x() array:"
30      format(a,$)
        do i = 1,nparam
          write(*,20) x_param(i)
20        format(1x,f10.5,$)
          write(iu,25) x_param(i)
25        format(1x,d27.20,$)        ! a maximal precision is needed
        enddo
        write(*,*)

        close(iu)
      endif

c read observed SFD's (only once)
      if (i1st.eq.0) then
        do i = 1,OBS_n
          OBS_numb(i) = OBS_DATA_MAX
          call read_ascii(OBS_file(i), OBS_numb(i), OBS_diam(1,i),
     :      OBS_data(1,i))

c sort ascending
          call sort2(OBS_numb(i), OBS_diam(1,i), OBS_data(1,i)) 
        enddo

        i1st = 1
      endif

c call boulder! (and save the result to arrays)
      call boulder_call(x_param, SFD_n, SFD_time, SFD_nbins, SFD_diam,
     :  SFD_data)

c what are the the indices corresponding to the time tend +- dt_chi2?
      k1 = 1
      eps = 1.d-8
      do while ((SFD_time(k1).lt.tend-dt_chi2-eps).and.(k1.lt.SFD_n))
        k1 = k1+1
      enddo
      k2 = k1
      do while ((SFD_time(k2).lt.tend+dt_chi2-eps).and.(k2.lt.SFD_n))
        k2 = k2+1
      enddo

c      write(*,*) '# SFD_time() = ', k1, k2, SFD_n, SFD_time(k1),
c     :  SFD_time(k2)	! dbg

c-----------------------------------------------------------------------
c
c calculate chi^2 for SFD's
c

      do k = k1, k2     ! scan a given RANGE of times!

        chi2 = 0.d0
        n_fit_SFD = 0

        do i = 1, OBS_n
          do j = 1, SFD_nbins(i,k)

            x_syn = SFD_diam(j,i,k)
            y_syn = SFD_data(j,i,k)

            if ((x_syn.ge.D_1(i)).and.(x_syn.lt.D_2(i))) then   ! take only bins within a prescribed range od D's

              y_obs = loginterp(OBS_diam(1,i),OBS_data(1,i),OBS_numb(i),
     :          x_syn,extra)                                    ! a simple linear interpolation is OK

              if (extra) then
                write(*,*) 'chi2_func(): extrapolation is NOT allowed!'
                stop
              endif

c              sigma2 = y_obs                     ! sigma = sqrt(y_obs) <- realistic sigma does NOT work!!!
              sigma2 = (y_obs/10.d0)**2.d0        ! THIS IS A PSEUDO-CHI^2!!! (actually, a 10 % relative error)

              chi2 = chi2 + (y_syn-y_obs)**2 / sigma2
              n_fit_SFD = n_fit_SFD + 1

            endif
          enddo   ! SFD_nbins
        enddo     ! OBS_n

        chi2_arr(k) = chi2
        n_SFD_arr(k) = n_fit_SFD

      enddo       ! SFD_time
c
c take the minimal chi^2 over the interval of times
c
      bestk = k1
      do k = k1,k2
        if (chi2_arr(k).lt.chi2_arr(bestk)) then
          bestk = k
        endif
c        write(*,*) '# chi2_arr(', k, '), SFD_time() = ', chi2_arr(k),
c     :    SFD_time(k)	! dbg
      enddo

      chi2 = chi2_arr(bestk)
      n_fit_SFD = n_SFD_arr(bestk)
      time_chi2 = SFD_time(bestk)

      chi2_SFD = w_sfd*chi2
c
c finally, create a datafile suitable for Gnuplot
c
      if (debug) then
        iu=10
        open(unit=iu,file="chi2_SFD.dat",status="unknown")

        k = bestk
        do i = 1, OBS_n
          do j = 1, SFD_nbins(i,k)

            x_syn = SFD_diam(j,i,k)
            y_syn = SFD_data(j,i,k)

            if ((x_syn.ge.D_1(i)).and.(x_syn.lt.D_2(i))) then

              y_obs = loginterp(OBS_diam(1,i),OBS_data(1,i),OBS_numb(i),
     :          x_syn,extra)
              sigma2 = (y_obs/10.d0)**2.d0

              write(iu,*) x_syn,y_syn,y_obs,sqrt(abs(sigma2)),i
            else
              write(iu,*) x_syn,y_syn,'?	?',i       ! plot also data points which are NOT used
            endif

          enddo   ! SFD_nbins
        enddo     ! OBS_n

        close(iu)
      endif       ! debug

c-----------------------------------------------------------------------
c
c chi^2 for the number of families (which fullfil given criteria)
c
      n_fit_FAM = 0
      chi2 = 0.d0
      chi2_FAM = 0.d0
      k = 1

      if (debug) then
        iu=10
        open(unit=iu,file="chi2_FAM.dat",status="unknown")
      endif

      do i = 1, OFM_n

c count the families in boulder output (and check criteria)
        n_fams = 0.d0
        do j = 1, FAM_n
          if ( (FAM_pop(j).eq.OFM_pop(i))
     :      .and.(FAM_dtarget(j).ge.OFM_dpb(i))
     :      .and.(FAM_dlf(j).ge.OFM_dlf(i))
     :      .and.(FAM_dproject(j).ge.OFM_dproject(i))
     :      .and.(FAM_m2nd_over_mtot(j).ge.OFM_lfpb_min(i))
     :      .and.(FAM_m2nd_over_mtot(j).le.OFM_lfpb_max(i)) ) then

            n_fams = n_fams + FAM_ncoll_kl(j)
          endif
        enddo

        sigma2 = OFM_numb(i)                               ! sigma = sqrt(n)          
        if (sigma2.lt.1.d0) sigma2 = 1.d0

        chi2 = (n_fams - OFM_numb(i))**2 / sigma2

        n_fit_FAM = n_fit_FAM + 1
        chi2_FAM = chi2_FAM + w_fams(i)*chi2

        if (debug) then
          write(iu,*) i,n_fams,OFM_numb(i),sqrt(abs(sigma2)),chi2
        endif

      enddo

      if (debug) then
        close(iu)
      endif
c
c output individual families too
c
      if (debug) then
        iu=10
        open(unit=iu,file="families.out",status="unknown")

        do j = 1, FAM_n
          write(iu,60)
     :      j,
     :      FAM_pop(j),
     :      FAM_time(j),
     :      FAM_dtarget(j),
     :      FAM_dproject(j),
     :      FAM_ncoll_kl(j),
     :      FAM_vcoll(j),
     :      FAM_dlr(j),
     :      FAM_dlf(j),
     :      FAM_m2nd_over_mtot(j)
60        format(i8.8,1x,i4,1x,f12.4,2(1x,f12.3),1x,f8.0,1x,f12.4,
     :      8x,2(1x,f12.4),1x,f12.6)
        enddo

        close(iu)
      endif     ! debug

c-----------------------------------------------------------------------

c sum all things up

      chi2 = chi2_SFD + chi2_FAM
      n_fit = n_fit_SFD + n_fit_FAM

c print and save debugging info
      if (debug) then
        write(*,30) "# chi^2 value: "
        write(*,40) time_chi2,
     :    n_fit_SFD, n_fit_FAM, n_fit,
     :    chi2_SFD, chi2_FAM, chi2
40      format(218x,f8.2,1x,3(i5,1x),3(f8.2,1x))

        open(unit=iu, file="chi2_func.tmp", access="append")
        write(iu,50) time_chi2,
     :    n_fit_SFD, n_fit_FAM, n_fit,
     :    chi2_SFD, chi2_FAM, chi2
50      format(5x,f8.2,1x,3(i5,1x),3(f8.2,1x))
        close(iu)
      endif

      chi2_func = chi2

      return
      end


