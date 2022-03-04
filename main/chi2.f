c chi2.f
c Calculate chi^2 for a single boulder run.
c Miroslav Broz (miroslav.broz@email.cz), Mar 4th 2011

      program chi2_

      include '../boulder/ucrm3.4.inc'
      include '../simplex/simplex.inc'
      include '../simplex/dependent.inc'

      real*8 x(ndimmax),chi2,probp,probq
      integer i
      character*80 str

      real*8 chi2_func,gammp,gammq

c  read input parameters
5     continue
        read(*,10,err=990,end=990) str
10      format(a)
      if (str(1:1).eq.'#') goto 5
      read(str,*,err=990,end=990) nparam

      read(*,*,err=990,end=990) (x_param(i), i=1,nparam)

      read(*,*,err=990,end=990) OBS_n
      if (OBS_n.gt.OBS_MAX) then
        write(*,*) 'Error: OBS_n > OBS_MAX = ', OBS_MAX
        stop
      endif
      do i = 1,OBS_n
        read(*,10,err=990,end=990) OBS_file(i)
      enddo
      read(*,*,err=990,end=990) w_sfd

      read(*,*,err=990,end=990) OFM_n
      if (OFM_n.gt.OFM_MAX) then
        write(*,*) 'Error: OFM_n > OFM_MAX = ', OFM_MAX
        stop
      endif
      do i = 1,OFM_n
        read(*,*,err=990,end=990)
     :    OFM_pop(i),
     :    OFM_numb(i),
     :    OFM_dpb(i),
     :    OFM_dlf(i),
     :    OFM_dproject(i),
     :    OFM_lfpb_min(i),
     :    OFM_lfpb_max(i),
     :    w_fams(i)
      enddo

      do i = 1,OBS_n
        read(*,*,err=990,end=990) D_1(i), D_2(i)
      enddo
      read(*,*,err=990,end=990) n_chi2
      read(*,*,err=990,end=990) dt_chi2
      read(*,*,err=990,end=990) Pint_tdepend
      read(*,*,err=990,end=990) mfactor
      read(*,*,err=990,end=990)
     :  family_lfpb_min,
     :  family_lfpb_max,
     :  family_dpb_treshold,
     :  family_dlf_treshold
      read(*,*,err=990,end=990) n_ceres_vesta
      if (n_ceres_vesta.gt.CERES_VESTA_MAX) then
        write(*,*) 'Error: n_ceres_vesta > CERES_VESTA_MAX = ',
     :    CERES_VESTA_MAX
        stop
      endif
      do i = 1,n_ceres_vesta
        read(*,*,err=990,end=990) ceres_vesta_diam(i),
     :    ceres_vesta_pop(i)
      enddo
      read(*,*,err=990,end=990) debug
      read(*,*,err=990,end=990) debug_boulder

c  write input parameters
      if (debug) then
        write(*,*) '# nparam = ', nparam
        do i = 1,nparam
          write(*,*) "# x_param(", i, ") = ", x_param(i)
        enddo
        write(*,*) '# OBS_n = ', OBS_n
        do i = 1,OBS_n
          write(*,*) '# OBS_file(', i, ') = ', OBS_file(i)
        enddo
        write(*,*) '# w_sfd = ', w_sfd
        write(*,*) '# OFM_n = ', OFM_n
        do i = 1,OFM_n
          write(*,*) '# OFM_pop(', i, ') = ', OFM_pop(i)
          write(*,*) '# OFM_numb(', i, ') = ', OFM_numb(i)
          write(*,*) '# OFM_dpb(', i, ') = ', OFM_dpb(i), ' km'
          write(*,*) '# OFM_dlf(', i, ') = ', OFM_dlf(i), ' km'
          write(*,*) '# OFM_dproject(', i, ') = ', OFM_dproject(i),' km'
          write(*,*) '# OFM_lfpb_min(', i, ') = ', OFM_lfpb_min(i)
          write(*,*) '# OFM_lfpb_max(', i, ') = ', OFM_lfpb_max(i)
          write(*,*) '# w_fams(', i, ') = ', w_fams(i)
        enddo
        do i = 1,OBS_n
          write(*,*) '# D_1(', i, ') = ', D_1(i)
          write(*,*) '# D_2(', i, ') = ', D_2(i)
        enddo
        write(*,*) '# n_chi2 = ', n_chi2
        write(*,*) '# dt_chi2 = ', dt_chi2
        write(*,*) '# Pint_tdepend = ', Pint_tdepend
        write(*,*) '# mfactor = ', mfactor
        write(*,*) '# family_lfpb_min = ', family_lfpb_min
        write(*,*) '# family_lfpb_max = ', family_lfpb_max
        write(*,*) '# family_dpb_treshold = ', family_dpb_treshold,' cm'
        write(*,*) '# family_dlf_treshold = ', family_dlf_treshold,' cm'
        write(*,*) '# n_ceres_vesta = ', n_ceres_vesta
        do i = 1,n_ceres_vesta
          write(*,*) '# ceres_vesta_diam(', i, ') = ',
     :      ceres_vesta_diam(i)
          write(*,*) '# ceres_vesta_pop(', i, ') = ',
     :      ceres_vesta_pop(i)
        enddo
        write(*,*) '# debug = ', debug
        write(*,*) '# debug_boulder = ', debug_boulder
      endif

c  calculate chi^2 and the corresponding probability
      ndim = nparam
      do i = 1, nparam
        x(i) = x_param(i)
        variable(i) = .TRUE.
      enddo

      chi2 = chi2_func(x)

      probp = gammp(dble(n_fit)/2.d0,chi2/2.d0)
      probq = gammq(dble(n_fit)/2.d0,chi2/2.d0)

c  write output
      if (debug) then
        write(*,*) '# chi2	probp	probq	nfit	',
     :  'n_SFD	n_FAM	',
     :  'chi2_SFD	chi2_FAM	',
     :  'time_chi2 [Myr]'
      endif
      write(*,*) real(chi2), real(probp), real(probq), n_fit,
     :  n_fit_SFD, n_fit_FAM,
     :  real(chi2_SFD), real(chi2_FAM),
     :  real(time_chi2)

      stop

c  error handlers
990   continue
      write(*,*) 'chi2: Error reading standard input.'

      end


