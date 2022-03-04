c simplex.f
c Run simplex to optimize SFD's with boulder.
c Miroslav Broz (miroslav.broz@email.cz), Sep 21st 2011

      program simplex

      include '../boulder/ucrm3.4.inc'
      include '../simplex/simplex.inc'
      include '../simplex/dependent.inc'
      include '../simplex/cb_itmax.inc'

      integer i,j,iter,mp,np
      integer id(ndimmax+1)
      real*8 x(ndimmax), e(ndimmax)
      real*8 p(ndimmax+1,ndimmax), y(ndimmax+1),
     :  ftol, xtry(ndimmax)
      real*8 e_param(ndimmax)
      character*80 str
      real*8 chi2_func

c  read input parameters
5     continue
        read(*,10,err=990,end=990) str
10      format(a)
      if (str(1:1).eq.'#') goto 5
      read(str,*,err=990,end=990) nparam

      read(*,*,err=990,end=990) (x_param(i), i=1,nparam)
      read(*,*,err=990,end=990) (e_param(i), i=1,nparam)
      read(*,*,err=990,end=990) (variable(i), i=1,nparam)

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
      read(*,*,err=990,end=990) ftol
      read(*,*,err=990,end=990) ITMAX
      read(*,*,err=990,end=990) debug
      read(*,*,err=990,end=990) debug_boulder

c  write input parameters
      if (debug) then
        write(*,*) '# nparam = ', nparam
        do i = 1,nparam
          write(*,*) "# x_param(", i, ") = ", x_param(i)
        enddo
        do i = 1,nparam
          write(*,*) "# e_param(", i, ") = ", e_param(i)
        enddo
        do i = 1,nparam
          write(*,*) "# variable(", i, ") = ", variable(i)
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
          write(*,*) '# OFM_dpb(', i, ') = ', OFM_dpb(i)
          write(*,*) '# OFM_dlf(', i, ') = ', OFM_dlf(i)
          write(*,*) '# OFM_dproject(', i, ') = ', OFM_dproject(i)
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
        write(*,*) '# family_dpb_treshold = ', family_dpb_treshold
        write(*,*) '# family_dlf_treshold = ', family_dlf_treshold
        write(*,*) '# n_ceres_vesta = ', n_ceres_vesta
        do i = 1,n_ceres_vesta
          write(*,*) '# ceres_vesta_diam(', i, ') = ',
     :      ceres_vesta_diam(i)
          write(*,*) '# ceres_vesta_pop(', i, ') = ',
     :      ceres_vesta_pop(i)
        enddo
        write(*,*) '# ftol = ', ftol
        write(*,*) '# ITMAX = ', ITMAX
        write(*,*) '# debug = ', debug
        write(*,*) '# debug_boulder = ', debug_boulder
      endif

c  resolve variable/fixed parameters
      ndim = 0
      do i = 1,nparam
        if (variable(i)) then
          ndim = ndim + 1
          x(ndim) = x_param(i)
          e(ndim) = e_param(i)
        endif
      enddo

      if (debug) then
        write(*,*) '# ndim = ', ndim
      endif

c  pass the rest in common block /dependent/ 

c  initialise the simplex (p array)
      do i = 1,ndim+1
        do j = 1,ndim
          if (j.eq.i) then
            p(i,j) = x(j) + e(i)
          else
            p(i,j) = x(j)
          endif
        enddo
      enddo

      if (debug) then
        write(*,*) "# initial p() array:"
        do i = 1,ndim+1
          write(*,30)
30        format('# ',$)
          do j = 1,ndim
            write(*,20) p(i,j)
20          format(1x,d27.20,$)
          enddo
          write(*,*)
        enddo
      endif

c  y() array
      do i = 1,ndim+1
        do j = 1,ndim
          xtry(j) = p(i,j)
        enddo
        y(i) = chi2_func(xtry)
      enddo

      if (debug) then
        write(*,*) "# initial y() array:"
        do i = 1,ndim+1
          write(*,30)
          write(*,*) y(i)
        enddo
      endif

      mp = ndimmax+1
      np = ndimmax
      iter = 0

c  run it!
      call amoeba(p,y,mp,np,ndim,ftol,chi2_func,iter)

c  sort the output according to chi^2
      call srtidx(ndim+1,y,id)

c  write the result of minimalisation

      if (debug) then
        write(*,*) "# iter = ", iter

        write(*,*) "# p() array:"
        do i = 1,ndim+1
          write(*,30)
          do j = 1,ndim
            write(*,20) p(id(i),j)
          enddo
          write(*,*)
        enddo

        write(*,*) "# y() array:"
        do i = 1,ndim+1
          write(*,30)
          write(*,*) y(id(i))
        enddo
      endif

      j = 0
      do i = 1,nparam
        if (variable(i)) then
          j = j+1
          write(*,20) p(id(1),j)
        else
          write(*,20) x_param(i)
        endif
      enddo
      write(*,*) y(id(1))

      stop

c  error handlers
990   continue
      write(*,*) '# simplex: Error reading standard input.'

      end

