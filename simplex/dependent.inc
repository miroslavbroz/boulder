c dependent.inc
c Include file for chi2 and chi2_func with dependent parameters
c   (i.e., not passed as the X array, but in a common block).
c Miroslav Broz (miroslav.broz@email.cz), Dec 7th 2009

c general parameters of the chi^2 calculation
      integer ndim, nparam
      integer n_fit, n_fit_SFD, n_fit_FAM
      real*8 chi2_SFD, chi2_FAM

c parameters related to the OBServed SFD's
      integer n_chi2
      real*8 dt_chi2, time_chi2
      real*8 w_sfd
      integer OBS_n
      character*80 OBS_file(OBS_MAX)
      real*8 D_1(OBS_MAX), D_2(OBS_MAX)

c parameters related to the Observed FaMilies
      integer OFM_n
      integer OFM_pop(OFM_MAX)
      real*8 OFM_numb(OFM_MAX)
      real*8 OFM_dpb(OFM_MAX)
      real*8 OFM_dlf(OFM_MAX)
      real*8 OFM_dproject(OFM_MAX)
      real*8 OFM_lfpb_min(OFM_MAX)
      real*8 OFM_lfpb_max(OFM_MAX)
      real*8 w_fams(OFM_MAX)

      real*8 x_param(NDIMMAX)
      logical variable(NDIMMAX)
      logical Pint_tdepend
      integer n_ceres_vesta
      real*8 ceres_vesta_diam(CERES_VESTA_MAX)
      integer ceres_vesta_pop(CERES_VESTA_MAX)
      logical debug, debug_boulder

      common /dependent/
     :  chi2_SFD, chi2_FAM,
     :  dt_chi2, time_chi2,
     :  w_sfd,
     :  OBS_file,
     :  D_1, D_2,
     :  OFM_pop,
     :  OFM_numb,
     :  OFM_dpb,
     :  OFM_dlf,
     :  OFM_dproject,
     :  OFM_lfpb_min,
     :  OFM_lfpb_max,
     :  w_fams,
     :  x_param,
     :  ceres_vesta_diam,
     :  ceres_vesta_pop,
     :  ndim, nparam,
     :  n_fit, n_fit_SFD, n_fit_FAM,
     :  n_chi2,
     :  OBS_n,
     :  OFM_n,
     :  variable,
     :  Pint_tdepend,
     :  n_ceres_vesta,
     :  debug, debug_boulder

