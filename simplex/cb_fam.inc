c cb_fam.inc
c Arrays suitable for the output of families.
c Miroslav Broz (miroslav.broz@email.cz), Sep 20th 2012

      integer FAM_MAX
      parameter(FAM_MAX = 10000)		! max number of families to save

      integer FAM_n
      integer FAM_pop(FAM_MAX)
      real*8 FAM_time(FAM_MAX)
      real*8 FAM_dtarget(FAM_MAX)
      real*8 FAM_dproject(FAM_MAX)
      real*8 FAM_ncoll_kl(FAM_MAX)
      real*8 FAM_vcoll(FAM_MAX)
      real*8 FAM_dlr(FAM_MAX)
      real*8 FAM_dlf(FAM_MAX)
      real*8 FAM_m2nd_over_mtot(FAM_MAX)

      common /cb_fam/ FAM_time, FAM_dtarget, FAM_dproject, FAM_ncoll_kl,
     :  FAM_vcoll, FAM_dlr, FAM_dlf, FAM_m2nd_over_mtot, FAM_n, FAM_pop


