c simplex.inc
c Include file for simplex with max dimensions of arrays.
c Miroslav Broz (miroslav.broz@email.cz), Dec 7th 2009

      integer NDIMMAX
      parameter(NDIMMAX=20)		! maximum number of dimensions

      integer OBS_MAX, OBS_DATA_MAX, OFM_MAX, CERES_VESTA_MAX
      parameter(OBS_MAX = 5)		! max number of observed SFD's
      parameter(OBS_DATA_MAX = 10000)	! max number of data points in the observed SFD file
      parameter(OFM_MAX = 10)		! max number of Observed FaMilies (criteria)
      parameter(CERES_VESTA_MAX = 100)	! max number of bodies to add to pop1 SFD

