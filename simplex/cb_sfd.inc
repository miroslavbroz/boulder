c cb_sfd.inc
c Arrays suitable for boulder output storage.
c Miroslav Broz (miroslav.broz@email.cz), Mar 7th 2011

      integer SFD_MAX
      parameter(SFD_MAX = 1000)			! maximum number of SFD's to save
      integer SFD_n
      integer SFD_nbins(Manuli,SFD_MAX)
      real*8 SFD_time(SFD_MAX)
      real*8 SFD_diam(0:BINMAX,Manuli,SFD_MAX)
      real*8 SFD_data(0:BINMAX,Manuli,SFD_MAX)

