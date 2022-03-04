c save_ob.f
c Save boulder output to the SFD_data array.
c Miroslav Broz (miroslav.broz@email.cz), Mar 7th 2011

      subroutine save_ob(Nanuli,t,nbins,marr,sarr,mpop,
     :   SFD_n,SFD_time,SFD_nbins,SFD_diam,SFD_data)

      include '../boulder/ucrm3.4.inc'
      integer Nanuli, nbins(Manuli)
      real*8 t
      real*8 marr(Manuli,BINNEG:BINMAX), sarr(Manuli,BINNEG:BINMAX) 
      real*8 mpop(Manuli,0:BINMAX)
      include 'cb_sfd.inc'

      integer j,jj
      real*8 N_gt_D

      SFD_n = SFD_n + 1
      if (SFD_n.gt.SFD_MAX) then
        write(*,*) '# save_ob(): SFD_n .gt. SFD_MAX = ', SFD_MAX
        stop
      endif

c get differential SFD's from marr,sarr,mpop arrays
      SFD_time(SFD_n) = t/1.d6                          ! yr -> Myr
      do j=1,Nanuli
        SFD_nbins(j,SFD_n) = nbins(j)
        do jj=1,nbins(j)
          SFD_diam(jj,j,SFD_n) = 2.d0*sarr(j,jj)/1.d5   ! cm -> km
          SFD_data(jj,j,SFD_n) = mpop(j,jj)/marr(j,jj)  ! total mass/one body mass
        enddo
      enddo

c convert them to cumulative SFD's
      do j=1,Nanuli
        N_gt_D = 0.d0
        do jj=nbins(j),1,-1
          N_gt_D = N_gt_D + SFD_data(jj,j,SFD_n)
          SFD_data(jj,j,SFD_n) = N_gt_D
        enddo
      enddo

      return
      end


