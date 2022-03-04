c add_ceres_vesta.f
c Add (1) Ceres and (4) Vesta to the SFD.
c Miroslav Broz (miroslav.broz@email.cz), Sep 20th 2011

      subroutine add_ceres_vesta(n, diam, pop, nbins, marr, mpop, sarr,
     :  npop)

      include '../boulder/ucrm3.4.inc'   ! all parameters file
     
      integer n
      real*8 diam(n)
      integer pop(n)
      integer nbins
      real*8 rho_bulk 
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX) 
      real*8 mpop(Manuli,0:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)

      integer i,j,jj
      real*8 logrfactor,radius,tmp

      logrfactor = 0.5d0*log10(mfactor**(1.d0/3.d0))

      do i=1,n
        radius = diam(i)*1.d5/2.d0    ! cm
        j = pop(i)                    ! to which population to add the body
        do jj=1,nbins
          tmp = abs(log10(sarr(j,jj))-log10(radius))
          if (tmp.lt.logrfactor) then
            mpop(j,jj) = mpop(j,jj) + marr(j,jj)
            npop(j,jj,1) = npop(j,jj,1) + 1.d0

            if (dbg) then
              write(*,*) '# added a ', diam(i), ' km body:	', 
     :          marr(j,jj), sarr(j,jj), mpop(j,jj)
            endif

          endif
        enddo
      enddo

      return
      end

