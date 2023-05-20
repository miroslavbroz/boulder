c pop_merge.f
c Merge 1 bin with other bins.
c Miroslav Broz (miroslav.broz@email.cz), Jul 11th 2022

      subroutine pop_merge(k,kk,nbins,marr,sarr,npop,mpop,m,n)

      include 'ucrm3.4.inc'

      integer k,kk
      integer nbins(Manuli)
      real*8 marr(Manuli,BINNEG:BINMAX)
      real*8 sarr(Manuli,BINNEG:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 mpop(Manuli,0:BINMAX)
      real*8 m,n

      real*8 mtot

c add bin if needed
      if (m.gt.marr(k,nbins(k))) then
        if (nbins(k).ge.BINMAX) then
          write(*,*) 'Error: BINMAX = ',BINMAX,' too small'
          stop
        endif
        nbins(k) = nbins(k)+1
        kk = nbins(k)
        marr(k,kk) = m
        npop(k,kk,1) = n
        mpop(k,kk) = m*n
        return
      endif

c find the closest bin
      do while (m.ge.marr(k,kk))
        kk = kk+1
      enddo
      if (m/marr(k,kk-1).lt.marr(k,kk)/m) then
        kk = kk-1
      endif

c merge bins
      mtot = npop(k,kk,1)*marr(k,kk) + m*n
      npop(k,kk,1) = npop(k,kk,1) + n
      marr(k,kk) = mtot/npop(k,kk,1)
      sarr(k,kk) = (marr(k,kk)/(4.d0/3.d0*pi*rho(k)))**(1.d0/3.d0)
      if (kk.ge.0) mpop(k,kk) = mtot

      return
      end


