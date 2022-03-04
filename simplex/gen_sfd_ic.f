c gen_sfd_ic.f
c Generate SFD's and prepare boulder input arrays - based on gen_sfd and gen_ic programs (by David Nesvorny).
c modified by Miroslav Broz (miroslav.broz@email.cz), Mar 10th 2011

      subroutine gen_sfd_ic(rho_bulk, d1, d2, qa, qb, qc,
     :  dmax, d_norm, n_norm,
     :  j,axe,delta_a,nbins,marr,mpop,sarr,ecc,inc,npop)

      include '../boulder/ucrm3.4.inc'   ! all parameters file

      real*8 rho_bulk, d1, d2, qa, qb, qc, dmax, d_norm, n_norm
      integer j
      real*8 axe(Manuli),delta_a(Manuli)
      integer nbins(Manuli)
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX) 
      real*8 mpop(Manuli,0:BINMAX)
      real*8 ecc(Manuli,BINNEG:BINMAX),inc(Manuli,BINNEG:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata) ! the Ndata vector is for additional data

      integer i,ibin(BINMAX),nbin
      real*8 dfactor,onethird,mass_earth,rfactor
      real*8 rad(BINMAX),numb,cumul
      real*8 r1,r2,rmax,r_norm,rmin
      real*8 lrbin,lrbin_start,lrbin_end,nstart,nend
      real*8 rstart,rend
      real*8 numb_int(BINMAX),ncum1(BINMAX),ncum2(BINMAX)
      real*8 dpi,totm,dm

      integer jj
      real*8 nnx,nnall,num1,num2
      real*8 trash_bin_marr
      integer nbinneg(Manuli)

c---------------------------------------------------------------

c gen_sfd part

c... Constants
      onethird = 1.d0/3.d0
      mass_earth = 5.9742d27
      dpi = acos(-1.d0)

c... Input differential SFD
      r1 = d1/2.d0
      r2 = d2/2.d0
      rmax = dmax/2.d0
      r_norm = d_norm/2.d0

c... Check on things
      if(r_norm.lt.r1) then
         write(*,*)'# gen_sfd_ic: r_norm must be .ge. r1. Exiting...'
         stop
      end if

      rfactor = log10(mfactor**onethird)
      rmin = 0.005d0     ! km
      rmax = 4000.d0    ! km
      nbin = int( (log10(rmax)-log10(rmin)) / rfactor )

c... Change things to cumulative
      qa = qa + 1.d0
      qb = qb + 1.d0  
      qc = qc + 1.d0

c...  Calculate number of objects in each bin
      do i=1,nbin
         lrbin = log10(rmin) + dble(i-1)*rfactor
         lrbin_start = lrbin - 0.5d0*rfactor
         lrbin_end = lrbin + 0.5d0*rfactor

         rad(i) = 10.d0**lrbin 
         rstart = 10.d0**lrbin_start
         rend = 10.d0**lrbin_end

         nstart = cumul(rmax,r1,r2,qa,qb,qc,r_norm,n_norm,rstart)
         nend = cumul(rmax,r1,r2,qa,qb,qc,r_norm,n_norm,rend)
         numb = nstart - nend
c...  Round things
         numb_int(i) = dnint(numb)
      end do

c...  Calculate cumulative
      ncum1(nbin) = numb_int(nbin)
      do i=nbin-1,1,-1
         ncum1(i) = ncum1(i+1) + numb_int(i)
      end do
      do i=1,nbin
         ncum2(i) = cumul(rmax,r1,r2,qa,qb,qc,r_norm,n_norm,rad(i))
      end do

c... Calculate total mass
      totm=0.d0
      do i=1,nbin
         dm = numb_int(i)*(4.d0/3.d0)*dpi*((rad(i)*1.d5)**3)*rho_bulk
         totm = totm + dm 
      end do
      totm=totm/mass_earth

c gen_sfd.f:
c      write(*,8848)nbin,rho,rad(1),totm
c 8848 format(i5,3(1x,f15.7))
c      do i=1,nbin
c         write(*,*)i,rad(i),numb_int(i)
c      end do

c---------------------------------------------------------------

c gen_ic part
      
c zeroing of all arrays (neccesarry for repeated calls by simplex)
      do jj=BINNEG,BINMAX
        sarr(j,jj)=0.d0
        marr(j,jj)=0.d0
        npop(j,jj,1)=0.0d0
      enddo
      do jj=0,BINMAX
        mpop(j,jj)=0.d0
      enddo

c... Read, convert and write into Boulder input file
      trash_bin_marr=4.d0/3.d0*dpi*rho_bulk*(rmin*1.d5)**3 / mfactor
      marr(j,0) = trash_bin_marr
      mpop(j,0) = 0.d0

      nbins(j) = nbin
      do jj=1,nbins(j)
        sarr(j,jj)=rad(jj)*1.d5 
        marr(j,jj)=4.d0/3.d0*dpi*rho_bulk*sarr(j,jj)**3 
        mpop(j,jj)=marr(j,jj) * numb_int(jj)
        ecc(j,jj) = 0.1d0    ! Dummy eccentricity and inclination; keep at 0.1 or weird things happen
        inc(j,jj) = 0.1d0
        axe(j) = 42.d0
        delta_a(j) = 0.01d0
        npop(j,jj,1)=mpop(j,jj)/marr(j,jj)

c        write(*,*) j,jj,marr(j,jj),sarr(j,jj),mpop(j,jj)	! dbg
      enddo
c      write(*,*)	! dbg

      return
      end

c gen_ic.f:
c      open(unit=31,file=boulder_name, status = 'unknown')
c      write(31,*) Nanuli
c      do j=1,Nanuli
c         open(11,file=filename(j),status='old')
c         read(11,*)nbins(j),rho(j),rmin
c         trash_bin_marr=4.d0/3.d0*dpi*rho(j)*(rmin*1.d5)**3 / mfactor
c         write(31,8488)axe,daxe,nbins(j),trash_bin_marr,0.
c 8488    format(2(1x,f15.7),1x,i5,2(1x,e15.7))
c         do jj=1,nbins (j)
c            read(11,*)iden,radius,nnall
c            sarr(j,jj)=radius*1.d5 
c            marr(j,jj)=4.d0/3.d0*dpi*rho(j)*sarr(j,jj)**3 
c            mpop(j,jj)=marr(j,jj) * nnall 
c            write(31,*)marr(j,jj),sarr(j,jj),mpop(j,jj)
c            write(31,*) 0.1,0.1    ! Dummy eccentricity and inclination; keep at 0.1 or weird things happen
c         enddo
c         close(11)
c      enddo
c      close(31)

c read_ib.f:
c      open(unit=1,file=boulder_name,status='old', iostat = ios)
c      if (ios .ne. 0) then
c       write (6, *) ' No input SFD file, error ', ios, ' Try again '
c       stop
c      endif
c      read(1,*)Nanuli
c      do i=1,Nanuli
c       read(1,*)axe(i),delta_a(i),nbins(i),marr(i,0),mpop(i,0)
c       do jj=1,nbins(i)
c        read(1,*)marr(i,jj),sarr(i,jj),mpop(i,jj)
c        read(1,*)ecc(i,jj),inc(i,jj)
c        npop(i,jj,1)=mpop(i,jj)/marr(i,jj)
c       enddo
c       ecc(i,0)=ecc(i,1)
c       inc(i,0)=inc(i,1)
c      enddo
c      close(1)

