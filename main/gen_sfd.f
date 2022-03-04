c gen_sfd.f
c Generate a DIFFERENTIAL size-frequency distribution.
c original version by David Nesvorny, modified by MB

      program gen_sfd

      include '../boulder/ucrm3.4.inc'

      integer i,ibin(BINMAX),nbin
      real*8 rho1,dfactor,onethird,mass_earth,rfactor
      real*8 rad(BINMAX),numb,cumul
      real*8 r1,r2,qa,qb,qc,rmax,r_norm,n_norm,rmin
      real*8 lrbin,lrbin_start,lrbin_end,nstart,nend
      real*8 rstart,rend
      real*8 numb_int(BINMAX),ncum1(BINMAX),ncum2(BINMAX)
      real*8 dpi,totm,dm
      real*8 d1,d2,dmin,dmax,d_norm

      integer iargc
      character*256 str

c... Boulder's mfactor
c      mfactor = 2.0d0

c... Constants
      onethird = 1.d0/3.d0
      mass_earth = 5.9742d27
      dpi = acos(-1.d0)

c... Read input parameters
      if (iargc().eq.0) then
        read(*,*,end=990,err=990) rho1
        read(*,*,end=990,err=990) d1
        read(*,*,end=990,err=990) d2
        read(*,*,end=990,err=990) qa
        read(*,*,end=990,err=990) qb
        read(*,*,end=990,err=990) qc
        read(*,*,end=990,err=990) dmin
        read(*,*,end=990,err=990) dmax
        read(*,*,end=990,err=990) d_norm
        read(*,*,end=990,err=990) n_norm
        read(*,*,end=990,err=990) mfactor
      else if (iargc().eq.11) then
        call getarg(1, str)
        read(str,*,err=990,end=990) rho1
        call getarg(2, str)
        read(str,*,err=990,end=990) d1
        call getarg(3, str)
        read(str,*,err=990,end=990) d2
        call getarg(4, str)
        read(str,*,err=990,end=990) qa
        call getarg(5, str)
        read(str,*,err=990,end=990) qb
        call getarg(6, str)
        read(str,*,err=990,end=990) qc
        call getarg(7, str)
        read(str,*,err=990,end=990) dmin
        call getarg(8, str)
        read(str,*,err=990,end=990) dmax
        call getarg(9, str)
        read(str,*,err=990,end=990) d_norm
        call getarg(10, str)
        read(str,*,err=990,end=990) n_norm
        call getarg(11, str)
        read(str,*,err=990,end=990) mfactor
      else
        write(*,*) 'Usage: gen_sfd  rho d1 r2 qa qb qc ',
     :    'dmin dmax d_norm n_norm boulder_mfactor'
        stop
      endif

c... Input differential SFD
       r1 = d1/2.d0
       r2 = d2/2.d0
       rmin = dmin/2.d0
       rmax = dmax/2.d0
       r_norm = d_norm/2.d0

c      rho = 1.d0
c      r1 = 2.d0
c      r2 = 0.05d0
c      qa = -5.d0
c      qb = -1.d0
c      qc = -3.5d0
c      rmax = 500.d0
c      r_norm = 50.d0
c      n_norm = 45.0d4
    
c... Check on things
      if(r_norm.lt.r1) then
         write(*,*)'Error: r_norm must be .ge. r1. Exiting...'
         stop
      end if

c      nbin = 48
c      rfactor = rfactor * 48.d0/dble(nbin)
c      rmax = rmin * (10**rfactor)**nbin

c	rmin = 0.0005d0  ! km - now it is read from input
c	rmax = 40000.d0	 ! km	
c      rmin = 0.05d0     ! km 
c      rmax = 4000.d0    ! km
      rfactor = log10(mfactor**onethird)
      nbin = int( (log10(rmax)-log10(rmin)) / rfactor )

      if (nbin.gt.BINMAX) then
        write(*,*) 'Error: nbin must be .le. BINMAX. Exiting...'
        stop
      end if

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
         dm = numb_int(i)*(4.d0/3.d0)*dpi*((rad(i)*1.d5)**3)*rho1
         totm = totm + dm 
      end do
      totm=totm/mass_earth

c...  Write things out 
      write(*,8848)nbin,rho1,rad(1),totm
 8848 format(i5,2(1x,f15.7),1x,e15.7)
      do i=1,nbin
c         write(*,*)rad(i),ncum1(i),ncum2(i)
         write(*,*)i,rad(i),numb_int(i)
      end do

      stop

c... Error handlers
990   continue
      write(*,*) '# Error reading standard input.'
      stop

      end


