c gen_ic.f
c Generate ib.data file (i.e. input for boulder) from gen_sfd.out.
c original version by David Nesvorny, modified by MB

      program gen_ic
      include '../boulder/ucrm3.4.inc'   ! all parameters file

      integer j,jj
      real*8 radius,iden,nnx,nnall,num1,num2
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX) 
      real*8 mpop(Manuli,0:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata) ! the Ndata vector is for additional data
      real*8 ecc(Manuli,BINNEG:BINMAX),inc(Manuli,BINNEG:BINMAX)
      real*8 axe,daxe,trash_bin_marr,dpi,rmin
      integer Nanuli,nbins(Manuli),nbinneg(Manuli),ios
      character*80 filename(Manuli)
      character*80 boulder_name

c... Constants (not used in the fragmentation version)
      axe = 42.d0     ! this one IS used in ../boulder/focussing_factor.f !
      daxe = 0.01d0
      
      dpi = acos(-1.d0)

c... Read file names and axe
      read(5,*) Nanuli
      do j=1,Nanuli
         read(5,*) filename(j)
      end do
      read(5,*) boulder_name
      read(5,*) mfactor
c      read(5,*,iostat=ios) axe

c... Read, convert and write into Boulder input file
      open(unit=31,file=boulder_name, status = 'unknown')
      write(31,*) Nanuli
      do j=1,Nanuli
         open(11,file=filename(j),status='old')
         read(11,*)nbins(j),rho(j),rmin
         trash_bin_marr=4.d0/3.d0*dpi*rho(j)*(rmin*1.d5)**3 / mfactor
         write(31,8488)axe,daxe,nbins(j),trash_bin_marr,0.
 8488    format(2(1x,f15.7),1x,i5,2(1x,e22.16))
         do jj=1,nbins (j)
            read(11,*)iden,radius,nnall
            sarr(j,jj)=radius*1.d5 
            marr(j,jj)=4.d0/3.d0*dpi*rho(j)*sarr(j,jj)**3 
            mpop(j,jj)=marr(j,jj) * nnall 
            write(31,*)marr(j,jj),sarr(j,jj),mpop(j,jj)
            write(31,*) 0.1,0.1    ! Dummy eccentricity and inclination; keep at 0.1 or weird things happen
         enddo
         close(11)
      enddo
      close(31)

      end
