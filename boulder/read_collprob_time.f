c read_collprob_time.f
c Read time-dependent collisional probabilities and velocities.
c Miroslav Broz (miroslav.broz@email.cz), Mar 3rd 2011

      subroutine read_collprob_time(Nanuli,n,Pint_time,Pint_arr,
     :    vrel_arr)

      include 'ucrm3.4.inc'
      integer Nanuli,n
      real*8 Pint_time(Pint_max)
      real*8 Pint_arr(Pint_max,Manuli,Manuli)
      real*8 vrel_arr(Pint_max,Manuli,Manuli)

      integer i,j,k,ierr
      character*80 filename

      filename = 'collprob_time.dat'
      open(unit=1,file=filename,status='old')
      call skip(1)

      k=0
5     continue
        k=k+1
        if (k.gt.Pint_max) then
          write(*,*) '# read_collprob_time: index k = ', k,
     :      ' .gt. Pint_max = ', Pint_max
          stop
        endif

        read(1,*,iostat=ierr,end=15) Pint_time(k),
     :    ( (Pint_arr(k,i,j), j=1,Nanuli), i=1,Nanuli ),
     :    ( (vrel_arr(k,i,j), j=1,Nanuli), i=1,Nanuli )

        if (ierr.ne.0) then
          write(*,*) '# read_collprob_time: error reading ', filename
          stop
        endif
      goto 5

15    continue
      close(1)
      n=k-1

      do i=1,Nanuli
        do j=1,Nanuli
          do k=1,n
            Pint_arr(k,i,j)=Pint_arr(k,i,j)/(1.d5)**2       ! notice, we convert to cm^2 here (cm/s; cm^-2 y^-1)
            vrel_arr(k,i,j)=vrel_arr(k,i,j)*1.d5
          enddo
        enddo
      enddo

      write(*,*) '# read_collprob_time: n = ', n

      return
      end


