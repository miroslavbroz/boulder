
c*********************************************************************************
c                                UCRM_ROUND
c*********************************************************************************      
c     Round off any input vector < Input_min to a statistical value
c     Alan Stern		                          30 Mar 2003
c     c/o OK                                          30 Mar 2003
c     translated to Fortran, D. Nesvorny              16 Aug 2003
c ********************************************************************************

      subroutine ucrm_round(x,dt,i1st)

      include 'ucrm3.4.inc'
c     input
      real*8 dt
c     input/output
      real*8 x
      integer i1st
c     internals
      real*8 chance,xin,x0,xr
      real*8 ran3
c      integer idum0	! moved to ucrm3.4.inc!
c     save iseed
c     external ran

c... Calcs
      x = x*dt                                ! we multiply by the time now
      x0  = dint(x)                           ! Mantissa
      xr  = x - x0                            ! Remainder

      if (xr .ge. 1.d0) then
       write(*,*) 'Error in ucrm_round. Exiting...'
       stop           ! Error (could never happen, I guess)
      end if
         
c...  Uniform random number 0 to 1
      chance=ran3(idum0)
ccc      chance=grnd()

c... Add one to the floor if random is "thumbs up" 
      if (chance .le. xr) x = dint(x) + 1.d0 ! round up (if condition is true) 
      if (chance .gt. xr) x = dint(x)        ! else use floor

c... Done
      return
      end                       ! ucrm_round
c---------------------------------------------------------------------------------

