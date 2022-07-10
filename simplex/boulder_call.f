c boulder_call.f
c Generate initial SFD's and call boulder code to get evolved SFD's.
c Miroslav Broz (miroslav.broz@email.cz), Mar 4th 2011

      subroutine boulder_call(Y, SFD_n, SFD_time, SFD_nbins, SFD_diam,
     :  SFD_data)

      include '../boulder/ucrm3.4.inc'   ! all parameters file
      include '../boulder/yarko.inc'
      include '../boulder/trans.inc'
      include 'simplex.inc'
      include 'dependent.inc'
      real*8 Y(ndim)
      include 'cb_sfd.inc'
      include 'cb_fam.inc'
      
c mass array, size array
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX) 
c total mass array
      real*8 mpop(Manuli,0:BINMAX)
c number of particles array
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)! the Ndata vector is for 
c additional data. Ndata=1 is the populaiton in each bin. the only
c currently used orbital parameters arrays
      real*8 ecc(Manuli,BINNEG:BINMAX),inc(Manuli,BINNEG:BINMAX)
      real*8 axe(Manuli),delta_a(Manuli)
c number of anuli, mass bins, positive and negative
      integer Nanuli,nbins(Manuli),nbinneg(Manuli)
c bin number for isolated objects (WS93)
      integer l_iso(Manuli),l_iso_tmp
c changes in npop and mpop
      real*8 npop_change(2,0:BINMAX,Ndata)
      real*8 mpop_change(2,0:BINMAX)
      real*8 npop_change_tot(Manuli,0:BINMAX,Ndata)
      real*8 mpop_change_tot(Manuli,0:BINMAX)
c changes in eccentricity and inclination
      real*8 fe_tot(Manuli,0:BINMAX),fi_tot(Manuli,0:BINMAX)
      real*8 fe(2,0:BINMAX),fi(2,0:BINMAX)
c timesteps for the different parts of the algorithm
      real*8 dtcol,dt_est,dtstern,dtstern_tmp,dt_up

c auxiliary stuff
      integer i,j,jj,i1st,iwr,ifail,npart_n
      real*8 total_mass_new,total_mass_old,fme
      real*8 tend,t0,dtout,time,tout,dtmax
      real*8 te
      real*8 npart_init(Manuli,BINNEG:BINMAX)
      real*8 npart_end(Manuli,BINNEG:BINMAX)
      integer i1stcol,ierr
      real*8 xi
      real*8 ran3

c mira's additions
      real*8 tpart(npart_max),npart(npart_max,Manuli)   ! population decay data
      real*8 Pint_time(Pint_max)                        ! time-dependent Pint and vrel
      real*8 Pint_arr(Pint_max,Manuli,Manuli)
      real*8 vrel_arr(Pint_max,Manuli,Manuli)
      integer Pint_n

      real*8 rho_bulk_1, d1_1, d2_1, qa_1, qb_1, qc_1,
     :  dmax_1, d_norm_1, n_norm_1
      real*8 rho_bulk_2, d1_2, d2_2, qa_2, qb_2, qc_2,
     :  dmax_2, d_norm_2, n_norm_2
      integer idum,i2nd,iu
      real*8 q_fact_1, q_fact_2
      real*8 q_fact_DUMMY(Manuli), tend_MAX, dtcol_INI

      character*80 outname

      data i2nd /0/
      save i2nd
      save dtcol_INI,t0,tend_MAX,dtout
      save Pint_n,Pint_time,Pint_arr,vrel_arr
      save npart_n,tpart,npart

c ...................................................................

c Get parameters from Y() array
      rho_bulk_1 = Y(1)
      d1_1       = Y(2)
      d2_1       = Y(3)
      qa_1       = Y(4)
      qb_1       = Y(5)
      qc_1       = Y(6)
      dmax_1     = Y(7)
      n_norm_1   = Y(8)

      rho_bulk_2 = Y(9)
      d1_2       = Y(10)
      d2_2       = Y(11)
      qa_2       = Y(12)
      qb_2       = Y(13)
      qc_2       = Y(14)
      dmax_2     = Y(15)
      n_norm_2   = Y(16)

      q_fact_1   = Y(17)
      q_fact_2   = Y(18)
      idum       = int(Y(19))
      tend       = (Y(20) + dt_chi2)*1.d6

c some parameters have to be non-negative
      if (idum.gt.0) then
        idum = -idum - 1
      endif
      if (tend.lt.0.d0) then
        tend = -tend
      endif
      if (q_fact_1.lt.0.d0) then
        q_fact_1 = -q_fact_1
      endif
      if (q_fact_2.lt.0.d0) then
        q_fact_2 = -q_fact_2
      endif

c Generate initial SFD's
      idum0=idum
c      write(*,*) '# idum0 = ', idum0	! dbg
      xi=ran3(idum0)
c      write(*,*) '# xi = ', xi		! dbg
      i1st=0
      dbg=debug_boulder
      write_rather_save=.false.
      Nanuli = 2
      d_norm_1   = d1_1	! this cannot be free (due to condition d_norm .ge. d1)
      d_norm_2   = d1_2

      call gen_sfd_ic(rho_bulk_1, d1_1, d2_1, qa_1, qb_1, qc_1,
     :  dmax_1, d_norm_1, n_norm_1,
     :  1,axe,delta_a,nbins,marr,mpop,sarr,ecc,inc,npop)

      call gen_sfd_ic(rho_bulk_2, d1_2, d2_2, qa_2, qb_2, qc_2,
     :  dmax_2, d_norm_2, n_norm_2,
     :  2,axe,delta_a,nbins,marr,mpop,sarr,ecc,inc,npop)

c add Ceres and Vesta to the initial SFD's
      call add_ceres_vesta(n_ceres_vesta,ceres_vesta_diam,
     :  ceres_vesta_pop,nbins,marr,mpop,sarr,npop)

c output initial populations
      if (debug) then
        iu=10
        open(unit=iu, file="gen_sfd_ic.tmp", status="unknown")
        do j = 1,Nanuli
          do jj = 0,nbins(j) 
            write(iu,*) t0, j, marr(j,jj), sarr(j,jj), mpop(j,jj)
          enddo
        enddo
        close(iu)
      endif

c Read data files (only once)
      if (i2nd.eq.0) then
        call read_param(dtcol_INI,dtmax,t0,tend_MAX,dtout)

        call read_phys_par(Nanuli)

        if (Pint_tdepend) then
          call read_collprob_time(Nanuli,Pint_n,Pint_time,Pint_arr,
     :      vrel_arr)
        else
          call read_collprob(Nanuli)
        endif

        call read_pop_decay(Nanuli,npart_n,tpart,npart)
        call read_yarko(Nanuli,yarko_n,yarko_r,yarko_tau)
        call read_trans(Nanuli,trans_m)
        i2nd=1
      endif

      if (tend.gt.tend_MAX) then
        tend = tend_MAX
      endif
      dtcol = dtcol_INI
      q_fact(1) = q_fact_1
      q_fact(2) = q_fact_2

c now do the interpolation to find population at t0
      call pop_decay_interp(Nanuli,tpart,npart,npart_n,t0,npart_init)

      if (Pint_tdepend) then
        call coll_rate_interp(Nanuli,Pint_time,Pint_arr,vrel_arr,
     :    Pint_n,t0)
      endif

c check the input, generate empty arrays etc.
      do j = 1,Nanuli
       call ucrm_update_arrs(i1st,j,nbinneg(j),nbins(j),marr,sarr,
     &                       mpop,npop,axe,delta_a,ecc,inc)
      enddo
      total_mass_new = 0.d0
      do j = 1,Nanuli
       do jj = 0,nbins(j) 
        total_mass_new = total_mass_new + mpop(j,jj)
       enddo
      enddo

c Output of the initial conditions to an array
      SFD_n = 0
      FAM_n = 0
      call save_ob(Nanuli,t0,nbins,marr,sarr,mpop,
     :   SFD_n,SFD_time,SFD_nbins,SFD_diam,SFD_data)

      if (debug_boulder) then
        outname = "ob.data"
        call write_ob(outname,Nanuli,t0,nbins,marr,sarr,mpop)
      endif

ccc
c Big time loop <- THE SAME CODE AS IN BOULDER.F .................(begin)
ccc
      time=t0
      tout=t0+dtout

1000  continue                           ! timestep loop

        if (debug_boulder) then
          write(*,*) '# boulder_call(): t = ', time/1.d6, ' Myr, ',
     :      ' dt = ', dtcol/1.d6, ' Myr'
        endif

c timestep initializations:
      dtstern=1.d10
c zero the global arrays that accumulate debris exchange and e,i changes 
c and also the binaries arrays into fake ones...
      do j=1,Nanuli
       do jj = 0,nbins(j)
        mpop_change_tot(j,jj) = 0.d0
        npop_change_tot(j,jj,1) = 0.d0
        fe_tot(j,jj) = 0.d0
        fi_tot(j,jj) = 0.d0
       enddo
      enddo
            
c ... Loop over all pairs of tracers: i is the impactor population
c     j is the target population ..........................(begin)
      do i=1,Nanuli  
       do j=1,Nanuli
            
c ... Attempt one collisional step with timestep dtcol on the i,j pair  
c ... here we OUTPUT families too
        call ucrm_step_ij(i1st,time,dtcol,i,j,nbinneg,nbins,
     &                   marr,sarr,mpop,npop,axe,ecc,inc,delta_a,
     &                   mpop_change,npop_change,fe,fi,dt_est,l_iso_tmp)
        if (i.eq.j) l_iso(j)=l_iso_tmp ! save bin number for isolated objects
c ... first check. Internal subroutine step needs to be not smaller than dtcol
        if (dt_est.lt.dtcol) then
         dtcol=dt_est
         if (debug_boulder) then
           write(*,*)'# boulder_call: dt_est smaller than dt_coll: ',
     :       dt_est, dtcol
         endif
         goto 1000
        endif
c  update population numbers and masses. 
c  Returns  appropriate timestep
        call pop_change(i,j,nbins,mpop,npop,
     &     mpop_change,npop_change,mpop_change_tot,npop_change_tot,
     &     fe,fi,fe_tot,fi_tot,dtcol,dtstern,ifail)
         if (ifail.eq.1) goto 1000
       enddo                 ! over j
      enddo                  ! over i
c ... Loop over all pairs of tracers: i is the impactor population
c     j is the target population ............................(end)

c ... Check that no population or mass bin is negative. 
c     Return appropriate timestep
      call pop_check(Nanuli,nbins,mpop,npop,
     &     mpop_change_tot,npop_change_tot,dtcol,dtstern,ifail)
      if (ifail.eq.1) goto 1000

c ... Now we can update all variables; the step has been accepted
      call update_pops(Nanuli,nbins,mpop,npop,ecc,inc,
     &                 mpop_change_tot,npop_change_tot,fe_tot,fi_tot)

c ... Recompute marrs and sarrs from new mass and mpop (moving bins)
c ... this subroutine makes sure that npop is integer.
c ... Also checks for problems and re-arrange marrs 
      do j = 1,Nanuli
       call ucrm_update_arrs(i1st,j,nbinneg(j),nbins(j),marr,sarr,
     &                       mpop,npop,axe,delta_a,ecc,inc)
      enddo   

c ... Advance time after step has been accepted
      te=time+dtcol

ccc
c Now the dynamical decay is handled ............................(begin)
ccc
      call pop_decay_interp(Nanuli,tpart,npart,npart_n,te,npart_end)

      call pop_decay(Nanuli,npart_init,npart_end,npop,mpop,marr)
ccc
c Now the dynamical decay is handled ..............................(end)
ccc

c ... decay due to the Yarkovsky effect
      call yarko_decay(Nanuli,yarko_n,yarko_r,yarko_tau,
     :  trans_m,dtcol,nbins,marr,sarr,npop,mpop)

c ... Calculate new Pint and vrel if neccessary
      if (Pint_tdepend) then
        call coll_rate_interp(Nanuli,Pint_time,Pint_arr,vrel_arr,
     :    Pint_n,te)
      endif

c ... Calculate the new total mass   
      total_mass_new = 0.d0
      do j = 1,Nanuli
       do jj = 0,nbins(j)
        total_mass_new = total_mass_new + mpop(j,jj)
       enddo
      enddo
      time=te ! time which we are at

c Output data to an array
      if (time.ge.tout) then  ! here we OUTPUT

       call save_ob(Nanuli,time,nbins,marr,sarr,mpop,
     :   SFD_n,SFD_time,SFD_nbins,SFD_diam,SFD_data)

       if (debug_boulder) then
         call write_ob(outname,Nanuli,time,nbins,marr,sarr,mpop)
       endif

c       write(*,*) 'time = ', time/1.d6, ' Myr'	! dbg

c ... Set the new output time
       tout=time+dtout        
      endif

c stop
      if(time.ge.tend)then
        if (debug_boulder) then
          write(*,*) '# boulder_call(): t >= tend = ', time/1.e6, ' Myr'
        endif
        return
      else
c select the minimum of three stepsizes
       dtcol=min(dt_est/1.5d0,dtstern)
       if (dtcol.gt.dtmax) dtcol = dtmax ! added by Miroslav Broz, Jun 20th 2013
       if (time+dtcol.gt.tout) dtcol=tout-time
       goto 1000
      endif
ccc
c Big time loop ....................................................(end)
ccc
      return
      end


