cc   code BOULDER
cc
cc   simulates the accretion and collisional fragmentation of planetesimals
cc   multi-anulus version
cc
cc   current version: Morby Jan 27, 2008
cc
cc   this version includes dynamics of KBO binaries
cc   version Jul 26, 2010 by David V
cc
cc   split to subroutines, time-dependent Pint, vrel
cc   by Miroslav Broz, Mar 2nd 2011
cc
cc   Yarkovsky effect (i.e. size-dependent decay)
cc   by Miroslav Broz, Jun 20th 2013
cc

      include '../boulder/ucrm3.4.inc'   ! all parameters file
      include '../boulder/yarko.inc'
      include '../boulder/trans.inc'
      
c mass array, size array
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX) 
c total mass array
      real*8 mpop(Manuli,0:BINMAX)
c number of particles array
c the Ndata vector is for additional data
c Ndata=1 is the populaiton in each bin
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)
c orbital parameters arrays
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
      real*8 dtmax

c auxiliary stuff
      integer i,j,jj,i1st,iwr,ifail,npart_n
      real*8 total_mass_new,total_mass_old,fme
      real*8 tend,t0,dtout,time,tout
      real*8 te
      real*8 npart_init(Manuli,BINNEG:BINMAX)
      real*8 npart_end(Manuli,BINNEG:BINMAX)
      integer i1stcol,ierr
      real*8 xi
      real*8 ran3

c bottke's changes
      integer ios
      character*80 boulder_name, outname

c mira's additions
      real*8 tpart(npart_max),npart(npart_max,Manuli)   ! population decay data
      real*8 Pint_time(Pint_max)                        ! time-dependent Pint and vrel
      real*8 Pint_arr(Pint_max,Manuli,Manuli)
      real*8 vrel_arr(Pint_max,Manuli,Manuli)
      integer Pint_n
      logical Pint_tdepend

c Read job parameters from 'boulder.in' file ...................(begin)

      open(1,file='boulder.in',status='old')
c (a) input SFDs filename
      read(1,*,iostat=ios) boulder_name
      if (ios.ne.0) then
        write(6,*) ' Error reading boulder_name'
        stop
      endif
      write(*,*) '# boulder_name = ', boulder_name
c (b) output SFDs filename
      read(1,*,iostat=ios) outname
      if (ios.ne.0) then
        write(6,*) ' Error reading outname'
        stop
      endif
      write(*,*) '# outname = ', outname
c (c) random seed
      read(1,*,iostat=ios) idum0
      if ((ios.ne.0).or.(idum0.ge.0)) then
        write (6,*) ' Error reading idum0 (negative)'
        stop
      endif
      write(*,*) '# idum0 = ', idum0
c (d) whether to use time-dependent Pint and vrel
      read(1,*,iostat=ios) Pint_tdepend
      if (ios.ne.0) then
        write(6,*) ' Error reading Pint_tdepend'
        stop
      endif
      write(*,*) '# Pint_tdepend = ', Pint_tdepend
c (e) whether to use time-dependent Pint and vrel
      read(1,*,iostat=ios) mfactor
      if (ios.ne.0) then
        write(6,*) ' Error reading mfactor'
        stop
      endif
      write(*,*) '# mfactor = ', mfactor
c (f) what disruptions to output explicitly
      read (1,*,iostat=ios) family_lfpb_min, family_lfpb_max,
     :  family_dpb_treshold, family_dlf_treshold
      if (ios.ne.0) then
        write (6,*) ' Error reading family_lfpb_min family_lfpb_max',
     :    ' family_dpb_treshold family_dlf_treshold'
        stop
      endif
      write(*,*) '# family_lfpb_min = ', family_lfpb_min
      write(*,*) '# family_lfpb_max = ', family_lfpb_max
      write(*,*) '# family_dpb_treshold = ', family_dpb_treshold, ' cm'
      write(*,*) '# family_dlf_treshold = ', family_dlf_treshold, ' cm'
c (g) what craters to output
      read (1,*,iostat=ios) crater_dpb_min, crater_dc_min
      if (ios.ne.0) then
        write(6,*) ' Error reading crater_dpb_min crater_dc_min'
        stop
      endif
      write(*,*) '# crater_dpb_min = ', crater_dpb_min, ' cm'
      write(*,*) '# crater_dc_max = ', crater_dc_min, ' cm'
      close(1)

c Read job parameters from 'boulder.in' file .....................(end)

c ... some initializations and checks of input/output files
      i1st=0
      dbg=.true.
      write_rather_save=.true.
      xi=ran3(idum0)
      write(*,*) "# xi = ", xi
      write(*,*) "# idum0 = ", idum0

      open(unit=1,file=outname,status='new',iostat = ios)
      if (ios.ne.0) then
        write (6,*) ' Error pre-existing output file!'
        stop
      endif
      close(1)

c Read data files ...............................................(begin)

      call read_ib(boulder_name, Nanuli,axe,delta_a,nbins,marr,mpop,
     :  sarr,ecc,inc,npop)

      call read_param(dtcol,dtmax,t0,tend,dtout)

      call read_phys_par(Nanuli,rho,Q0,a_benz,BB,b_benz,rho_bas,q_fact)

      if (Pint_tdepend) then
        call read_collprob_time(Nanuli,Pint_n,Pint_time,Pint_arr,
     :    vrel_arr)
      else
        call read_collprob(Nanuli)
      endif

      call read_pop_decay(Nanuli,npart_n,tpart,npart)

      call read_yarko(Nanuli,yarko_n,yarko_r,yarko_tau)

      call read_trans(Nanuli,trans_m)

c Read data files .................................................(end)

c Interpolation to find population at t0
      call pop_decay_interp(Nanuli,tpart,npart,npart_n,t0,npart_init)

      if (Pint_tdepend) then
        call coll_rate_interp(Nanuli,Pint_time,Pint_arr,vrel_arr,
     :    Pint_n,t0)
      endif

c Check the input, generate empty arrays, etc.
      do j = 1,Nanuli
        call ucrm_update_arrs(i1st,j,nbinneg(j),nbins(j),marr,sarr,
     :    mpop,npop,axe,delta_a,ecc,inc)
      enddo
      total_mass_new = 0.d0
      do j = 1,Nanuli
        do jj = 0,nbins(j) 
          total_mass_new = total_mass_new + mpop(j,jj)
        enddo
      enddo

c Output of the initial conditions
      call write_ob(outname,Nanuli,t0,nbins,marr,sarr,mpop)

c Big time loop .................................................(begin)

      time = t0
      tout = t0+dtout

1000  continue
      write(*,*) 'time  = ', time, ' y = ', time/1.d6, ' My'
      write(*,*) 'dtcol = ', dtcol, ' y = ', dtcol/1.d6, ' My'

c timestep initialization
      dtstern=1.d10

c Zero the global arrays that accumulate debris exchange and e,i changes 
c and also the binaries arrays into fake ones.
      do j = 1,Nanuli
        do jj = 0,nbins(j)
          mpop_change_tot(j,jj) = 0.d0
          npop_change_tot(j,jj,1) = 0.d0
          fe_tot(j,jj) = 0.d0
          fi_tot(j,jj) = 0.d0
        enddo
      enddo
            
c Loop over all pairs of tracers ................................(begin)

c i .. impactor population
c j .. target population
c k .. target bins
c l .. impactor bins

      do i = 1,Nanuli  
        do j = 1,Nanuli
            
c ... Attempt one collisional step with timestep dtcol on the i,j pair.
c ... Here we OUTPUT families too!
          call ucrm_step_ij(i1st,time,dtcol,i,j,nbinneg,nbins,
     :      marr,sarr,mpop,npop,axe,ecc,inc,delta_a,
     :      mpop_change,npop_change,fe,fi,dt_est,l_iso_tmp)

c ... Save bin number for isolated objects.
          if (i.eq.j) l_iso(j)=l_iso_tmp

c ... First check. Internal step needs to be NOT smaller than dtcol.
          if (dt_est.lt.dtcol) then
            dtcol=dt_est
            if (dbg) then
              write(*,*) '# dt_est = ', dt_est, ' y too small'
            endif
            goto 1000
          endif

c ... Update population numbers and masses. 
c ... Returns appropriate timestep.
          call pop_change(i,j,nbins,mpop,npop,
     :      mpop_change,npop_change,mpop_change_tot,npop_change_tot,
     :      fe,fi,fe_tot,fi_tot,dtcol,dtstern,ifail)

          if (ifail.eq.1) goto 1000

        enddo  ! over j
      enddo  ! over i

c Loop over all pairs of tracers ..................................(end)

c ... Check that no population or mass bin is negative. 
c ... Returns appropriate timestep.
      call pop_check(Nanuli,nbins,mpop,npop,
     :  mpop_change_tot,npop_change_tot,dtcol,dtstern,ifail)

      if (ifail.eq.1) goto 1000

c ... Now we can update all variables; the step has been accepted.
      call update_pops(Nanuli,nbins,mpop,npop,ecc,inc,
     :  mpop_change_tot,npop_change_tot,fe_tot,fi_tot)

c ... Recompute marrs and sarrs from new mass and mpop (moving bins).
c ... This subroutine makes sure that npop is integer.
c ... Also checks for problems and re-arrange marrs.
      do j = 1,Nanuli
        call ucrm_update_arrs(i1st,j,nbinneg(j),nbins(j),marr,sarr,
     :    mpop,npop,axe,delta_a,ecc,inc)
      enddo   

c ... Advance time after step has been accepted.
      te = time+dtcol

c ... Now dynamical decay is handled.
      call pop_decay_interp(Nanuli,tpart,npart,npart_n,te,npart_end)

      call pop_decay(Nanuli,nbinneg,nbins,npart_init,npart_end,
     :  npop,mpop,marr)

c ... Decay due to the Yarkovsky effect.
      call yarko_decay(Nanuli,yarko_n,yarko_r,yarko_tau,
     :  trans_m,dtcol,nbinneg,nbins,marr,sarr,npop,mpop)

c ... Calculate new Pint and vrel (if neccessary),
      if (Pint_tdepend) then
        call coll_rate_interp(Nanuli,Pint_time,Pint_arr,vrel_arr,
     :    Pint_n,te)
      endif

c ... Calculate the new total mass.
      total_mass_new = 0.d0
      do j = 1,Nanuli
        do jj = 0,nbins(j)
          total_mass_new = total_mass_new + mpop(j,jj)
        enddo
      enddo

      time = te

c ... Output data to file.
      if (time.ge.tout) then
        call write_ob(outname,Nanuli,time,nbins,marr,sarr,mpop)
        tout = time+dtout        
      endif

c ... Select the minimum of three stepsizes.
      if (time.lt.tend) then
        dtcol=min(dt_est/1.5d0,dtstern)
        if (dtcol.gt.dtmax) dtcol = dtmax
        if (time+dtcol.gt.tout) dtcol=tout-time
        write(*,*) '.'
        goto 1000
      endif

c Big time loop ...................................................(end)

c Done!
      write(*,*) 'time = ', time, ' y = ', time/1.d6, ' My (exiting)'
      close(11)
      stop

      end


