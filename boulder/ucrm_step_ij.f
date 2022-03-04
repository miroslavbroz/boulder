
c*************************************************************************
c                          UCRM_STEP_IJ                              
c*************************************************************************
c Main subroutine. Does one step of the algorithm by interacting anulus i
c with anulus j, over all their respective mass bins.
c Returns the changes in population mass, number, ecc and inc.
c Moreover it returns the bin number for isolated bodies and the suggested
c timestep so that the individual collision probability of the bodies
c is not too high

      subroutine ucrm_step_ij(i1st,time,dtcol,i,j,nbinneg,nbins,
     &     marr,sarr,mass,npop,a,e,inc,delta_a,
     &     mass_change,npop_change,fe,fi,dt_est,l_iso)

      include 'ucrm3.4.inc'
      
c Inputs
      real*8 marr(Manuli,BINNEG:BINMAX),sarr(Manuli,BINNEG:BINMAX)
      real*8 mass(Manuli,0:BINMAX),npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 time
      integer i,j,nbins(Manuli),nbinneg(Manuli)
      real*8 a(Manuli),e(Manuli,BINNEG:BINMAX),inc(Manuli,BINNEG:BINMAX)
      real*8 delta_a(Manuli)
      real*8 dtcol              ! collisional timestep 
      integer i1st
c Outputs
      real*8 npop_change(2,0:BINMAX,Ndata),mass_change(2,0:BINMAX)
      real*8 fe(2,0:BINMAX),fi(2,0:BINMAX),dt_est
      integer l_iso

c Internals
      real*8 macc,mf,mbin,mvel,h,h0,v,v0,HH,space,rhill,vhill,beta
      integer jj,k,l,lmax,new_k,kk
      real*8 npop_remove(2,BINNEG:BINMAX,Ndata) ! npop loss to erosion
      real*8 mass_remove(2,BINNEG:BINMAX) ! mass loss from bins due to collisions
      real*8 fe_rem(2,BINNEG:BINMAX),fi_rem(2,BINNEG:BINMAX) ! mass weighted eccentricity and inclination loss from bins due to collisions
      real*8 mass_debris(0:BINMAX)   ! mass gained from debris creation and accretion
      real*8 fe_deb(0:BINMAX),fi_deb(0:BINMAX)   ! eccentricity and inclination of the debris creation and accretion
      real*8 num_debris(0:BINMAX)   ! number of objects gained from debris creation and accretion
      real*8 fe_s(0:BINMAX),fi_s(0:BINMAX) ! debris eccentricity and incl.     
      real*8 deb_spect_m(0:BINMAX),deb_spect_n(0:BINMAX) ! debris mass and num.
      real*8 mlr,mej,mlf,mtot 
      real*8 ejr,diff,min_diff,kk_colls,rfg
      real*8 ncoll,ncoll_kl,vcoll,fg,ncoll_kl_ind,kest,lest
      real*8 mratio,qf,ncollmin,pcollmin,vcollmin,fgmin
      real*8 totmass,totm,totmM
      real*8 dm,dn,de,di
      real*8 Q_Qstar

c .......................................................................................
c...  Zero the arrays that accumulate debris exchange 
c...  Index 1 stands for target, index 2 stands for impactor 
      do jj = 0,nbins(j)
       mass_change(1,jj)=0.d0
       npop_change(1,jj,1)=0.d0
       fe(1,jj)=0.d0
       fi(1,jj)=0.d0
       mass_debris(jj) = 0.d0
       num_debris(jj) = 0.d0
       fe_deb(jj)=0.d0
       fi_deb(jj)=0.d0
       mass_remove(1,jj) = 0.d0 
       npop_remove(1,jj,1) = 0.d0
       fe_rem(1,jj)=0.d0
       fi_rem(1,jj)=0.d0
      enddo
      do jj = -nbinneg(i),nbins(i)
       mass_remove(2,jj) = 0.d0 
       npop_remove(2,jj,1) = 0.d0
       fe_rem(2,jj)=0.d0
       fi_rem(2,jj)=0.d0
       if (jj.ge.0) then
        mass_change(2,jj)=0.d0
        npop_change(2,jj,1)=0.d0
        fe(2,jj)=0.d0
        fi(2,jj)=0.d0
       endif
      enddo         
      npop_change(1,0,1)=0.d0
      npop_change(2,0,1)=0.d0
      fe(1,0)=0.d0
      fi(1,0)=0.d0
      fe(2,0)=0.d0
      fi(2,0)=0.d0
      dt_est=1.d10
      
c...  Loop over target bins
      do k = 1,nbins(j) 
       if(npop(j,k,1).lt.one) goto 200 ! label 200 = next_k,  
                                ! no collisions allowed on a target until 
                                ! there is a "real" object (i.e., npop>1)
c ...  Determine the largest impactor marr that is smaller than target marr  
       lmax = 0
       do l = 1,nbins(i) 
        if(marr(i,l).le.marr(j,k)) then  
         lmax = l 
        endif
       enddo
c isolated bodies in the same anulus cannot collide with each other
c        if(i.eq.j)then
c           call iso_bodies(i,nbins,marr,npop,a,delta_a,e,l_iso)
c           lmax=min(lmax,l_iso-1)
c        endif

c...  Loop over impactor bins 
       do l = -nbinneg(i),lmax 
        if(npop(i,l,1).lt.one) goto 100 ! label 100 = next_l, 
                                ! no collisions allowed unless there is 
                                ! a real projectile
        if((i.ne.j).or.(k.ne.l).or.(npop(i,l,1).ge.two)) then ! if target and projectile bins are the same there must be at least two bodies

c...  Calculate the collisional rate between i and j pops
c     Output variables are ncoll = the number of collisions of i on j per year 
c                      and vcoll = mean (weighted) collision velocity in cgs   
c     assuming one body in each population, PI cm^2 cross-section and
c     neglecting effects of gravitational focussing
c HERE IS WHERE THE UNITS BECOME IMPORTANT. NCOLL IS THE COLLISION PROBABILITY 
C PER CM2 PER YEAR. VCOLL IS IN CM/S. THIS MEANS THAT TIME MUST BE IN YR AND
C SIZE IN CM AND Q* IN ERG/S. semi major axis must be in AU
         call coll_rate(i,j,ncoll,vcoll)
         
c...  If the two tori do not intersect, return; patch! <- commented by MB
c         if (k.ge.l_iso) then
c          if (e(i,l).lt.(marr(j,k)/3.d0/msun)**0.33333333d0) ncoll=0.d0
c         endif
c                                               end patch
         if (ncoll.eq.0.d0) goto 100

c  for i.ne.j we divide by 2 because this k-l collision will be treated 
c  twice [because of if(marr(i,l).GE.marr(j,k)) condition]. 
c  For i=j we divide by 2 because otherwise each collision is counted 
c  twice [N(N-1)/2... see also KL98]               
c         if (marr(i,l).eq.marr(j,k)) then	! comparison of real*8 number is dangerous! MB
         if ((i.eq.j).and.(l.eq.k)) then
           ncoll=0.5d0*ncoll
         endif
c...  Calculate fg, the focusing factor for i,l on j,k collisions 
         call focusing_factor(marr(j,k),sarr(j,k),marr(i,l),
     &                        sarr(i,l),a(j),vcoll,fg)

c...  Calculate the number of collisions by one l on one k per yr. sarr must be in cm.             
         ncoll_kl = ncoll*fg*((sarr(j,k)+sarr(i,l))**2)
c...  The total number of collisions by l's on k's per yr    
         ncoll_kl_ind = ncoll_kl*npop(i,l,1) ! per target body
         ncoll_kl = ncoll_kl_ind*npop(j,k,1) ! total number of colls.
c...  Returns ncoll_kl=ncoll_kl*dtcol as an integer value if < 1d4, 
c...  rounded up or down based on a normally distributed random number
c...  based and the size of the remainder (don't see where < 1d4 condition is;
c...  seems like rounding always)
         call ucrm_round(ncoll_kl,dtcol,i1st)
         
         if (ncoll_kl.eq.0.d0) goto 100

c...  Calculate the debris mass (mej) and the mass fraction in the 
c...  largest fragment (mlf) and the slope of the fragment SFD
c...  due to one collision of k on l. MLR=mtot-mej is > marr(j) for accretional
c...  collisions.
         mtot = marr(j,k) + marr(i,l)
         call ucrm_ejecta_mass(j,i,marr(j,k),marr(i,l),vcoll,
     &                         mlr,mej,mlf,qf,Q_Qstar)

c...  Calculate the debris spectrum (mass, number) for a single 
c     l on k collision (i.e., per collision)
c...  debris spectrum INCLUDES the largest remnant  
         call ucrm_partition(nbins(j),j,k,e(j,k),inc(j,k),marr,
     &              mtot,mej,mlf,qf,deb_spect_m,deb_spect_n,fe_s,fi_s)

c output the families (to a file or to an array)
         if (write_rather_save) then
           call write_families(j,time,sarr(j,k),sarr(i,l),ncoll_kl,
     :       vcoll,mtot,mlr,mej,mlf,qf,Q_Qstar,rho(j),nbins(j),
     :       deb_spect_m,deb_spect_n,marr)
         else
           call save_families(j,time,sarr(j,k),sarr(i,l),ncoll_kl,
     :       vcoll,mtot,mlr,mej,mlf,qf,rho(j))
         endif

c output the craters
c j ... target pop.
c i ... projectile
c k ... bin
c l ... another bin
         call write_craters(j,time,sarr(j,k),sarr(i,l),ncoll_kl,
     :       rho(j),rho(i),vcoll,npop(j,k,1))

         totmass=0.d0
         do jj=0,nbins(j)
          totmass=totmass+deb_spect_m(jj)
         enddo
c  bunch of security checks.
         if(abs(totmass-mtot)/mtot.gt.1.d-14)then
          write(*,*)'Warning: poor mass conserv.'
          write(*,*)l,k,(totmass-mtot)/mtot
          write(*,*)marr(i,l),marr(j,k)
          write(*,*)mej,mtot-mej,sqrt(marr(j,0)*marr(j,1)),mlf,qf
          totmass=totmass-deb_spect_m(0)
          if(totmass.gt.mtot)then
           write(*,*)'Fatal: Creating mass'
           write(*,*)l,k,totmass,mtot
           stop
          endif
         endif
         totmass=totmass-deb_spect_m(0)
         if(l.le.0.and.totmass.gt.marr(j,k)) goto 100 ! avoid introduction of mass from the extrapolation tail.
c collisions with a projectile at the end of the tail do remove objects from the target bin. Thus the tail maight be too short
         if((deb_spect_n(k).eq.0.d0).and.(l.eq.-nbinneg(i))
     %      .and.(mtot-mej.lt.marr(j,k)))then
           if (dbg) then
             write(*,*)'# WARNING!: TAIL MIGHT BE TOO SHORT ',
     :         i,l,j,k,sarr(j,k)
           endif
         endif
c we want that the probability that a target is removed from its bin in one timestep is less than pthr. this sets a constraint on the timestep.
         if(deb_spect_n(k).eq.0.d0.and.l.gt.0)then
          if(dt_est.gt.pthr/ncoll_kl_ind)then
           dt_est=pthr/ncoll_kl_ind
          endif
         endif
c...  Sum the net per mass changes owing to debris generation 
c...  mass_debris is a variable that accumulates debris_spectrum into bins 
c...  over all combinations of l's and k's. Compute also a 
c     mass weighted mean e and i of debris.  
         totm=0.d0
         totmM=0.d0
c    target's j-bins
         do jj = 0,nbins(j)   
          if (jj.eq.k) then
           dm=deb_spect_m(jj)-marr(j,k)
           dn=deb_spect_n(jj)-1.d0
           de=fe_s(jj)*deb_spect_m(jj)-marr(j,k)*e(j,k)
           di=fi_s(jj)*deb_spect_m(jj)-marr(j,k)*inc(j,k)
           mass_change(1,jj)=mass_change(1,jj)+ncoll_kl*dm
           npop_change(1,jj,1)=npop_change(1,jj,1)+ncoll_kl*dn
           fe(1,jj)=fe(1,jj)+ncoll_kl*de
           fi(1,jj)=fi(1,jj)+ncoll_kl*di
          else
           mass_change(1,jj)=mass_change(1,jj)+
     %                       ncoll_kl*deb_spect_m(jj)
           npop_change(1,jj,1)=npop_change(1,jj,1)+
     %                         ncoll_kl*deb_spect_n(jj)
           fe(1,jj)=fe(1,jj)+ncoll_kl*fe_s(jj)*deb_spect_m(jj)
           fi(1,jj)=fi(1,jj)+ncoll_kl*fi_s(jj)*deb_spect_m(jj)
          endif
         enddo
c    projectile l-bin
         if(l.gt.0)then
          npop_change(2,l,1) = npop_change(2,l,1)-ncoll_kl
          mass_change(2,l) = mass_change(2,l)-ncoll_kl*marr(i,l)
c    mass weighted e and i of particles removed from projectile and target bins
          fe(2,l)=fe(2,l)-ncoll_kl*marr(i,l)*e(i,l)
          fi(2,l)=fi(2,l)-ncoll_kl*marr(i,l)*inc(i,l)
         else
          mass_change(2,0) = mass_change(2,0)-ncoll_kl*marr(i,l)
         endif
        endif  ! collision computation needed
 100    continue
       enddo                 ! (l)
 200   continue
      enddo                  ! (k)

      return
      end                    ! ucrm_step_ij
c----------------------------------------------------------------------------
