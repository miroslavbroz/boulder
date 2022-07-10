
      implicit none ! you got it baby	

c   parameters related to dimensions of vectors
      integer Manuli,BINMAX,BINNEG,Ndata
c   Manuli is the maximum number of anuli 
c   BINMAX is the maximum number of positive mass bins in each anulus
c   BINNEG is the maximum number of `negative' mass bins used in the extrapolation of the size distribution at the small-size end to avoid truncature effects
c   Ndata is the number of informations carried by the npop vector
      parameter(Manuli=4,BINMAX=200,BINNEG=-30,Ndata=1) 

c   some physical constants
      real*8 msun               ! G*M_sun in Swift units
      parameter (msun=1.989d33) ! in grams
      real*8 one,two,tiny
      parameter (one=0.8d0,two=1.8d0,tiny=.1d-8) 
      real*8 AU,YEAR
      parameter (AU=149598000.d5,YEAR=31557600.d0) ! in cm and sec.
      real*8 pi
      parameter (pi=3.141592653589793d0)
      real*8 G
      parameter (G=6.67300d-8) ! in cm3/g/s2	

c some physical parameters for the planetesimal population
      real*8 rho(Manuli)     ! material density in cgs units one for each annulus  
      real*8 Q0(Manuli),a_benz(Manuli),BB(Manuli),
     $	b_benz(Manuli),rho_bas(Manuli), q_fact (Manuli) ! parameters for Qstar (one per annulus)
c collision probability and impact velocity matrix
      real*8 Pint(Manuli,Manuli),vrel(Manuli,Manuli)

c  parameters related to the simulation scheme
      real*8 mfactor
c      parameter (mfactor=2.d0) ! mass ratio between mass bins	! this is now read from input
      real*8 pthr
      parameter (pthr=0.1d0) ! 10% tolerance on population change (for timestep calculation in subroutine _step) 
      integer Nmul
      parameter (Nmul=10) !multiple of marr(0) for which interpolation is used when a new bin is created in subroute update_arrs
      real*8 tollS
      parameter (tollS=.8d0) ! tollerance used in subroutines pop_change and pop_check. The suggested timestep is for the population to decrease no more than toll.

c  parameters related to the activity of the code
      integer ifrag
      parameter (ifrag=1) ! if 1/0 collision fragmentation is taken or not into account 

c stuff related to e,i evolution. Not used in collision only code

c some physical parameters concerning the environment
      real*8 diffvel
      parameter (diffvel=6000.d0) ! cm/s relative velocity between gas and particle
      real*8 fgas
      parameter (fgas=1.d0) ! equivalent to 1.4 MMSN
      real*8 gamma
      parameter (gamma=0.001d0) ! strength of turbulence
      real*8 rho_gas
      parameter (rho_gas=1.18d-9) ! density of gas at 1 AU (in g/cm3)

      real*8 toll
      parameter (toll=0.1d0) ! we allow only a 10% change in e or i. used to cmpute the stepsize in subroutine update_ei
      real*8 Vlv
      parameter (Vlv=1.d0) ! low velocity limit in Hill velocity (KL99 p. 1117) used in subroutine update_ei
      integer icd
      parameter (icd=1) ! if 1/0 collisional damping is taken or not into account 
      integer igd
      parameter (igd=1) !if 1/0 gas drag is taken or not into account 
      integer iturb
      parameter (iturb=1)  !if 1/0 turbulent stirring is taken or not into account 
      integer itide
      parameter (itide=0)  !if 1/0 tidal damping is taken or not into account 

c additions by Mira Broz, Mar 2nd 2011
      integer npart_max			! maximum number of data in pop_decay.dat
      parameter (npart_max=300000)
      integer Pint_max			! maximum number of data in collprob_time.dat
      parameter (Pint_max=300000)
      real*8 family_lfpb_min, family_lfpb_max,
     :   family_dpb_treshold, family_dlf_treshold	! which disruptions to output
      real*8 crater_dpb_min, crater_dc_min
      integer idum0			! iseed for random-number generator
      logical dbg			! debugging?
      logical write_rather_save		! write to file or save data to an array?
c      real*8 cnt	! dbg

c Yarkovsky effect (i.e. size-dependent decay)
      integer yarko_max
      parameter (yarko_max=1000)

c common blocks here
      common/phys_par/rho,Q0,a_benz,BB,b_benz,rho_bas,q_fact
      common/impacts/Pint,vrel
      common/debug/dbg
      common/seed/idum0
      common/cb_mfactor/mfactor
      common/cb_write/write_rather_save
      common/families/family_lfpb_min, family_lfpb_max,
     :  family_dpb_treshold, family_dlf_treshold
      common/craters/crater_dpb_min, crater_dc_min
c      common/cb_cnt/cnt	! dbg

