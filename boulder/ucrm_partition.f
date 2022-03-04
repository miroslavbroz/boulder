
c--------------------------------------------------------------------------
c  distributes the fragments according to the existent mass bins
      subroutine ucrm_partition(nbins,j,k,e1,inc1,marr,
     &              mtot,mej,Mlf,qf,
     %              deb_spect_m,deb_spect_n,fe_s,fi_s)  

      include 'ucrm3.4.inc'
c     Input
      integer nbins,j,k
      real*8 e1,inc1,marr(Manuli,BINNEG:BINMAX),mtot,mej,Mlf,qf
c     mej=ejected mass
c     mlf=mass of the largest fragment (not remnant: the beginning of the SFD)
c     qf=cumulative slope of SFD
c     e1,inc1=e and i of progenitor
c     Output
      real*8 deb_spect_m(0:BINMAX),deb_spect_n(0:BINMAX)
      real*8 fe_s(0:BINMAX),fi_s(0:BINMAX)
c     Internals
      real*8 mlr,minimum,ratio,r_lf,r
      integer mlloc,mllocMlr
      real*8 targ,val,der,dx,r_break,r_low,r_high
      real*8 ftp,onethird,Pref,Prefprim,totmass,totnum,DonMass
      real*8 totmasstmp,totnumtmp
      integer jj,kk

c initialization stuff....
      ftp=4.d0/3.d0*pi*rho(j) ! NOTE we assume all fragments have the same
                              ! physical properties (rho etc.) of the target
      onethird=1.d0/3.d0

      do jj=0,nbins     ! this was changed from 1,nbins to 0,nbins by MB
         deb_spect_m(jj)=0.d0
         deb_spect_n(jj)=0.d0
         fe_s(jj)=0.d0
         fi_s(jj)=0.d0
      enddo
      totmass=0.d0
      totnum=0.d0

      Mlr=Mtot-Mej ! mass of largest remnant
c...  Find the bin of largest debris object, i.e "marr(jj)" that is closest to "ml"

      if(Mlr.gt.sqrt(marr(j,0)*marr(j,1)))then
         minimum = marr(j,nbins)/marr(j,0)  
         mlloc = 0
         do jj = 1,nbins
            ratio=mlr/marr(j,jj)
            if(ratio.lt.1.d0)ratio=1.d0/ratio
            if(ratio.lt.minimum) then
               minimum = ratio
               mlloc = jj       ! mlloc is the index location of the largest debris object bin
            end if
         end do
         deb_spect_m(mlloc)=mlr
         deb_spect_n(mlloc)=1.d0
         fe_s(mlloc)=e1   ! largest remanant inherits e and i of progenitor
         fi_s(mlloc)=inc1
         totmass=Mlr
         mllocMLr=mlloc
      else
c         deb_spect_m(0)=Mlr
c         return
      endif
  

c     now build the SFD starting from the larges fragment DOWNWARDS
      r_LF=(Mlf/ftp)**onethird
      Pref=-qf/r_LF**qf ! normalizing factor in the diff distribution so that the integral from r_LF to infinity of the differential distribution is 1
      qf=qf-1.d0 ! now it is the differential slope
      

      if(Mlf.lt.sqrt(marr(j,0)*marr(j,1)))then
         deb_spect_m(0)=mtot-totmass ! even the largest fragment is in the dust bin
         return
      endif
c search for bin of the largest fragment
      minimum = marr(j,nbins)/marr(j,0)  
      mlloc = 0
      do jj = 1,nbins
         ratio=Mlf/marr(j,jj)
         if(ratio.lt.1.d0)ratio=1.d0/ratio
         if(ratio.lt.minimum) then
            minimum = ratio
            mlloc = jj     
         end if
      end do
c  compute mass and number of objects in that bin (adding Mlf and 1 because the integration starts from r_LF)
      r_low=(sqrt(marr(j,mlloc)*marr(j,mlloc-1))/ftp)**onethird
      if(qf.ne.-1.d0)then
         deb_spect_n(mlloc)=Pref*(r_LF**(1.d0+qf)
     %        -r_low**(1.d0+qf))/(1.d0+qf)+1.d0
         if(deb_spect_n(mlloc).lt.1.d4)then ! too big for an integer otherwise
            deb_spect_n(mlloc)=int(deb_spect_n(mlloc))
         endif
         if(deb_spect_n(mlloc).gt.1.d0)then ! compute the radius at which the integral is integer and equal to deb_spect_n(mlloc)-1.d0
            r_low=(r_LF**(1.d0+qf)-(deb_spect_n(mlloc)-1.d0)
     %           *(1.d0+qf)/Pref)**(1.d0/(1.d0+qf))
         else
            r_low=r_LF
         endif
      else ! same stuff, but with logarithmic integrand
         deb_spect_n(mlloc)=Pref*(log(r_LF)-log(r_low))+1.d0
         if(deb_spect_n(mlloc).lt.1.d4)then ! too big for an integer otherwise
            deb_spect_n(mlloc)=int(deb_spect_n(mlloc))
         endif
         if(deb_spect_n(mlloc).gt.1.d0)then
            r_low=exp(log(r_LF)-(deb_spect_n(mlloc)-1.d0)/Pref)
         else
            r_low=r_LF
         endif
      endif
c     now compute the mass. If there is only one body we add the real mass. Otherwise we use the integral of the SFD
      if(deb_spect_n(mlloc).eq.2.d0)then
         deb_spect_m(mlloc)=Mlf+ftp*r_low**3
      elseif(qf.ne.-4.d0)then
         deb_spect_m(mlloc)=Pref*ftp*(r_LF**(4.d0+qf)-r_low**(4.d0+qf))
     %        /(4.d0+qf)+Mlf
      else
         deb_spect_m(mlloc)=Pref*ftp*(log(r_LF)-log(r_low))+Mlf
      endif
      totmass=totmass+deb_spect_m(mlloc)
c      if(totmass.gt.mtot)then	! this test was problematic if the numbers are almost equal; modified by M. Broz (Apr 2nd 2012)
      if ((totmass/mtot-1.d0).gt.1.d-12) then
         write(*,*) 'problem: too much mass in Mlf'
         write(*,*) 'totmass/mtot-1 = ', totmass/mtot-1.d0
         write(*,*) 'mlloc = ', mlloc
         write(*,*) 'deb_spect_m-Mlf = ', deb_spect_m(mlloc)-Mlf
         write(*,*) 'Mlf     = ', Mlf
         write(*,*) 'Mlr     = ', Mlr
         write(*,*) 'Mlr+Mlf = ', Mlr+Mlf
         write(*,*) 'totmass = ', totmass
         write(*,*) 'mtot    = ', mtot
         stop
      endif
      totnum=deb_spect_n(mlloc)+1.d0
      if(mlloc.eq.mllocMlr)then
         deb_spect_m(mlloc)=deb_spect_m(mlloc)+Mlr
         deb_spect_n(mlloc)=deb_spect_n(mlloc)+1.d0
      endif
      fe_s(mlloc)=e1     ! all fragments inherit the e and i of the progenitor
      fi_s(mlloc)=inc1

c now do the remaining bins. Same algorithms as above for Mlf
      do jj=mlloc-1,1,-1
         r_high=r_low
         r_low=(sqrt(marr(j,jj)*marr(j,jj-1))/ftp)**onethird
 76      if(qf.ne.-1.d0)then
            deb_spect_n(jj)=Pref*(r_high**(1.d0+qf)-r_low**(1.d0+qf))
     %           /(1.d0+qf)
            if(deb_spect_n(jj).lt.1.d4)then ! too big for an integer otherwise
               deb_spect_n(jj)=int(deb_spect_n(jj))
            endif
            if(deb_spect_n(jj).ne.0.d0)then ! compute the radius at which the integral is integer and equal to deb_spect_n(jj)
               r_low=(r_high**(1.d0+qf)-deb_spect_n(jj)
     %              *(1.d0+qf)/Pref)**(1.d0/(1.d0+qf))
            else
               r_low=r_high
            endif
         else ! same as before but with logarithmic integrand
            deb_spect_n(jj)=Pref*(log(r_high)-log(r_low))
            if(deb_spect_n(jj).lt.1.d4)then ! too big for an integer otherwise
               deb_spect_n(jj)=int(deb_spect_n(jj))
            endif
            if(deb_spect_n(jj).ne.0.d0)then
               r_low=exp(log(r_high)-deb_spect_n(jj)/Pref)
            else
               r_low=r_high
            endif
         endif
c now add the mass
         if(deb_spect_n(jj).eq.0.d0)then
            deb_spect_m(jj)=0.d0
         else
            if(qf.ne.-4.d0)then
               deb_spect_m(jj)=Pref*ftp*(r_high**(4.d0+qf)
     %              -r_low**(4.d0+qf))/(4.d0+qf)
            else
               deb_spect_m(jj)=Pref*ftp*(log(r_high)-log(r_low))
            endif
c if the mean mass of the fragment is not in the bin add the masses one by one 
            if(deb_spect_m(jj)/deb_spect_n(jj).gt.
     %           sqrt(marr(j,jj)*marr(j,jj+1)).or.
     %         deb_spect_m(jj)/deb_spect_n(jj).lt.
     %           sqrt(marr(j,jj-1)*marr(j,jj)))then !count objects one by one
               deb_spect_m(jj)=0.d0
               do kk=1,int(deb_spect_n(jj))
c  compute the size of each object (for which the integrand is unity)
                  if(qf.ne.-1.d0)then
                     r=(r_high**(1.d0+qf)-kk
     %                    *(1.d0+qf)/Pref)**(1.d0/(1.d0+qf))
                  else
                     r=exp(log(r_high)-kk/Pref)
                  endif
                  deb_spect_m(jj)=deb_spect_m(jj)+ftp*r**3 ! add its mass
               enddo
            endif
         endif
c     total mass and number at this point of the SFD. It is a tmp variable because we might have integrated the qf slope too far. So, we need to check that if we had an SFD with a Donhany slope from this point down to 0, the total mass would be smaller than mej. Otherwise we need to take a step back and use the Donhany slope already in the last bin.
         totmasstmp=totmass+deb_spect_m(jj) 
         totnumtmp=totnum+deb_spect_n(jj)
c   normalization coefficient if the SFD were Donahny from this point downwards
         Prefprim=2.5d0*totnumtmp/r_low**(-2.5d0)
         DonMass=Prefprim*ftp*2.d0*sqrt(r_low) ! total mass of a Donhany distribution down to r=0
         if(totmasstmp+DonMass.gt.mtot.and.qf.ne.-3.5d0)then ! use Donhani from this last bin
            qf=-3.5d0
            Pref=2.5d0*totnum/r_high**(-2.5d0)
            goto 76
         elseif(totmasstmp+DonMass.gt.mtot.and.qf.eq.-3.5d0)then ! not enough mass even for Donhany. Skip bin
            r_low=(sqrt(marr(j,jj)*marr(j,jj-1))/ftp)**onethird
            Pref=2.5d0*totnum/r_low**(-2.5d0)
            deb_spect_m(jj)=0.d0
            deb_spect_n(jj)=0.d0
            goto 77
         endif
c step succesful. we don't need to change slope yet
         fe_s(jj)=e1
         fi_s(jj)=inc1      
c total mass and number at this point of the SFD
         totmass=totmass+deb_spect_m(jj)
         totnum=totnum+deb_spect_n(jj)
 77      continue
      enddo

c     now put in deb_spect_m(0) the rest
      deb_spect_m(0)=mtot-totmass
      if(deb_spect_m(0).lt.0.d0)then
         write(*,*)'problem, deb_spect_m(0)<0:',deb_spect_m(0)
         stop
      endif
      
      return
      end

