ccccccccccccccccccccccccccccccc
      subroutine iso_bodies(i,nbins,marr,npop,a,delta_a,e,l_iso)

c compute bin number for isolated bodies
c isolated bodies: those for which the sum of mutual Hill radii*2sqrt(3)+2ae
c is less than the bin width (from WS93)
      include 'ucrm3.4.inc'

c     inputs
      integer i,nbins(Manuli)
      real*8 marr(Manuli,BINNEG:BINMAX),npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 a(Manuli),delta_a(Manuli)
      real*8 e(Manuli,BINNEG:BINMAX)
c     output
      integer l_iso
c     internals
      real*8 space
      integer kk,l

      l_iso=nbins(i)
      do l=1,nbins(i)
         space=0.d0
         if(npop(i,l,1).gt.0.)then
            do kk=l,nbins(i)
               space=space+npop(i,kk,1)*
     %              (2.*sqrt(3.)*a(i)*(2.*marr(i,kk)/3./msun)**(1./3.)
     %              +2.*a(i)*e(i,kk))
            enddo
            if(space.lt.delta_a(i))then
               l_iso=l
               return
            elseif(npop(i,l,1).eq.1..and.
     %          2.*sqrt(3.)*a(i)*(marr(i,l)/3./msun)**(1./3.)
     %              +2.*a(i)*e(i,l)
     %              .gt.delta_a(i))then
c we count as an isolated body also a single body whose Hill radius exceed the size of the anulus
               l_iso=l
               return
            endif
         endif
      enddo
      
      return
      end

