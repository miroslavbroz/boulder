c pop_decay.f

      subroutine pop_decay(Nanuli,nbinneg,nbins,npart_init,npart_end,
     :  npop,mpop,marr)

      include 'ucrm3.4.inc'
      integer Nanuli,nbinneg(Manuli),nbins(Manuli)
      real*8 npart_init(Manuli,BINNEG:BINMAX)
      real*8 npart_end(Manuli,BINNEG:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 mpop(Manuli,0:BINMAX)
      real*8 marr(Manuli,BINNEG:BINMAX)

      integer j,jj
      real*8 fdec, n, n_integer, n_remainder, chance
      real*8 ran3

c ... now force decay of the population by (npart_init-npart_end)/npart_init
      do j=1,Nanuli
        do jj=-nbinneg(j),nbins(j)

c .... compute fractional decay (fdec): different for each jj-bin
c      because npart_init may be different
          if (npart_init(j,jj) .gt. 0.d0) then
            fdec=(npart_init(j,jj)-npart_end(j,jj))/npart_init(j,jj) 
          else
            fdec=0.d0
          endif

c ... compute the number of bodies to be dynamically decayed
c ... and handle fractional probabity correctly (in order to
c ... account for really slooow decay rates)
          n = npop(j,jj,1)*fdec

          n_integer = dint(n)
          n_remainder = n - n_integer
          chance = ran3(idum0)
          if (chance .le. n_remainder) then
            n = n_integer + 1.d0
          else
            n = n_integer
          endif

          npop(j,jj,1) = dint(npop(j,jj,1) - n + 0.5d0)

c update npart_init by the applied fdec. The applied fdec may
c be different from original fdec because of the int operation
c WE MUST NOT (!) perform an operation like this if we already
c have a random-number generator working above...
c         ...

          npart_init(j,jj) = npart_init(j,jj)*(1.d0-fdec)

          if (jj.gt.0) mpop(j,jj)=marr(j,jj)*npop(j,jj,1) !done to preserve marr

        enddo  ! jj SFD bins
      enddo  ! j  (anuli)

      return
      end


