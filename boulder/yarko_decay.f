c yarko_decay.f
c Decay due to the Yarkovsky effect.
c Miroslav Broz (miroslav.broz@email.cz), Jun 20th 2013

c This decay similar to a radioactive one:
c
c   dN = -lambda N dt
c   N(t) = N0 exp(-lambda t)
c   lambda = 1/tau
c   tau = Delta_a / (da/dt)
c   N/N0 = exp(-t/tau)

      subroutine yarko_decay(Nanuli,yarko_n,yarko_r,yarko_tau,
     :  trans_m,dt,nbins,marr,sarr,npop,mpop)

      include 'ucrm3.4.inc'
      include 'yarko.inc'
      include 'trans.inc'

      integer Nanuli
      real*8 dt
      integer nbins(Manuli)
      real*8 marr(Manuli,BINNEG:BINMAX)
      real*8 sarr(Manuli,BINNEG:BINMAX)
      real*8 npop(Manuli,BINNEG:BINMAX,Ndata)
      real*8 mpop(Manuli,0:BINMAX)

      integer j,jj,k,nint
      real*8 tau, fdec, n, n_integer, n_remainder, chance
      logical extra

c functions
      real*8 linterp, ran3

      do j = 1, Nanuli
        do jj = BINNEG, BINMAX

          if (npop(j,jj,1).gt.0.d0) then

            tau = linterp(yarko_r(1,j),yarko_tau(1,j),yarko_n(j),
     :        sarr(j,jj),extra)
            if (extra) then
              write(*,*) '# yarko_decay: extrapolation is not allowed!'
              stop
            endif

            fdec = 1.d0-exp(-dt/tau)

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

            if (jj.gt.0) mpop(j,jj) = marr(j,jj)*npop(j,jj,1) ! to preserve marr

c ... transport bodies to other populations according to the matrix
c ... increase the number of bins if necessary!

            do k = 1,Nanuli
              if (trans_m(j,k).ne.0.d0) then
                npop(k,jj,1) = dint(npop(k,jj,1) + n*trans_m(j,k)+0.5d0)
                if (jj.gt.0) mpop(k,jj) = marr(k,jj)*npop(k,jj,1)
                if (jj.gt.nbins(k)) then
                  nbins(k) = jj
                  marr(k,jj) = marr(j,jj)
                  sarr(k,jj) = sarr(j,jj)
                endif
              endif
            enddo

          endif   ! npop.gt.0

        enddo     ! jj (bins)
      enddo       ! j  (anuli)

      return
      end


