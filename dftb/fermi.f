C This code was modified by Yasuhito Ohta   8/1/2008
      subroutine fermi(nel,telec,ndim,ev,dacc,occ,efermi)
C
C maxima definitions
C
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
      integer ndim
      real*8 nel,telec,efermi,dacc,ev(mdim),occ(mdim)
C
C Boltzmann constant in H/K
C
      parameter(ckbol = 3.16679e-6)
C
C degeneracy tolerance for T = 0 
C
      parameter(degtol = 1.0e-4)
      real*8 beta,etol,occdg,chleft,ef1,ef2,ceps,eeps,charge,fac
      integer i,nef1,nef2,nup,ndown,nocc2,ndeg,istart,iend
      logical tzero
      common /machine/ racc

c   parallel  trajectories execution 
      include 'mpiio.inc'     !ninfile, noutfile
C              machine accuracy
C start defining occupation numbers and their derivatives 
C
      do i = 1,ndim 
        occ(i) = 0.0
      enddo

      if(nel .lt. 1.0e-5)then
        goto 10
      endif

      if(nel .gt. 2*ndim)then
        write(istdout,*)'too many electrons'
        stop
      endif
C
C find energy interval for Fermi energy
C
      if(telec .gt. 5.0)then
        beta = 1.0/(ckbol*telec)
        etol = ckbol*telec*(log(beta)-log(racc))
        tzero = .false.
      else
        etol = degtol
        tzero = .true.
      endif
C
C find energy range for Fermi energy
C
      if(nel.gt.int(nel))then
        nef1 = (nel+2)/2
        nef2 = (nel+2)/2
      else
        nef1 = (nel+1)/2
        nef2 = (nel+2)/2
      endif

      efermi = 0.5*(ev(nef1)+ev(nef2))
      nup = nef1
      ndown = nef1

23010 if(nup .lt. ndim)then
        if(abs(ev(nup+1)-efermi) .gt. etol)then
          goto 23011
        endif
        nup = nup+1
        goto 23010
      endif
23011 continue

23014 if(ndown .gt. 0)then
        if(abs(ev(ndown) - efermi) .gt. etol)then
          goto 23015
        endif
        ndown = ndown-1
        goto 23014
      endif
23015 continue

      ndeg = nup - ndown
      nocc2 = ndown

      do i = 1,nocc2 
        occ(i) = 2.0
      enddo

      if(ndeg .eq. 0)then
       goto 10
      endif
C
C for T == 0 occupy orbitals as usual
C

      if(tzero)then
        occdg = ndeg
        occdg = (nel-2*nocc2)/occdg
        do i = nocc2+1,nocc2+ndeg 
          occ(i) = occdg
        enddo

      else
C
C for finite T, use Fermi distribution
C
      chleft = nel-2*nocc2
      istart = nocc2+1
      iend = istart+ndeg-1

      if(ndeg .eq. 1)then
        occ(istart) = chleft
        goto 10
      endif

C
C bracket and iterate Fermi energy by bisection
C
      ef1 = efermi-etol-degtol
      ef2 = efermi+etol+degtol
      ceps = dacc*chleft
      eeps = dacc*max(abs(ef1),abs(ef2))
23028 continue

      efermi = 0.5*(ef1+ef2)
      charge = 0.0

      do i = istart,iend 
        occ(i) = 2.0/(1.0+exp(beta*(ev(i)-efermi)))
        charge = charge+occ(i)
      enddo

      if(charge .gt. chleft)then
        ef2 = efermi
      else
        ef1 = efermi
      endif

23029 if(.not.((abs(charge-chleft) .lt. ceps) .or. 
     &        (abs(ef1-ef2) .lt. eeps))) goto 23028

23030 continue

C rescale occ(i) if the accuracy was not good enough 
C
        if(abs(charge-chleft) .gt. ceps)then
          fac = chleft/charge
          do i = istart,iend 
            occ(i) = occ(i)*fac
          enddo
        endif
      endif

10    continue
      end
