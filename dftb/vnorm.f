C Output from Public domain Ratfor, version 1.0
C     SUBROUTINE VNORM
C     =================
C
C     Copyright 1991 by Peter Blaudeck
C *********************************************************************
C
C     PROGRAM CHARACTERISTICS
C     -----------------------
C
C Scaling of velocities, total impulse = 0
C
C PARAMETERS:
C nn        i number of atoms 
C x         r coordinates, a.u.
C xold      r      "      deltat earlier
C t         r temperature, K
C deltat    r time stepwidth, a.u.
C xm        r atom masses, nucleon masses
C izp       i sort numbers
C
C *********************************************************************
C
      subroutine vnorm(nn,x,xold,t,deltat,xm,izp)
      implicit real*8 (a-h,o-z)
      integer izp(*)
      real*8 x(3,*), xold(3,*), xm(*)
C 3/2*k in H/K
      parameter (r032 = 1.5 * 3.16679e-6)
      gm = 0.0
C
C total mass
C
      do23000 j = 1,nn 
      gm = gm +xm(izp(j))
23000 continue
23001 continue
      rcdt = 1.0/deltat
      rcgm = 1.0/gm
C
C set total impuls = 0
C
      do23002 i = 1,3 
      puls = 0.0
      do23004 j = 1,nn 
      v = (x(i,j) - xold(i,j))*rcdt
      puls = puls+xm(izp(j))*v
23004 continue
23005 continue
      puls = puls*deltat*rcgm
      do23006 j = 1,nn 
      xold(i,j) = xold(i,j) + puls
23006 continue
23007 continue
C
C rescale velocities according to temperature
C
23002 continue
23003 continue
      ekin = 0.0
      do23008 i = 1,nn 
      do23010 j = 1,3 
      v = (x(j,i) - xold(j,i))*rcdt
      ekin = ekin+0.5*xm(izp(i))*v**2
23010 continue
23011 continue
23008 continue
23009 continue
      if(ekin .gt. 0.0)then
      fak = sqrt(nn*r032*t/(1822.887428*ekin))
      else
      fak = 0.0
      endif
      do23014 i = 1,3 
      do23016 j = 1,nn 
      xold(i,j) = (xold(i,j) - x(i,j))*fak + x(i,j)
23016 continue
23017 continue
23014 continue
23015 continue
      end
