C Output from Public domain Ratfor, version 1.0
C     SUBROUTINE XNORM
C     =================
C
C     Copyright 1991 by Peter Blaudeck
C *********************************************************************
C
C     PROGRAM CHARACTERISTICS
C     -----------------------
C
C XNORM moves the mass centre to the origin
C
C PARAMETERS:
C nn       i number of atoms
C xx       r coordinates
C xold     r      "      deltat earlier
C xm       r atom masses
C izp      i sort numbers
C
C *********************************************************************
C
      subroutine xnorm(nn,xx,xold,xm,izp)
      implicit real*8 (a-h,o-z)
      integer izp(*)
      real*8 xx(3,*), xold(3,*), xm(*)
      gmrc = 0.0
      do23000 i = 1,nn 
      gmrc = gmrc+xm(izp(i))
23000 continue
23001 continue
      gmrc = 1.0/gmrc
      do23002 i = 1,3 
      sum = 0.0
      do23004 j = 1,nn 
      sum = sum+xx(i,j)*xm(izp(j))
23004 continue
23005 continue
      sum = sum*gmrc
      do23006 j = 1,nn 
      xx(i,j) = xx(i,j)-sum
      xold(i,j) = xold(i,j)-sum
23006 continue
23007 continue
23002 continue
23003 continue
      end
