C Output from Public domain Ratfor, version 1.0
C
C Calculates Hamilton and Overlap matrix
C for Atoms i and j, the results are stored in ham,over
C
      subroutine slkmatrices(i,j,xat,ham,over)
      implicit none
      include 'maxima.inc'
      integer i,j, izp(nndim), nbeweg, l1,l2,nu,nv,nw, ind1,ind2,ind3
      real*8 xat(3,*), nel
      real*8 ham(ldim,ldim),over(ldim,ldim),ham_tmp(ldim,ldim),over_tmp(
     *ldim,ldim)
      real*8 boxsiz, xinvbox, xnullvec,skcut2
      logical period
      integer nlat
      common /box/ boxsiz(3,3),xinvbox(3,3),xnullvec(3), skcut2,nlat(3),
     *period
      integer nmax,nmin
      parameter (nmax = 10)
      parameter (nmin = 2)
      real*8 gmat,sumlat
      common /comgmat/gmat(4,2*nmax+1,2*nmax+1,2*nmax+1)
      common /comrmat/sumlat(3,2*nmax+1,2*nmax+1,2*nmax+1)
      common /izp/ izp, nel, nbeweg
      external skspar,skhpar
      integer izpi,izpj
      real*8 dif(3),dif0(3)

      dif0(1)=xat(1,j)-xat(1,i)
      dif0(2)=xat(2,j)-xat(2,i)
      dif0(3)=xat(3,j)-xat(3,i)
      izpi=izp(i)
      izpj=izp(j)

      do23000 l2 = 1,ldim 
      do23002 l1 = 1,ldim 
      ham(l1,l2) = 0.0d0
      over(l1,l2) = 0.0d0
23002 continue
23003 continue
23000 continue
23001 continue

      if(period)then
      do23006 nu =-nlat(1),nlat(1) 
      ind3=nu+nmax+1
      do23008 nv = -nlat(2),nlat(2) 
      ind2=nv+nmax+1
      do23010 nw = -nlat(3),nlat(3) 
      ind1=nw+nmax+1
      dif(1) = dif0(1)-sumlat(1,ind1,ind2,ind3)
      dif(2) = dif0(2)-sumlat(2,ind1,ind2,ind3)
      dif(3) = dif0(3)-sumlat(3,ind1,ind2,ind3)
      call slkode(dif,skcut2,izpi,izpj,ham_tmp,skhpar)
      call slkode(dif,skcut2,izpi,izpj,over_tmp,skspar)
      if(izpi .eq. -1969 )then
      izpi=izp(i)
      else
      do23014 l2 = 1,ldim 
      do23016 l1 = 1,ldim 
      ham(l1,l2) = ham(l1,l2) + ham_tmp(l1,l2)
      over(l1,l2) = over(l1,l2) + over_tmp(l1,l2)
23016 continue
23017 continue
23014 continue
23015 continue
      endif
23010 continue
23011 continue
23008 continue
23009 continue
23006 continue
23007 continue
      else
      call slkode(dif0,skcut2,izpi,izpj,ham,skhpar)
      call slkode(dif0,skcut2,izpi,izpj,over,skspar)
      endif
      end
      subroutine slkode(dum,skcut2,i,j,em,iovpar)
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
      real*8 dum(3),em(ldim,ldim)
      real*8 x(6),x2(6),dummy(ldim,ldim)
      external iovpar
      common /lmax/ lmax(maxtyp)
      r2=0.0
      do23018 l=1,3 
      x(l)=dum(l)
      x2(l)=x(l)*x(l)
      r2=r2+x2(l)
C!!!!!!!!!!!!!!!!!!!!!!
23018 continue
23019 continue
      if(r2.ge.skcut2)then
      i=-1969
      return
C!!!!!!!!!!!!!!!!!!!!!
      endif
      if(r2 .ge. 1.0e-8)then
      r2i=1.0/r2
      ri=sqrt(r2i)
      do23024 l=1,3 
      x(l)=x(l)*ri
      x(l+3)=x(l)
      x2(l)=x2(l)*r2i
      x2(l+3)=x2(l)
23024 continue
23025 continue
      maxmax = max(lmax(i),lmax(j))
      minmax = min(lmax(i),lmax(j))
C
C s Interaction
C
      call skss(x,x2,i,j,r2,iovpar,em,ldim)
      if(maxmax .le. 1)then
      return
C
C p Interaction
C
      endif
      if(minmax .ge. 2)then
      call skpp(x,x2,i,j,r2,iovpar,em(2,2),ldim)
      call sksp(x,x2,i,j,r2,iovpar,em(1,2),em(2,1),ldim)
      if(i .ne. j)then
      call sksp(x,x2,j,i,r2,iovpar,dummy,em(2,1),ldim)
      endif
      else
      if(lmax(j) .ge. 2)then
      call sksp(x,x2,i,j,r2,iovpar,em(1,2),em(2,1),ldim)
      else
      call sksp(x,x2,j,i,r2,iovpar,dummy,em(2,1),ldim)
      endif
      endif
      if(maxmax .le. 2)then
      return
C
C d Interaction
C
      endif
      if(minmax .eq. 3)then
      call skdd(x,x2,i,j,r2,iovpar,em(5,5),ldim)
      call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),ldim)
      call skpd(x,x2,i,j,r2,iovpar,em(2,5),em(5,2),ldim)
      if(i .ne. j)then
      call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),ldim)
      call skpd(x,x2,j,i,r2,iovpar,dummy,em(5,2),ldim)
      endif
      else
      if(lmax(i) .eq. 1)then
      call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),ldim)
      else
      if(lmax(i) .eq. 2)then
      call sksd(x,x2,i,j,r2,iovpar,em(1,5),em(5,1),ldim)
      call skpd(x,x2,i,j,r2,iovpar,em(2,5),em(5,2),ldim)
      else
      if(lmax(j) .eq. 1)then
      call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),ldim)
      else
      call sksd(x,x2,j,i,r2,iovpar,dummy,em(5,1),ldim)
      call skpd(x,x2,j,i,r2,iovpar,dummy,em(5,2),ldim)
      endif
      endif
      endif
      endif
      else
      do23046 k = 1,ldim 
      do23048 l = 1,ldim 
      em(k,l) = 0.
23048 continue
23049 continue
23046 continue
23047 continue
      if(i.ne.j)then
      return
      endif
      call selfs(i,j,r2,iovpar,em,ldim)
      if(lmax(i) .le. 1)then
      return
      endif
      call selfp(i,j,r2,iovpar,em(2,2),ldim)
      if(lmax(i) .le. 2)then
      return
      endif
      call selfd(i,j,r2,iovpar,em(5,5),ldim)
      endif
      end
