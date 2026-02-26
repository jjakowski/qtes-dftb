C 8/1/2008 This code was modified by Yasuhito Ohta
C**********************************************************************
      integer function skspar(i,j,r2,dd)
C**********************************************************************
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
      integer i,j
      real*8 r2,dd(13)
      integer dim(maxtyp,maxtyp),mxind
      real*8 skhtab(10,maxtab,maxtyp,maxtyp)
      real*8 skstab(10,maxtab,maxtyp,maxtyp)
      real*8 skself(3,maxtyp),dr(maxtyp,maxtyp),x0,x1,x2,f0,f1,f2
      real*8 xh,hl
      common /sktab/ skhtab,skstab,skself,dr,dim
      common /lmax/ lmax(maxtyp)

      skspar = 0
      maxmax = max(lmax(i),lmax(j))
      minmax = min(lmax(i),lmax(j))

      if(maxmax .le. 1)then
        inu = 10
      else

      if(maxmax .le. 2)then
        if(minmax .le. 1)then
          inu = 9
        else
          inu = 6
        endif
      else

      if(minmax .le. 1)then
        inu = 8
      else
        if(minmax .le. 2)then
         inu = 4
        else
         inu = 1
        endif
      endif
      endif
      endif
C
C mxind = Maximaler Index bis zu dem der Spline weitergefuehrt wird
C
      mxind=dim(i,j)+(0.3/dr(i,j)-1.0)
      r = sqrt(r2)
      ind = r/dr(i,j)+1.0
      if(r2 .lt. 1e-8)then
        do in = 1,3 
          dd(in+10) = 1.0
        enddo
      else

C
C neu einfuegen
C
      if(ind+2 .gt. dim(i,j))then
C
C free cubic spline
C
      if(ind .lt. dim(i,j))then
      x0=(dim(i,j)-3)*dr(i,j)
      x1=x0+dr(i,j)
      x2=x1+dr(i,j)
      xh=r-x1
      hl=x2-x1
      do in = inu,10 
        f0 = skstab(in,dim(i,j)-2,i,j)
        f1 = skstab(in,dim(i,j)-1,i,j)
        f2 = skstab(in,dim(i,j),i,j)
        dd(in)=cubicspline(f0,f1,f2,x0,x1,xh,hl,dr(i,j))
      enddo

      else
      if( (ind .ge. dim(i,j)) .and. (ind .lt. mxind) )then
C
C 5th degree spline
C
      x0=(dim(i,j)-3)*dr(i,j)
      x1=x0+dr(i,j)
      x2=x1+dr(i,j)
      xh=r-(mxind-1)*dr(i,j)
      do in = inu,10 
        f0 = skstab(in,dim(i,j)-2,i,j)
        f1 = skstab(in,dim(i,j)-1,i,j)
        f2 = skstab(in,dim(i,j),i,j)
        dd(in)=spline5th(f0,f1,f2,x0,x1,x2,xh,dr(i,j),mxind)
      enddo

      else
C
C zero
C
      do in = inu,10 
        dd(in) = 0.0
      enddo
C
C rem Ende
C
      endif
      endif

      else

      grdr = (r-(ind-1.)*dr(i,j))/dr(i,j)
      do in = inu,10 
        f0 = skstab(in,ind,i,j)
        f1 = skstab(in,ind+1,i,j)
        f2 = skstab(in,ind+2,i,j)
        dd(in) = f0 + (f1-f0)*grdr + (f2+f0-2.0*f1)*grdr*(grdr-1.0)/2.0
      enddo
      endif
      endif

      end
C
C looking for the hailton matrix parameters
C  use Newton's approximation for the calculation of elements
C  and in the end a spline (same routine as for ovelapp)
C
C**********************************************************************
      integer function skhpar(i,j,r2,dd)
C**********************************************************************
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
      integer i,j
      real*8 r2,dd(13)
      integer dim(maxtyp,maxtyp),mxind
      real*8 skhtab(10,maxtab,maxtyp,maxtyp) 
      real*8 skstab(10,maxtab,maxtyp,maxtyp)
      real*8 skself(3,maxtyp),dr(maxtyp,maxtyp),x0,x1,x2,f0,f1,f2
      real*8 xh,hl
      common /sktab/ skhtab,skstab,skself,dr,dim
      common /lmax/ lmax(maxtyp)
      skhpar = 0
      maxmax = max(lmax(i),lmax(j))
      minmax = min(lmax(i),lmax(j))

      if(maxmax .le. 1)then
        inu = 10
      else
        if(maxmax .le. 2)then
        if(minmax .le. 1)then
           inu = 9
        else
           inu = 6
        endif
        else
         if(minmax .le. 1)then
           inu = 8
         else
           if(minmax .le. 2)then
             inu = 4
           else
             inu = 1
           endif
         endif
      endif
      endif
C
C mxind = Maximaler Index bis zu dem der Spline weitergefuehrt wird
C
      mxind=dim(i,j)+(0.3/dr(i,j)-1.0)
      r = sqrt(r2)
      ind = r/dr(i,j)+1.0

      if(r2 .lt. 1e-8)then
        do in = 1,3 
          dd(in+10) = skself(in,i)
        enddo
      else

      if(ind+2 .gt. dim(i,j))then
C
C neu einfuegen
C
      if(ind .lt. dim(i,j))then
C
C free cubic spline
C
      x0=(dim(i,j)-3)*dr(i,j)
      x1=x0+dr(i,j)
      x2=x1+dr(i,j)
      xh=r-x1
      hl=x2-x1
      do in = inu,10 
        f0 = skhtab(in,dim(i,j)-2,i,j)
        f1 = skhtab(in,dim(i,j)-1,i,j)
        f2 = skhtab(in,dim(i,j),i,j)
        dd(in)=cubicspline(f0,f1,f2,x0,x1,xh,hl,dr(i,j))
      enddo

      else

      if( (ind .ge. dim(i,j)) .and. (ind .lt. mxind) )then
C
C 5th degree spline
C
      x0=(dim(i,j)-3)*dr(i,j)
      x1=x0+dr(i,j)
      x2=x1+dr(i,j)
      xh=r-(mxind-1)*dr(i,j)
      do in = inu,10 
        f0 = skhtab(in,dim(i,j)-2,i,j)
        f1 = skhtab(in,dim(i,j)-1,i,j)
        f2 = skhtab(in,dim(i,j),i,j)
        dd(in)=spline5th(f0,f1,f2,x0,x1,x2,xh,dr(i,j),mxind)
      enddo

      else
C
C zero
C
      do in = inu,10 
        dd(in) = 0.0
      enddo
C
C rem Ende
C
      endif
      endif
      else

      grdr = (r-(ind-1.0)*dr(i,j))/dr(i,j)
      do in = inu,10 
        f0 = skhtab(in,ind,i,j)
        f1 = skhtab(in,ind+1,i,j)
        f2 = skhtab(in,ind+2,i,j)
        dd(in) = f0 + (f1-f0)*grdr + (f2+f0-2.0*f1)*grdr*(grdr-1.0)/2.0
      enddo

      endif
      endif
      end
C**********************************************************************
      real*8 function cubicspline(f0,f1,f2,x0,x1,xh,hl,dr)
C**********************************************************************
      implicit none
      real*8 f0,f1,f2,x0,x1,xh,hl,dr
      real*8 f1abl,f2abl,a,b,c,d

      f2abl = (f2 + f0 - 2.0*f1)/(dr*dr)
      f1abl = (f1 - f0)/dr + 0.5*f2abl*(x1 - x0)

      a = f1
      b = f1abl
      c = f2abl/2.0
      d = (f2 - a)/(hl*hl*hl) - b/(hl*hl) - c/hl
      cubicspline = a + b*xh + c*xh*xh + d*xh*xh*xh
      end
C**********************************************************************
      real*8 function spline5th(f0,f1,f2,x0,x1,x2,xh,dr,mxind)
C**********************************************************************
      implicit none
      real*8 f0,f1,f2,x0,x1,x2,xh,dr
      integer mxind
      real*8 hl,f1abl,f2abl,a,b,c,d,hsp,isp,jsp

      f2abl=(f2+f0-2.0*f1)/(dr*dr)
      f1abl=(f1-f0)/dr+0.5*f2abl*(x1-x0)

      a=f1
      b=f1abl
      c=f2abl/2.0

      hl=x2-x1
      d=(f2-a)/(hl*hl*hl)-b/(hl*hl)-c/hl

      f1abl=b+2.0*c*hl+3.0*d*hl*hl
      f2abl=2.0*c+6.0*d*hl
      hl=x2-(mxind-1)*dr

      hsp=10.0*f2/(hl*hl*hl)-4.0*f1abl/(hl*hl)+f2abl/(2.0*hl)
      isp=-15.0*f2/(hl*hl*hl*hl)+7.0*f1abl/(hl*hl*hl)-f2abl/(hl*hl)
      jsp=6.0*f2/(hl*hl*hl*hl*hl)-3.0*f1abl/(hl*hl*hl*hl)+f2abl/(2.0*hl*
     *hl*hl)
      hl=xh*xh*xh
      spline5th = (hsp + isp*xh + jsp*xh*xh)*hl
      end
