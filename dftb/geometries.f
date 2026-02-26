C Output from Public domain Ratfor, version 1.0
C
C Calculate the inverse matrix (xinvbox) of the basis (boxsiz)
C
      subroutine inversebox
      implicit real*8 (a-h,o-z)
      real*8 xhlp
      real*8 boxsiz, xinvbox, xnullvec,skcut2
      logical period
      integer nlat
      common /box/ boxsiz(3,3),xinvbox(3,3),xnullvec(3), skcut2,nlat(3),
     *period
      xhlp=-boxsiz(1,3)*boxsiz(2,2)*boxsiz(3,1)+boxsiz(1,2)*boxsiz(2,3)*
     *boxsiz(3,1)
      xhlp=xhlp+boxsiz(1,3)*boxsiz(2,1)*boxsiz(3,2)
      xhlp=xhlp-boxsiz(1,1)*boxsiz(2,3)*boxsiz(3,2)
      xhlp=xhlp-boxsiz(1,2)*boxsiz(2,1)*boxsiz(3,3)
      xhlp=xhlp+boxsiz(1,1)*boxsiz(2,2)*boxsiz(3,3)
      xinvbox(1,1)=(-boxsiz(2,3)*boxsiz(3,2)+boxsiz(2,2)*boxsiz(3,3))/xh
     *lp
      xinvbox(2,1)=(boxsiz(2,3)*boxsiz(3,1)-boxsiz(2,1)*boxsiz(3,3))/xhl
     *p
      xinvbox(3,1)=(-boxsiz(2,2)*boxsiz(3,1)+boxsiz(2,1)*boxsiz(3,2))/xh
     *lp
      xinvbox(1,2)=(boxsiz(1,3)*boxsiz(3,2)-boxsiz(1,2)*boxsiz(3,3))/xhl
     *p
      xinvbox(2,2)=(-boxsiz(1,3)*boxsiz(3,1)+boxsiz(1,1)*boxsiz(3,3))/xh
     *lp
      xinvbox(3,2)=(boxsiz(1,2)*boxsiz(3,1)-boxsiz(1,1)*boxsiz(3,2))/xhl
     *p
      xinvbox(1,3)=(-boxsiz(1,3)*boxsiz(2,2)+boxsiz(1,2)*boxsiz(2,3))/xh
     *lp
      xinvbox(2,3)=(boxsiz(1,3)*boxsiz(2,1)-boxsiz(1,1)*boxsiz(2,3))/xhl
     *p
      xinvbox(3,3)=(-boxsiz(1,2)*boxsiz(2,1)+boxsiz(1,1)*boxsiz(2,2))/xh
     *lp
      end
C
C find the lenght of the shortest vertex in your supercell 
C
      real*8 function shortvertice(boxsiz)
      real*8 yhlp,boxsiz(3,3),testbox(6)
      testbox(1)=boxsiz(1,1)**2+boxsiz(1,2)**2+boxsiz(1,3)**2
      testbox(2)=boxsiz(2,1)**2+boxsiz(2,2)**2+boxsiz(2,3)**2
      testbox(3)=boxsiz(3,1)**2+boxsiz(3,2)**2+boxsiz(3,3)**2
      testbox(4)=(boxsiz(1,1)-boxsiz(2,1))**2+(boxsiz(1,2)-boxsiz(2,2))*
     **2+  (boxsiz(1,3)-boxsiz(2,3))**2
      testbox(5)=(boxsiz(1,1)-boxsiz(3,1))**2+(boxsiz(1,2)-boxsiz(3,2))*
     **2+ (boxsiz(1,3)-boxsiz(3,3))**2
      testbox(6)=(boxsiz(3,1)-boxsiz(2,1))**2+(boxsiz(3,2)-boxsiz(2,2))*
     **2+  (boxsiz(3,3)-boxsiz(2,3))**2
      yhlp=min(testbox(1),testbox(2),testbox(3),testbox(4),testbox(5),te
     *stbox(6))
      shortvertice=sqrt(yhlp)
      end
C
C Put Atom at (x,z,z) back into supercell
C
      subroutine coordback(x,y,z,xo,yo,zo)
C extended to tranform the old xyz coordinates (xo,yo,zo) with the new one,
C to avoid TOTALLY ERRONEOUS  velocities 
      implicit real*8 (a-h,o-z)
      real*8 x,y,z, xo, yo, zo, tx, ty , tz
      real*8 xx1,xy1,xz1
      real*8 boxsiz, xinvbox, xnullvec,skcut2
      logical period
      integer nlat
      common /box/ boxsiz(3,3),xinvbox(3,3),xnullvec(3), skcut2,nlat(3),
     *period
      xx1=(x-xnullvec(1))*xinvbox(1,1)+(y-xnullvec(2))*xinvbox(2,1)+(z-x
     *nullvec(3))*xinvbox(3,1)
      xy1=(x-xnullvec(1))*xinvbox(1,2)+(y-xnullvec(2))*xinvbox(2,2)+(z-x
     *nullvec(3))*xinvbox(3,2)
      xz1=(x-xnullvec(1))*xinvbox(1,3)+(y-xnullvec(2))*xinvbox(2,3)+(z-x
     *nullvec(3))*xinvbox(3,3)
      if(xx1.gt.0.5)then
      xx1=xx1-1.0
      endif
      if(xx1.lt.-0.5)then
      xx1=xx1+1.0
      endif
      if(xy1.gt.0.5)then
      xy1=xy1-1.0
      endif
      if(xy1.lt.-0.5)then
      xy1=xy1+1.0
      endif
      if(xz1.gt.0.5)then
      xz1=xz1-1.0
      endif
      if(xz1.lt.-0.5)then
      xz1=xz1+1.0
      endif
      tx=xx1*boxsiz(1,1)+xy1*boxsiz(2,1)+xz1*boxsiz(3,1)+xnullvec(1)
      xo=xo+(tx-x)
      x=tx
      ty=xx1*boxsiz(1,2)+xy1*boxsiz(2,2)+xz1*boxsiz(3,2)+xnullvec(2)
      yo=yo+(ty-y)
      y=ty
      tz=xx1*boxsiz(1,3)+xy1*boxsiz(2,3)+xz1*boxsiz(3,3)+xnullvec(3)
      zo=zo+(tz-z)
      z=tz
      end
C
Cget the three summation limits for the matrix construction in
Cthe gamma point approximation
Coutput: nlat(3)
C
      subroutine gamma_summind(slkcutoff)
      implicit none
      real*8 slkcutoff
      real*8 u(3),v(3),w(3),helpv(3),l,lu,lv,lw
      integer nu,nv,nw
      real*8 boxsiz, xinvbox, xnullvec,skcut2
      logical period
      integer nlat
      common /box/ boxsiz(3,3),xinvbox(3,3),xnullvec(3), skcut2,nlat(3),
     *period
Cget vectors and lengths
      u(1) = boxsiz(1,1)
      u(2) = boxsiz(1,2)
      u(3) = boxsiz(1,3)
      v(1) = boxsiz(2,1)
      v(2) = boxsiz(2,2)
      v(3) = boxsiz(2,3)
      w(1) = boxsiz(3,1)
      w(2) = boxsiz(3,2)
      w(3) = boxsiz(3,3)
      lu = sqrt(u(1)**2 + u(2)**2 + u(3)**2)
      lv = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      lw = sqrt(w(1)**2 + w(2)**2 + w(3)**2)
Csee, whether length of u is shorter than projected side of v
      call cross(u,v,helpv)
      l = min(lu,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lv)
Csee, whether length of u is shorter than projected side of w
      call cross(u,w,helpv)
      l = min(lu,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lw)
Cset nlat(1) according to slkcutoff. At least 1.
      nlat(1) = max(int(2*slkcutoff/l),1)
Cfor v,w do the same as for u...
      call cross(v,u,helpv)
      l = min(lv,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lu)
      call cross(v,w,helpv)
      l = min(lv,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lw)
      nlat(2) = max(int(2*slkcutoff/l),1)
      call cross(w,u,helpv)
      l = min(lw,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lu)
      call cross(w,v,helpv)
      l = min(lw,sqrt(helpv(1)**2 + helpv(2)**2 + helpv(3)**2)/lv)
      nlat(3) = max(int(2*slkcutoff/l),1)
      end
