c=============================================================================
c
c      gamma resulting from exact Coulomb interaction of normalized
c      exp(-a*r) charge distribution
c      Attention: subroutine gamsub needed
c
c      input:  r:     distance
c         uhub1: Hubbard parameter orbital 1
c         uhub2: Hubbard parameter orbital 2
c      output: gval:  gamma12 function value
c
c=============================================================================

       subroutine gam12(r,uhub1,uhub2,gval)
       IMPLICIT NONE
       REAL*8 zero
       parameter(zero=1.0d-4)
       REAL*8 gval,a1,a2,src,avg,uhub1,uhub2,rrc,rrc3
       REAL*8 val12,val21,drv12,drv21,r,fac,fac2,efac
c     neu fuer H-bonds
       real*8 fhbond,uhubh,cut1,cut2,rcuth
       real*8 kl1,kl2,kl3,kkk
c        open(111,file='switch') 
c        read(111,*)  kl1
c        close(111)
c   parallel  trajectories execution 
      include 'mpiio2.inc'     !ninfile, noutfile
        kl1 = 4.0
        uhubh= 0.4195007d0
        fhbond = 1.0
c        if((uhub1.eq.uhubh).and.(uhub2.ne.uhubh)) then
c         write(*,*) 'uhubh, gamma on',uhubh,uhub1,uhub2
        if ((uhub1.eq.uhubh).or.(uhub2.eq.uhubh)) then
         fhbond=  exp(-(((uhub1+uhub2)/2)**kl1)*r**2)
        endif
c        if (uhub2.eq.uhubh) then
c        if((uhub2.eq.uhubh).and.(uhub1.ne.uhubh)) then
c         fhbond=  exp(-(((uhub1+uhub2)/2)**kl1)*r**2)
c        endif
c  end neu hbond
c gamma  besteht aus einem 1/r Term, und etwas, 
c was fuer r=0 gegen Hubbard geht:
c multipliziere einfach den zweiten Term mit fhbond!
       gval= 0.0
       a1= 3.2*uhub1
       a2= 3.2*uhub2
       IF (a1+a2 .lt. zero) THEN
        RETURN
       ENDIF
       src= 1.0/(a1+a2)
       fac= a1*a2*src
       avg= 1.6*(fac+fac*fac*src)
       IF (r .lt. zero) THEN
           gval= 0.3125*avg
       ELSE
         rrc= 1.0/r
         rrc3= rrc*rrc*rrc
         IF (abs(a1-a2) .lt. 1.0d-5) THEN
           fac= avg*r
           fac2= fac*fac
           efac= exp(-fac)/48.0
           gval= (1.0-fhbond*(48.0+33*fac+fac2*(9.0+fac))*efac)*rrc
         ELSE
           call gamsub(a1,a2,r,rrc,val12,drv12)
           call gamsub(a2,a1,r,rrc,val21,drv21)
           gval= rrc-fhbond*val12-fhbond*val21
         ENDIF
       ENDIF
       RETURN
       END



c=============================================================================
c
c      derivative of gamma resulting from exact Coulomb interaction 
c      of normalized  exp(-a*r) charge distribution
c      Attention: subroutine gamsub needed
c
c      input:  r:     distance
c         uhub1: Hubbard parameter orbital 1
c         uhub2: Hubbard parameter orbital 2
c      output: gdrv:  ((d gamma12)/(d r)) * 1/r
c
c=============================================================================

       subroutine gam121(r,uhub1,uhub2,gdrv)
       IMPLICIT NONE
       REAL*8 zero
       parameter(zero=1.0d-4)
       REAL*8 gdrv,a1,a2,src,avg,uhub1,uhub2,rrc,rrc3
       REAL*8 val12,val21,drv12,drv21,r,fac,fac2,efac
c     neu fuer H-bonds

       real*8 fhbond,uhubh,cut1,cut2,rcuth,cutab(5,5)
       real*8 dhbond
       real*8 kl1,kl2,kl3,kkk
c         open(111,file='switch') 
c         read(111,*)  kl1
c        close(111)
        kl1=4.0
        uhubh = 0.4195007d0
         fhbond = 1.0
         dhbond=0.0
         if((uhub1.eq.uhubh).or.(uhub2.eq.uhubh)) then
c        if((uhub1.eq.uhubh).and.(uhub2.ne.uhubh)) then
        fhbond= exp(-(((uhub1+uhub2)/2)**kl1)*r**2)
        dhbond= -2*fhbond*(((uhub1+uhub2)/2)**kl1)/r
        endif
c         if(uhub2.eq.uhubh) then
c        if((uhub2.eq.uhubh).and.(uhub1.ne.uhubh)) then
c        fhbond= exp(-(((uhub1+uhub2)/2)**kl1)*r**2)
c        dhbond= -2*fhbond*(((uhub1+uhub2)/2)**kl1)/r
c        endif
c  end neu hbond
c gamma  besteht aus einem 1/r Term, und etwas,
c was fuer r=0 gegen Hubbard geht:
c multipliziere einfach den zweiten Term mit fhbond!

       gdrv= 0.0
       a1= 3.2*uhub1
       a2= 3.2*uhub2
       IF (a1+a2 .lt. zero) THEN
        RETURN
       ENDIF
       src= 1.0/(a1+a2)
       fac= a1*a2*src
       avg= 1.6*(fac+fac*fac*src)
       IF (r .lt. zero) THEN
       ELSE
         rrc= 1.0/r
         rrc3= rrc*rrc*rrc
         IF (abs(a1-a2) .lt. 1.0d-5) THEN
           fac= avg*r
           fac2= fac*fac
           efac= exp(-fac)/48.0
           gdrv= 
     c      -(1.0-fhbond*(48.0+48*fac+fac2*(24.0+7*fac+fac2))*efac)*rrc3
     c      + dhbond*(48.0+33*fac+fac2*(9.0+fac))*efac*rrc
         ELSE
           call gamsub(a1,a2,r,rrc,val12,drv12)
           call gamsub(a2,a1,r,rrc,val21,drv21)
           gdrv= -rrc3-fhbond*(drv12+drv21)*rrc
     c        -dhbond*(val12+val21)
         ENDIF
       ENDIF
       RETURN
       END



c===========================================================================
c
c      auxiliary routine needed by gam12 and gam121
c
c      input   a:    alpha1
c              b:    alpha2
c              r:    distance
c              rrc:  1/distance
c      output: gval: function value
c              gdrv: function derivative
c
c===========================================================================

       subroutine gamsub(a,b,r,rrc,gval,gdrv)
       IMPLICIT NONE
       REAL*8 a,a2,b,b2,b4,b6,drc,drc2,r,efac,rrc,fac,gval,gdrv
       a2= a*a
       b2= b*b
       b4= b2*b2
       b6= b4*b2
       drc= 1.0/(a2-b2)
       drc2=drc*drc
       efac= exp(-a*r)
       fac= (b6-3*a2*b4)*drc2*drc*rrc
       gval= efac*(0.5*a*b4*drc2-fac)
       gdrv= -a*gval+efac*fac*rrc
       RETURN
       END
