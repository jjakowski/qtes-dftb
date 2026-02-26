c==============================================================================
c      get the hubbard contribution to the H matrix elements
c
c      INPUT:
c      INTEGER nat           number of atoms in cell
c      REAL*8  qmat(*)      mullikan charges
c      REAL*8  qzero(*)      neutral atomar charges
c      REAL*8  uhub(*)     list of hubbard parameters
c      REAL*8  basis(3,3)    basis of cell
c      REAL*8  rat(3,*)      atomic coordinates
c      LOGICAL period        .TRUE. if periodic boundary conditions
c      INTEGER izp           
c      INTEGER niter         step in scf cycle if == 1: build up gammamat          
c      OUTPUT:
c      REAL*8    gammamat(*,*)  matrix containing the Ewald & short range terms
c      REAL*8    shift(NNDIM) array contains shifts for hamilton matrix elements
c==============================================================================

       SUBROUTINE HAMILSHIFT(nat,rat,atomtype,basis,period,uhubb,
     &                       qmat,qzero,niter,gammamat,shift)
       IMPLICIT NONE
       INCLUDE 'maxima.inc'
       integer nat,niter,i,j
       integer atomtype(NNDIM)
       real*8 shift(NNDIM),qmat(NNDIM),qzero(MAXTP),uhubb(MAXTP)
       real*8 rat(3,*),basis(3,3),gammamat(NNDIM,NNDIM)
       logical period
       EXTERNAL GAM12
       real*8 Qdiff(NNDIM)
       real*8 tmpvalue
       
       common /varhelp/ Qdiff
       

c
c       build up gammamatrix before the first iteration. This matrix is
c       then reused at the following iterations
c
        IF (niter .eq. 1) THEN
        CALL gammamatrix(nat,rat,atomtype,basis,period,uhubb,gammamat)
        ENDIF

c       set net charges Qdiff = qmat - qzero

       DO i=1,nat
        qdiff(i) = qmat(i) - qzero(atomtype(i))       
       END DO

     
c
c      Calculate atomic hamilton shift (=sum over gamma * charges)
c      

       DO i=1,nat 

       shift(i) = 0.0
       DO j=1,nat
         if (j .gt. i) then
          tmpvalue = gammamat(j,i)
         else
          tmpvalue = gammamat(i,j)
         endif
         shift(i) = shift(i) + Qdiff(j)*tmpvalue
       END DO

       END DO
  
       END



c============================================================================
c      get the gamma contribution to the gradient
c
c
c      INPUT:
c      INTEGER nat           number of atoms in cell
c      INTEGER nmov          number of movable atoms in cell
c      REAL*8    qmat(*)     mullikan charges
c      REAL*8    qzero(*)    neutral atomar charges
c      REAL*8    basis(3,3)  basis of cell
c      REAL*8    rat(3,*)    atomic coordinates
c      INTEGER atomtype(*)   list of atomic types             
c      LOGICAL period        .TRUE. if periodic boundary conditions
c      REAL*8    derivx(*,*)  scratch space matrix for x derivative 
c      REAL*8    derivy(*,*)  scratch space matrix for y derivative 
c      REAL*8    derivz(*,*)  scratch space matrix for z derivative 
c
c      OUTPUT:
c      REAL*8    gmgrd(3,*)gamma contribution to the gradient
c
c============================================================================

       SUBROUTINE GAMMAGRAD(nat,nmov,rat,atomtype,basis,period,uhubb,
     &                      qmat,qzero,derivx,derivy,derivz,gmgrd)
       IMPLICIT NONE
       INCLUDE 'maxima.inc'
       integer nat,nmov,atomtype(*)
       real*8 qmat(NNDIM),qzero(MAXTP),qdiff(NNDIM)
       real*8 uhubb(MAXTP),rat(3,*),basis(3,3)
       real*8 derivx(NNDIM,NNDIM),derivy(NNDIM,NNDIM)
       real*8 derivz(NNDIM,NNDIM),deriv(3),tmpderiv(3)
       logical period

       integer i,k,l
       real*8 gmgrd(3,*)
       common /varhelp/ Qdiff
 
       DO k=1,nat
        DO i=1,3
         gmgrd(i,k) = 0.0
        END DO
       END DO

     
c      construct gamma derivatives
       CALL gammamatrix1(nat,rat,atomtype,basis,period,uhubb,
     &                   derivx,derivy,derivz)
        
c      set net charges Qdiff = qmat - qzero

       DO i=1,nat
        qdiff(i) = qmat(i) - qzero(atomtype(i))       
       END DO

        
       DO k=1,nmov 
     
       deriv(1) = 0
       deriv(2) = 0
       deriv(3) = 0

       DO i=1,nat
         if (i .gt. k) then
          tmpderiv(1) = derivx(i,k)
          tmpderiv(2) = derivy(i,k)
          tmpderiv(3) = derivz(i,k)
         else if (i .eq. k) then
          tmpderiv(1) = 0.0
          tmpderiv(2) = 0.0
          tmpderiv(3) = 0.0
         else
c         shortrange1(rrmuind - ri)  = -shortrange1(ri - rrmuind)
          tmpderiv(1) = -derivx(k,i)
          tmpderiv(2) = -derivy(k,i)
          tmpderiv(3) = -derivz(k,i)
         endif

         deriv(1) = deriv(1) + Qdiff(i)*tmpderiv(1)
         deriv(2) = deriv(2) + Qdiff(i)*tmpderiv(2)
         deriv(3) = deriv(3) + Qdiff(i)*tmpderiv(3)
       END DO

      DO l = 1,3 
        gmgrd(l,k) = gmgrd(l,k) + Qdiff(k)*deriv(l)
      END DO

      END DO

      END
