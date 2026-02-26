C Output from Public domain Ratfor, version 1.0
C     SUBROUTINE OUT_EIGENVECTORS
C     =================
C
C     Copyright 1997 by Adam Gali
C
C *********************************************************************
C
C     PROGRAM CHARACTERISTICS
C     -----------------------
C
C print eigenvectors prints out the eigenvalues and the contributed    
C eigenvectors in ascending order in file 'EVC.DAT'. 
C No effects to the calculation!!
C
C PARAMETERS:
C a      r() the matrix of eigenvectors
C ev     r() the matrix of eigenvalues
C occ    r() vector of occupation number
C ind    i() index of orbitals possessed by one atom
C nn       i number of atoms
C
C *********************************************************************
C
      subroutine outeigenvectors(a,ev,occ,ind,nn)
C################################
C
C maxima definitions
C
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
C
      integer ind(nndim+1), izp(nndim), nn
      integer i,j,k,l
      integer neig,filled,hfilled
C  number of eigenstates, filled levels and half occupied states
      real*8 a(mdim,mdim), ev(mdim), occ(mdim), nel
      character*4 orbital(9)
      common /izp/ izp, nel, nbeweg
C              atom types, number of electrons, number of movable atoms
      data orbital/'S   ', 'Px  ', 'Py  ', 'Pz  ', 'Dxx ', 'Dxy ', 'Dyy 
     *', 'Dyz ', 'Dzz '/
CThe last index points to the number of last eigenstates
C
      neig = ind(nn+1)
      filled = 0
      hfilled = 0
      do i = 1,neig 
        if(occ(i) .gt. 1)then
          filled = filled+1
        else
         if(occ(i) .gt. 0)then
           hfilled = hfilled+1
         endif
        endif
      enddo

      open (1,file='EVC.DAT',status='unknown')
      open (2,file='MOCOEFF.DAT',status='unknown')
      rewind 1
      write (1,'(20x,a)') 'THE ONE-ELECTRON EIGENVALUES AND EIGENVECTORS
     *'
      write (1,*)
      write (1,*)
      write (1,*) 'NUMBER OF FILLED LEVELS: ',filled

      if(hfilled .gt. 0)then
        write (1,*) 'NUMBER OF HALF OCCUPIED STATES: ',hfilled
      endif
      write (1,*)
      write (1,*)

      do i=1,neig
        write(1,13) 'THE ',i,'. EIGENVALUE:',ev(i),
     &              ' H',ev(i)*27.2116,' eV'
        write(1,*) 'Atom No.   Atom Type'
        do j=1,nn
          write(1,*) '  ',j,'         ',izp(j)
          k=0
          do l=ind(j)+1, ind(j+1)
            k=k+1
            write(1,'(a,f8.5)') orbital(k),a(l,i)
            write(2,'(a,f15.10)') orbital(k), a(l,i)
          enddo
        enddo
       write(1,*)
      enddo

      endfile 1
      close (1)
      close (2)

C format specification
C
13    format(1x,a,i4,a,f15.8,a,4x,f12.5,a)
      end

C====================================================================================
      subroutine outspec(nn,nmov,ndim,ind,lmax,izp,ev,occ,efermi,  qmat,
     *qmulli,qtot,dipol,dipabs,iprint,iMDstep)
C
C maxima definitions
C
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
C
      integer nn,ndim,ind(nndim+1),lmax(maxtyp),izp(nndim)
      real*8 ev(mdim),occ(mdim),efermi,qmat(nndim),qmulli(mdim)
      real*8 qtot,dipol(3),dipabs
      integer i

C
C   dump  MO energies and occupations
C


      if (iprint.ge.3) then 
        write(125,*)' M.O. energies & occup. at  MD step' , iMDstep,' :'
        write (125,'(" Fermi level= ",f20.12)') efermi
        write (125,'(2(1x,f20.12))') (ev(i),occ(i),i = 1,ndim)
        endif



      if (iprint.ge.1) then 
        write(123,*)
        write(123,*)'#================================================'
        write(123,*)' MD step : ',iMDstep
C
C print diple moments
C
        write (123,'(a)') '#Total charge'
        write (123,'(a)') '#------------'
        write (123,'(1x,f20.12)') qtot
        write (123,'(a)') ' '
        write (123,'(a)') '#Mulliken dipole moment [Debye]'
        write (123,'(a)') '#------------------------------'
        write (123,'(3(1x,f20.12))') (dipol(i), i= 1,3)
        write (123,'(a)') ' '
        write (123,'(a)') '#Norm of dipole moment [Debye]'
        write (123,'(a)') '#------------------------------'
        write (123,'(1x,f20.12)') dipabs
        write (123,'(a)')
        endif
C
C  print Mulliken  charges
C
      
        write (123,'(a)') '#MULLIKEN CHARGE'
        write (123,'(a)') '#CHARGE PER ATOM / PER ATOMIC STATES'
        write (123,'(a)') '#ATOM No.'
        write (123,'(a,2x,i4)') '#',nmov
        write (123,'(a,2x,i4)') '#',nn
        do23014 n = 1,nn 
          ind1 = ind(n)+1
          ind2 = ind1+lmax(izp(n))**2-1
          write (123,98) n,qmat(n),(qmulli(i), i= ind1,ind2)
23014 continue
       


98    format(1x,i4,11(1x,f12.6))
      end

C=================================================================================
C
C Output of coordinates , velocities and locked atoms
C
C Yasuhito Ohta: xold ---> velo
C
      subroutine outcoord(name, namev, nn, period, ntype, izp, x, velo, 
     *boxsiz, atnames,mode,iMDstep)
C
C maxima definitions
C
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
C
      integer mode
      character*65 name,namev
      character*70 atnames
      character tmpchr
      logical period
      integer nn, ntype, izp(nndim)
C
C Yasuhito Ohta: xold(3*NNDIM) -> velo(3*NNDIM)
C
      real*8 x(3,nndim),boxsiz(3,3),velo(3*nndim)
      integer i,j
      real*8 xhlp(3)
C      ME: for constraints like in CHARMM
      integer lcc,jcc(nndim),atcc1(nndim),atcc2(nndim)
      real*8 rccval,kccval,rccinc
      common /cconstr/ lcc,jcc,atcc1,atcc2,rccval,kccval,rccinc
C
C     print velocities 
c
      !open (1,file=name,status='unknown')
      !rewind 1
      n3 = 3*nn
      if( mode .eq.1 .or. mode .eq.2 )then
        !open(105,file=namev,status='unknown')
        !rewind 105
        write(122,'("MD step: ", I10)')iMDstep
        do 23018 i=1,nn 
          !write(105,'(3F12.6)') velo(3*i-2),velo(3*i-1),velo(3*i)
          write(122,'(3F12.6)') velo(3*i-2),velo(3*i-1),velo(3*i)
23018     continue
        !close (105)
        endif
      if(period)then
        tmpchr='S'
        if(lcc.ne.0) tmpchr='SM'
      else
        tmpchr='C'
        if(lcc.ne.0) tmpchr='M'
        endif
c
c  print coordinates
c


!1000   format(A, I6,"  swer") 
      !write(1,'(I4,X,A1)') nn,tmpchr
      !write(1,'(A)') atnames
      write(121,'(I4,X,A1)') nn,tmpchr
      write(121,'(A50," MD step: ", I10)') atnames,iMDstep
      do 23026 i=1,nn 
        xhlp(1)=x(1,i)*0.529177
        xhlp(2)=x(2,i)*0.529177
        xhlp(3)=x(3,i)*0.529177
        write(121,14) i,izp(i),(xhlp(j),j=1,3)
23026   continue



      if(period)then
        write(1,15) 0.0,0.0,0.0
        do 23032 i=1,3 
          do 23034 j=1,3 
            xhlp(j)=boxsiz(i,j)*0.529177
23034     continue
          !write(1,15)(xhlp(j),j=1,3)
          write(121,15)(xhlp(j),j=1,3)
23032     continue
C     c ME: constraints like in CHARMM
        endif
      if(lcc.ne.0)then
        write(121,'(I3,3F12.6)') lcc,kccval/0.529177d0,
     $           rccval*0.529177d0,rccinc*0.529177d0
        write(121,'(100I3)')(jcc(i),atcc1(i),atcc2(i),i=1,lcc)
        endif
14    format(i4,1x,i3,3(1x,f12.6))
15    format(3(1x,f12.6))
      !endfile 1
      !close (1)

      end


C================================================================================
C  print forces now
C


      subroutine outforces(nbeweg,izp,gr,xm)
      implicit real*8 (a-h,o-z)
      integer nbeweg,izp(*)
      real*8 gr(*),xm(*)
      integer i,j
C
C Output of energy, old coordinates and forces for vibrational modes
C
C FD=124
!      open (2,file='FRC.DAT',status='unknown')
!      rewind 2
      write(124,'(A)') '# INTERATOMIC FORCES AT ATOM N0. I'
      write(124,'(A)') '# ATOM N0.'
      write(124,'(A)') '# FORCE_X, FORCE_Y, FORCE_Z'
      write(124,'(A)') '#'
      write(124,'(A)') '#'
      do23038 i = 1,nbeweg 
        write (124,'(i5,3(1x,e19.10))') i,(-gr(3*i+j-3), j = 1,3)
23038   continue
23039 continue
      write(124,'(A)') 'rms/d:   0.0000   0.0000'
!      close (2)
C
C
C Output of number of atoms and masses for vibrational modes
C
!      open (2,file='MAS.DAT',status='unknown')
!      rewind 2
!      write (2,'(1x,i5)') nbeweg
!C Yasuhito Ohta: xm(izp(i))*1822.887428d0 is added.
!      do23040 i = 1,nbeweg 
!      write (2,'(1x,f24.15,1x,f24.15)') xm(izp(i)),xm(izp(i))*1822.88742
!     *8d0
!23040 continue
!23041 continue
!      close (2)
      end
c==============================================================================
