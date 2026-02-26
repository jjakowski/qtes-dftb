C
C      PROGRAM DYLCAO
C     ================
C
C     Copyright 1991 by Peter Blaudeck, Dirk Porezag 
C
C *********************************************************************
C
C maxima definitions
C
!      program dftb
!      subroutine pbdylcao(ibohm,natoms,iat,CC,Frc,Energy,Emermin,
!     $      blabel,V)
      subroutine do_dftb(ibohm,natoms,iat,CC,Frc,Energy,Emermin,
     $      blabel,ntype,qnel,scfhelp,scftol, telec,V,IV)
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
C
      integer izp(nndim), lmax(maxtyp)
      integer dim(maxtyp,maxtyp), conat(nndim+1)
      integer numint(maxtyp,maxtyp)
      integer nexc,nrem
      real*8 x(3*nndim), xparam(3*nndim), xold(3*nndim), gr(3*nndim)
      real*8 norm
      real*8 skhtab(10,maxtab,maxtyp,maxtyp) 
      real*8 skstab(10,maxtab,maxtyp,maxtyp)
      real*8 skself(3,maxtyp), sr(maxtyp,maxtyp), nel, convec(3,nndim)
      real*8 coeff(6,maxint,maxtyp,maxtyp),xr(2,maxint,maxtyp,maxtyp)
      real*8 efkt(3,maxtyp,maxtyp),cutoff(maxtyp,maxtyp)
      real*8 xm(maxtyp), xmrc(maxtyp), rmax(maxtyp,maxtyp)
      real*8 ar(10,maxtyp,maxtyp),qmat(nndim), slkcutoff
      real*8 espin(maxtyp),ekin
      real*8 V(*)
      integer IV(*)

      character*65 infile, outfile, name, tmpchr, namev
      character*70 atnames
      character adummy5*5
      !logical mdtemp, mdfree, stdc, cgrmin, atomic, evector, scfhelp
      logical cgrmin, atomic, evector, scfhelp
      logical converge, mdcorr, dorelax, doecalc
      logical chrr,constr,ex
      logical lbfgs,lzmat,lgrad
      logical lexc
      external skhpar, skspar
C Yasuhito Ohta ----------------------
      parameter(avog=6.0221367D+23,togrm=9.1093897D-28)
      parameter(amu2au = 1.0d0/avog/togrm)
      parameter(pi = 3.1415926535897932d0)
c see http://en.wikipedia.org/wiki/Atomic_units
      parameter(au2sec=2.418884326505d-17)
      parameter(au2fs=au2sec*1.D15,fs2au=1.0d0/au2fs)
      parameter(au2ps=au2sec*1.D12,ps2au=1.d0/au2ps)
      parameter(rKb=3.16683D-6)
      real*8 amass(nndim),velo(3*nndim)
      real*8 gamma_lgv,random_frc
!------------
!     1: S,Hc
!     2: 2e
!     3: diag   
!     4: frc      
!     5:  total eglcao      
      real*8 timetab(10)
      
C-------------------------------------
CCHARMM
      real*8 xe(3,4*nndim),ze(4*nndim)
Cinteger nE,nEiter
      character*2 ext
C      ME: for constraints like in CHARMM
      integer lcc,jcc(nndim),atcc1(nndim),atcc2(nndim)
c      real*8 rccval,kccval,rccin
CDISPERSION
C integer subsys(NNDIM)
C real*8 pol(MAXTYP,4),ion(MAXTYP,4),Edis,C6(NNDIM)
      real*8 edis
      logical dispers
      common /disper/ edis, dispers
CEND DISPERSION
CCHARMM
      common /cconstr/ lcc,jcc,atcc1,atcc2   !,rccval,kccval,rccinc
      common /extchr/ xe, ze, ne, neiter ,icycle, ext
      common /sktab/ skhtab,skstab,skself,sr,dim
      common /spin/ espin
      common /spltab/ coeff,xr,efkt,cutoff,numint
      common /concgr/ convec, constr, conat
C
C##### 
C These data are read in by gettab. 
C Every pair of atoms requires a table containing the Slater-Koster 
C data. The first line consists of the distance stepwidth and the
C number of distance values. If the file describes homonuclear
C interactions, an additional second line contains the energies
C of the valence d, p and s electrons (diagonal elements of the
C hamiltonian). 
C All other lines contain 10 hamiltonian and 10 overlap numbers
C in the order dds ddp ddd pds pdp pps ppp sds sps sss.
C The first of these lines (r = 0) contains (in this order):
C element 01:    mass of the atom (0 for heteroatomic files)
C         02-09: polynomial coefficients for repulsive part
C         10:    range of repulsive part
C         11:    cutoff radius for next neighbour determination
C All energies in H = Hartrees = 4.358d-11 g*cm**2/s**2
C All masses in nucleon masses = 1.6603d-24 g
C
      common /lmax/ lmax
      common /izp/ izp, nel, nbeweg
      common /rdaten/ ar
      common /machine/ racc
      real*8 boxsiz, xinvbox, xnullvec,skcut2
      logical period
      integer nlat
      common /box/ boxsiz(3,3),xinvbox(3,3),xnullvec(3), skcut2,nlat(3),
     *period
      common /mcharge/ qzero(maxtyp), uhubb(maxtyp)
C
c   parallel  trajectories execution 
      real*8  CC(3,*), Frc(*)
      integer iat(*)
      logical  LOPF           ! if file is open  
      character*70 blabel     !  print label in GEOMS.DAT comment line
      include 'mpiio.inc'     ! ninfile, noutfile

C
C Get the machine accuracy
C
      racc = 1.0
23000 if((1.0+racc) .gt. 1.0)then
      racc = racc*0.5
      goto 23000
      endif
      racc = 2*racc
c----  JACEK's comments: 
c  DFTB init:
c      ntype : argument
c      lmax  : passed through common block from dftb_init
c      scfhelp,scftol, telec,qnel:  input arguments 
c      nel, izp,nbweg  :  common block (we overwrite it here)
c      nel             : float   passed  through  qnel               
c
       nbeweg   = 0        !  number of  movable atoms  
                           ! /=0: No constraints; >0: some atoms do not move/
       nel =  qnel         ! set total  number of electrons  nel is float
       chrr    =  .false.  ! read  initial charges (guess=read?) 
       ext     =  'NO'     ! external charges ?
       iprint  =   3
       dispers =  .false.  !  I  
!-------------------------------

C
C Enter relaxation mode and fmax and scftol. Valid relaxation modes are
C 6 ... Mulliken analysis and atomic energy calculation
      cgrmin = .false.
      atomic = .false. 
      evector = .false.
      constr =.false.
      lbfgs = .false.
      lzmat = .false. 
      lgrad = .true. 
      lexc = .false.
c
c  iprint =0   : print only:   geoms, velocs  
c          1   : frc, chr, 
c          2   : SPE
     
      if(scfhelp)then
         scftol=abs(scftol)
      else
         scftol=-abs(scftol)
      endif

        atomic = .true.

!      if(nbeweg.eq.0)then
!        write(istdout,*) 'No constraints set.'
!      else
!        write(istdout,*) 'Number of moveable atoms set to ',nbeweg
!      endif


c
c  copy geometry  to dftb from input arguments
C All lengths in a.u. = 5.292d-9 cm
c     
      if (ibohm.gt.0) then
         nn= natoms 
         ii=0
         do i =1,nn
           x(ii+1) = CC(1,i)
           x(ii+2) = CC(2,i)
           x(ii+3) = CC(3,i)
           ii=ii+3
           izp(i) = iat(i)
           enddo
         period =.false.
         atnames=blabel
        endif   


      if(nn .gt. nndim)then
        write(istdout,*)'dftb: nn > ', nndim
        goto 9
        endif
      if(nbeweg.eq.0) nbeweg=nn
      n3 = 3*nn
      nbew3 = 3*nbeweg

!      if (ibohm.eq.0) close (1)
C
C Check if ntype in limits of MAXTYP 
C
      ntype = 1
      do n = 1,nn 
        if(izp(n) .gt. ntype)then
         ntype = izp(n)
        endif
      enddo

      if(ntype .gt. maxtyp)then
        write(istdout,*)'dftb: ntype >', maxtyp
        goto 9
      endif

C
C Read in whether atom shall have a virtual atom at the same place,
C i.e. wether you want an additional s (or p) orbital or not
C
      n3 = 3*nn
      nbew3 = 3*nbeweg

C
C Read in Slater-Koster-table 
C
!      call gettab(ntype)
C
C nel is the net charge
C nel=1, one electron missing
C nel =-1 one electron more

      nel = -nel
      do i = 1,nn 
        nel = nel+qzero(izp(i))
      enddo

      !write(istdout,*) 'total number of electrons:', nel
C
C set initial charges

          do i = 1,nn 
            qmat(i) = qzero(izp(i))
          enddo

C
C Get masses, data for repulsive part and next neighbours 
C
      do i = 1,ntype 
        ar(1,i,i)=skhtab(1,1,i,i)
        do j = 1,ntype 
          do k = 2,10 
            ar(k,i,j) = skhtab(k,1,i,j)
          enddo
            rmax(i,j) = skstab(1,1,i,j)
        enddo

      enddo

CDISPERSION
      if(dispers)then
        call dispersionread(nn,ntype,izp,x)
      endif
CEND DISPERSION

C Check boxsize according to the cutoff in SLK-data:
C    the shortest distance of two vertices has to be twice the shortest
C    cutoff
C
      slkcutoff=0.0d0

      do i = 1,ntype 
        do j = 1,ntype 
         yhlp=sr(i,j)*dim(i,j)+0.3
         if(yhlp.gt.slkcutoff)then
           slkcutoff=yhlp
         endif
        enddo
      enddo

      skcut2=slkcutoff*slkcutoff

!      if(period)then
!        call gamma_summind(slkcutoff)
!        write(istdout,*) 'Number of lattice sums for 
!     &Gammapoint matrix:', nlat(1), nlat(2), nlat(3)
!        call initphi(boxsiz)
!      endif


C
C Do the work (total energy, gradient)
C   if idump > 0   then print
c
c   iprint:   input for print    control   
c   iprint1:  pritns only at certain MD steps      
    
        idump=1
        iprint1=iprint 
        iprint1= 0
        iMDstep =  1   ! this is obsolete, used by Jacek to print in eglcalio 

      e = 0.0
      eel = 0.0

      call eglcao(n3,telec,x,scftol,atomic,e,eel,gr,niter,
     &       evector,qmat,lgrad,lexc,nexc,nrem,entropy,timetab,iprint1,
     &       iMDstep)

      if  (ibohm.gt.0)then 
        Energy = E
        Emermin = E- entropy* telec
        do i=1,Natoms*3  
          Frc(i) = -gr(i)
          enddo
        endif

9     continue

      end
