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
!     lmax : passsed through common       
!      subroutine dftb_init_mpi(ntype,V,MDV)
      subroutine dftb_init_mpi(ntype)
      implicit real*8 (a-h,o-z)
      include 'mpif.h'
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
      real*8 skself(3,maxtyp), dr(maxtyp,maxtyp), nel, convec(3,nndim)
      real*8 coeff(6,maxint,maxtyp,maxtyp),xr(2,maxint,maxtyp,maxtyp)
      real*8 efkt(3,maxtyp,maxtyp),cutoff(maxtyp,maxtyp)
      real*8 xm(maxtyp), xmrc(maxtyp), rmax(maxtyp,maxtyp)
      real*8 ar(10,maxtyp,maxtyp),qmat(nndim), slkcutoff
      real*8 espin(maxtyp),ekin

      parameter (MDV=10000000)
      real*8   V(MDV)             !scratch array / buffer

      character*65 infile, outfile, name, tmpchr, namev
      character*70 atnames
      character adummy5*5
      !logical mdtemp, mdfree, stdc, cgrmin, atomic, evector, scfhelp
      logical cgrmin, atomic, evector, scfhelp
      logical converge, mdcorr, dorelax, doecalc
      logical constr,ex
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
      real*8 rccval,kccval,rccinc
CDISPERSION
C integer subsys(NNDIM)
C real*8 pol(MAXTYP,4),ion(MAXTYP,4),Edis,C6(NNDIM)
      real*8 edis
      logical dispers
      common /disper/ edis, dispers
CEND DISPERSION
CCHARMM
      common /cconstr/ lcc,jcc,atcc1,atcc2,rccval,kccval,rccinc
      common /extchr/ xe, ze, ne, neiter ,icycle, ext
      common /sktab/ skhtab,skstab,skself,dr,dim
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
!      real*8  CC(3,*), Frc(*)
!      integer iat(*)
      logical  LOPF           ! if file is open  
      character*70 blabel     !  print label in GEOMS.DAT comment line
      include 'mpiio.inc'     ! ninfile, noutfile
      
c  get MPI   id      
      call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
      write(*,*)'dftb_init_mpi, myid= ' , myid

C
C Read in Slater-Koster-table  for myid=0
C
      if (myid.eq.0)  then
        write(istdout,*) 'How manny different atoms: Ntype =?'
        read(istdin,*) ntype      
        write(istdout,*) ntype
        read(istdin,*)(lmax(i),i = 1,ntype)

C       now read files
        call gettab(ntype)
        write(*,*)myid,':  done with dftb_init'
        endif 


c prepare buffer  for  MPI transfer 

      maxt2 = maxtyp*maxtyp     ! currently: 36=6*6

c  pointers  
      intype   =  1                       ! Int
      ilmax    = intype    +     1        ! int
      idimm    = ilmax     +    maxtyp    ! Int
      inumint  = idimm     +    maxt2     ! Int
      idr      = inumint   +    maxt2
      iskself  = idr       +    maxt2
      iespin   = iskself   +  3*maxtyp
      iuhubb   = iespin    +    maxtyp
      iqzero   = iuhubb    +    maxtyp
      iskhtab  = iqzero    +    maxtyp
      iskstab  = iskhtab   + 10*maxt2 *maxtab
! now splines
      icutoff  = iskstab   + 10*maxt2 *maxtab
      iefkt    = icutoff   +    maxt2
      ixr      = iefkt     +  3*maxt2
      icoeff   = ixr       +  2*maxt2 *maxint
      iend     = icoeff    +  6*maxt2 *maxint

      Mints    = idr  - intype                    ! how many  integers we store ??
      Mfloats  = iend - idr                       ! how many  floats

c  check if enough space to copy      
      call  tstcor(iend,mdv)

c now copy everything  to buffer
      if (myid.eq.0) then 
!        write(*,*)' pointers: intype, ilmax, idimm,inumint,idr'       
!        write(*,*)' pointers:', intype, ilmax, idimm,inumint,idr       
!
!        write(*,*)' pointers: iskself, iespin, iuhubb,iqzero ,iskhtab'  
!        write(*,*)' pointers:',iskself, iespin, iuhubb,iqzero ,iskhtab  
!
!        write(*,*)' pointers: iskstab, icutoff,iefkt, ixr, icoeff, iend'
!        write(*,*)' pointers:',iskstab, icutoff,iefkt, ixr, icoeff, iend


        call  cp_i2i(ntype,  V(intype),  1)
        call  cp_i2i(lmax,   V(ilmax),   maxtyp)
        call  cp_i2i(dim,    V(idimm),   maxt2)
        call  cp_i2i(numint, V(inumint), maxt2)
        call  cp_d2d(dr,     V(idr),     maxt2)
        call  cp_d2d(skself, V(iskself), 3*maxtyp)
        call  cp_d2d(espin,  V(iespin),  maxtyp)
        call  cp_d2d(uhubb,  V(iuhubb),  maxtyp)
        call  cp_d2d(qzero,  V(iqzero),  maxtyp)
        call  cp_d2d(skhtab, V(iskhtab), 10*maxt2 *maxtab)
        call  cp_d2d(skstab, V(iskstab), 10*maxt2 *maxtab)
        call  cp_d2d(cutoff, V(icutoff), maxt2)
        call  cp_d2d(efkt,   V(iefkt),   3*maxt2)
        call  cp_d2d(xr,     V(ixr),     2*maxt2 *maxint)
        call  cp_d2d(coeff,  V(icoeff),  6*maxt2 *maxint)

        endif

c  distribute it
        call MPI_Bcast(V(1),iend,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $               ierr)

c  ...  and finally copy from   buffer  to the right places after the transfer
      if (myid.gt.0) then
        call  cp_i2i(V(intype),  ntype,  1)
        call  cp_i2i(V(ilmax),   lmax,   maxtyp)
        call  cp_i2i(V(idimm),   dim,    maxt2)
        call  cp_i2i(V(inumint), numint, maxt2)
        call  cp_d2d(V(idr),     dr,     maxt2)
        call  cp_d2d(V(iskself), skself, 3*maxtyp)
        call  cp_d2d(V(iespin),  espin,  maxtyp)
        call  cp_d2d(V(iuhubb),  uhubb,  maxtyp)
        call  cp_d2d(V(iqzero),  qzero,  maxtyp)
        call  cp_d2d(V(iskhtab), skhtab, 10*maxt2 *maxtab)
        call  cp_d2d(V(iskstab), skstab, 10*maxt2 *maxtab)
        call  cp_d2d(V(icutoff), cutoff, maxt2)
        call  cp_d2d(V(iefkt),   efkt,   3*maxt2)
        call  cp_d2d(V(ixr),     xr,     2*maxt2 *maxint)
        call  cp_d2d(V(icoeff),  coeff,  6*maxt2 *maxint)
        endif   





C        
C   also these paramters needs to be set and are now hardwired in
C   derivs2.f:
c      :  Telec, scfhelp, scftol,qnel  
c      
!      call gettab(ntype)
!      write(*,*)'.... done with dftb_init'

c 
c   copy data  into buffer for  transfer
c      ntype,lmax
c      qzero,uhub, dr, dim, skself,espin, uhubbh,qzeroh,skhtab,
c      numint,cutoff,ppp, xr, efkt, coeff,
      
!      call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
!      write(*,*)'dftb_init_mpi, myid= ' , myid

      return

      end
c----------------------
      subroutine cp_d2d(A,B,N)
      real*8 A(N),B(N)
      integer i,N

      do i =1,N
          B(i) = A(i)
          enddo
      end

c----------------------
      subroutine cp_i2i(A,B,N)
      integer A(N),B(N)
      integer i,N

      do i =1,N
          B(i) = A(i)
          enddo
      end

c--------------------
c testwhether enough space allocated      
      subroutine  tstcor(iend,MDV)
      integer iend,mdv

      if (iend.gt.mdv)  then
         write(*,*) 'stop: not enough space in gettab'

         call  mpi_abort()
         stop
         endif

      end

c----------------------
      subroutine set_d2d(A,N)
      real*8  A(N)
      integer i,N

      do i =1,N
           A(i) = i 
          enddo
      end
      
c----------------------
      subroutine set_i2i(A,N)
      integer A(N)
      integer i,N

      do i =1,N
           A(i) = i 
          enddo
      end
      
c----------------------
      subroutine prn_d2d(A,N)
      integer A(N)
      integer i,N
           write(*,100)'->:',(A(i),i =1,N) 
100   format(' D:->',12(' ',E14.7))           
      end
      
c----------------------
      subroutine prn_i2i(A,N)
      integer A(N)
      integer i,N
           write(*,100)'->:',(A(i),i =1,N) 
100   format(' I:->',12(' ',I6))           
      end
      
