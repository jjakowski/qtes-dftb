!-------------------------------------------------------------
!  Old built-in DFTB engine (dftb/libdftb.a) on one geometry, set up
!  and called exactly as qtes.x does it (main program + derivs in
!  qtes.f), for comparison with tools/dftbplus/dftbp_driver.
!
!  Reads INdyn (line 1: nat), INxyz (6 parameter lines, then nat
!  atoms "sym x y z itype" in Angstrom) and IN_dftb + *.spl.
!  Prints energy and gradient at the INxyz geometry and with atom 1
!  moved by +0.1 bohr in z; writes old_result.dat.
!  "old_engine series FILE": instead, energy + gradient for every
!  geometry of a qt.dat file (positions in bohr on every 2nd line),
!  in file order; writes old_series_result.dat (index, energy, time)
!  and old_series_grad.dat (gradient per geometry).
!  Run on 1 MPI rank.
!-------------------------------------------------------------
        program old_engine
        implicit real*8 (a-h,o-z)
        include 'mpif.h'
        parameter(maxat=1000)
        real*8 c(3,maxat),c0(3,maxat),g(3,maxat),g0(3,maxat)
        integer iat(maxat)
        character*2 sym(maxat)
        character*80 line
        character*256 arg,sfile
        integer*8 ic0,ic1,irate
        include '../../dftb/mpiio.inc'

        call mpi_init(ierr)
!------- as in the qtes.f main program -----------------------
        open(istdin,file='IN_dftb',status='old')
        open(istdout,file='out.dftb',status='unknown')
        call dftb_init_mpi(ntype)
!------- geometry, converted as in read_input_mpi.f ----------
        open(9,file='INdyn',status='old')
        read(9,*) nat
        close(9)
        open(8,file='INxyz',status='old')
        do i=1,6
         read(8,'(a)') line
        enddo
        i=0
10      read(8,'(a)') line
        if(len_trim(line).eq.0) goto 10
        i=i+1
        read(line,*) sym(i),(c0(k,i),k=1,3),iat(i)
        if(i.lt.nat) goto 10
        close(8)
        do i=1,nat
         do k=1,3
          c0(k,i)=c0(k,i)/0.52918d0
         enddo
        enddo
!------- series mode -----------------------------------------
        arg=' '
        if(command_argument_count().ge.1)
     $     call get_command_argument(1,arg)
        if(trim(arg).eq.'series') then
         call get_command_argument(2,sfile)
         call system_clock(count_rate=irate)
         open(30,file=trim(sfile),status='old')
         open(31,file='old_series_result.dat',status='unknown')
         open(32,file='old_series_grad.dat',status='unknown')
         write(31,'(a)') '# index energy(Ha) time(s)'
         n=0
20       read(30,*,end=30) ((c(k,i),k=1,3),i=1,nat),w
         read(30,'(a1)') line
         n=n+1
         call system_clock(ic0)
         call oldcalc(nat,iat,c,e,g)
         call system_clock(ic1)
         write(31,'(i6,es24.16,es12.4)') n,e,dble(ic1-ic0)/dble(irate)
         write(32,'(1000es24.16)') ((g(k,i),k=1,3),i=1,nat)
         goto 20
30       close(30)
         close(31)
         close(32)
         write(*,*) 'series: ',n,' geometries'
         call mpi_finalize(ierr)
         stop
        endif
!------- 1. INxyz geometry -----------------------------------
        call oldcalc(nat,iat,c0,e0,g0)
        write(*,*)
        write(*,*) '== INxyz geometry'
        call report(nat,sym,e0,g0)
        open(20,file='old_result.dat',status='unknown')
        write(20,'(a)') '# Mermin free energy (Ha), then gradient'//
     $                  ' (Ha/bohr) per atom'
        write(20,'(es24.16)') e0
        do i=1,nat
         write(20,'(a2,3es24.16)') sym(i),(g0(k,i),k=1,3)
        enddo
        close(20)
!------- 2. atom 1 moved by +0.1 bohr in z -------------------
        do i=1,nat
         do k=1,3
          c(k,i)=c0(k,i)
         enddo
        enddo
        c(3,1)=c(3,1)+0.1d0
        call oldcalc(nat,iat,c,e,g)
        write(*,*)
        write(*,*) '== atom 1 moved by +0.1 bohr in z'
        call report(nat,sym,e,g)
        write(*,'(a,es14.6,a)') ' Energy change: ',e-e0,' Ha'
        call mpi_finalize(ierr)
        end
!-------------------------------------------------------------
        subroutine oldcalc(nat,iat,c,e,g)
!------- the body of derivs (qtes.f): same settings and calls --
        implicit real*8 (a-h,o-z)
        real*8 c(3,*),g(3,*),v(1)
        integer iat(*),iv(1)
        character*70 blabel
        logical scfhelp
        blabel=' '
        ibohm=1
        ntype=2
        qnel=0.0
        scfhelp=.true.
        scftol=1e-6
        telec=6000.0
        call do_dftb(ibohm,nat,iat,c,g,eelec,emermin,blabel,
     $       ntype,qnel,scfhelp,scftol,telec,v,iv)
        e=emermin
!------- do_dftb returns forces; derivs returns the gradient --
        do i=1,nat
         do k=1,3
          g(k,i)=-g(k,i)
         enddo
        enddo
        end
!-------------------------------------------------------------
        subroutine report(nat,sym,e,g)
        implicit real*8 (a-h,o-z)
        real*8 g(3,*)
        character*2 sym(*)
        write(*,'(a,f20.12,a)') ' Mermin free energy: ',e,' Ha'
        write(*,'(a)') ' Gradient (Ha/bohr):'
        do i=1,nat
         write(*,'(i4,1x,a2,3f18.12)') i,sym(i),(g(k,i),k=1,3)
        enddo
        end
