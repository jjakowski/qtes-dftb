!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine set_traj(nd,iat,np,pow,x,p,s,w,u,du,v,dv,
     $             myid,numprocs)
!------- read input file: ndx, ndq, npx are from bohm.inc ----
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        include 'mpif.h'
        real*8 vlocal(npx), dvlocal(ndx,npx) 
        character*2 asymbol(ndx)
        integer iat(ndx),idum(ndx)
        real*8 x(ndx,npx),p(ndx,npx),s(npx),w(npx),dx(ndx)
        real*8 v(npx),dv(ndx,npx),u(npx),du(nqx,npx),q(ndx),pq(ndx)
!  dver(0)       : ehrenfest averaged potential over the wavefunction
!  dver(1:ndx)   : ehnrefest averaged forces    over the wavefunction
        real*8 dver(0:ndx)
        real*8 ta,tb
        common/wf/a0(nqx,nqx),q0(ndx),p0(ndx),isr
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx),
     &          ipot,ipf,nqd,kmax,kout,krst
        common/molecule/asymbol,nat
        common/pot_xyz/ qc(3,ndx),ql,qls,dw,sc,esh,dch,zh,ch,bnd

        if(myid.eq.0) then

        write(*,*) 'READING TRAJ INFO via set_traj'

         call rdsav(nd,np,tc,x,p,s,w)

        endif  !<---- myid=0 is done with reading

!---------------------------------------------------------------
        call rankchk(np,numprocs,myid,1)
c distribute work
        nwork = (np + numprocs-1) /numprocs
        ip0 = myid*nwork  + 1 
        ip1 = myid*nwork  + nwork
        if(ip1.gt.np)  ip1 = np     ! the last one has the least work to do  
        
!---------------- compute classical potential ------------------
        nat1 = nd/3
         call mysecond(ta)
        do 112 ip =ip0,ip1
         !call derivs0(nd,iat,x(1,i),vlocal(i),dvlocal(1,i))
         !call derivs0(nd,iat,x(1,ip),vlocal(ip),dvlocal(1,ip))       !<--- old(fit)
         call derivs(nat1,iat,x(1,ip),vlocal(ip),dvlocal(1,ip))       !<--- new(dftb)
         vlocal(ip) = vlocal(ip) -esh                                 !<-- this probably should  be moved  to myid=0
         
112     continue
         call mysecond(tb)
        if (myid.eq.0) write(22,*) ta,tb,tb-ta
        call flush(22)
c
c  now collect all potentials and from  all Mpi tasks         
        nlocal = ip1-ip0+1
        call mpi_gather(vlocal(ip0),nlocal,MPI_DOUBLE_PRECISION,
     $        v,nlocal,MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)   

        call mpi_gather(dvlocal(1,ip0),ndx*nlocal,MPI_DOUBLE_PRECISION,
     $        dv,ndx*nlocal,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

        call mpi_barrier(MPI_COMM_WORLD,ierr)

        call mpi_Bcast(v,np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        
        call mpi_Bcast(dv,ndx*np,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $        ierr)
!--------------------------------------------------------
!--------------------- Ehrenfest averaging <<<<<<<<<<<<<<
!--------------------- added 09/24/2012 SG --------------
        if(ipf.eq.0) then
         do j=0,nd
          dver(j)=0d0
         enddo
         do i=1,np
          dver(0)=dver(0)+w(i)
          do j=nqd+1,nd
           dver(j)=dver(j)+dv(j,i)*w(i)
          enddo
         enddo
         do i=1,np
          do j=nqd+1,nd
           dv(j,i)=dver(j)/dver(0)
          enddo
         enddo
        endif
!------------------------->>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!-------- s(0) is absorbed into the sampling --------------------
        do 100 i=1,np
         ek=0d0
         do 108 n=1,nd
108      ek=ek+p(n,i)**2/2d0/am(n)
         epot=epot+v(i)*w(i)
         ekin=ekin+ek*w(i)
         anrm=anrm+w(i)
100     continue
         
        t=0d0
        call qpot(np,nd,t,w,x,s,u,du,qp,myid,numprocs)
        if (myid.eq.0) then 
          write(*,*) 'NORM=',anrm
          write(*,*) 'KIN energy = ',ekin/anrm
          write(*,*) 'POT energy = ',epot/anrm
          write(*,*) 'QP energy = ',qp/anrm
          write(*,*) 'TOTAL energy = ',(qp+epot)/anrm
          endif
        ener=(qp+epot)/anrm
!---------------------------------------------------------------------------
!----------- examine the initial point distribution ------------------------
        npr=nd
        if(npr.gt.40) npr=40
        if (myid.eq.0) then 
          open(28,file='wf0')
          do i=1,np
            write(28,1000) (x(j,i),j=1,npr)
          enddo
          close(28)
          endif 
!        write(*,*) 'LEAVE INITIALIZE'
!--------------------------------------------------------------------
1000    format(200(e12.5,1x))
1007    format(700(e14.7,1x))
        return
        end
!-------------------------------------------------------------
!-------------------------------------------------------------
        subroutine rdsav(nd,np,tc,x,p,s,w)
!------- read the restart file traj.sav on rank 0 (issue #9) -----
!------- optional header: '# np nd dt kout', checked against the --
!------- inputs; then the time and 2 lines per QT. Every QT must --
!------- be complete and the file must hold exactly np QTs. -------
        implicit real*8 (a-h,o-z)
        include 'bohm.inc'
        include 'mpif.h'
        real*8 x(ndx,npx),p(ndx,npx),s(npx),w(npx)
        character*256 line
        character*1 c1
        common/prop/dt,dt2,tau(ndx),tau2(ndx),am(ndx),
     &          ipot,ipf,nqd,kmax,kout,krst
        open(32,file='traj.sav',status='old',iostat=ios)
        if(ios.ne.0) then
         write(*,*) 'RESTART ERROR: cannot open traj.sav'
         call savfail
        endif
        read(32,'(a)',iostat=ios) line
        if(ios.ne.0) then
         write(*,*) 'RESTART ERROR: traj.sav is empty'
         call savfail
        endif
!------- header line written since 2026-11: '# np nd dt kout' -----
        ihd=0
        line=adjustl(line)
        if(line(1:1).eq.'#') then
         ihd=1
         read(line(2:),*,iostat=ios) nps,nds,dts,kos
         if(ios.ne.0) then
          write(*,*) 'RESTART ERROR: cannot read the traj.sav header'
          call savfail
         endif
         nerr=0
         if(nps.ne.np) then
          write(*,1001) nps,np,nps
          nerr=nerr+1
         endif
         if(nds.ne.nd) then
          write(*,1002) nds,nds/3,nd,nd/3
          nerr=nerr+1
         endif
         if(abs(dts-dt*kout).gt.1d-10*abs(dt*kout)) then
          write(*,1003) dts,dt*kout
          nerr=nerr+1
         endif
         if(kos.ne.kout) then
          write(*,1004) kos,kout
          nerr=nerr+1
         endif
         if(nerr.gt.0) call savfail
         read(32,*,iostat=ios) tc
        else
         read(line,*,iostat=ios) tc
        endif
        if(ios.ne.0) then
         write(*,*) 'RESTART ERROR: cannot read the time in traj.sav'
         call savfail
        endif
!------- the np QTs; a short or cut line ends in an error/EOF -----
        do i=1,np
         read(32,*,iostat=ios) (x(k,i),k=1,nd),w(i)
         if(ios.lt.0) then
          write(*,1005) i-1,np
          if(ihd.eq.0) write(*,1008) i-1
          if(ihd.eq.1) write(*,1009) nps
          call savfail
         endif
         if(ios.eq.0) read(32,*,iostat=ios) (p(k,i),k=1,nd),s(i)
         if(ios.ne.0) then
          write(*,1006) i
          call savfail
         endif
        enddo
!------- anything left means more QTs than np ---------------------
        nl=0
        do
         read(32,'(a1)',iostat=ios) c1
         if(ios.ne.0) exit
         nl=nl+1
        enddo
        if(nl.gt.0) then
         write(*,1007) np+nl/2,np
         if(ihd.eq.0) write(*,1008) np+nl/2
         if(ihd.eq.1) write(*,1009) nps
         call savfail
        endif
        close(32)
1001    format(' RESTART ERROR: np changed between runs: traj.sav ',
     $   i0,', INdyn line 9 ',i0/'   Set np in INdyn to ',i0,
     $   ' to continue that run.')
1002    format(' RESTART ERROR: number of coordinates changed: ',
     $   'traj.sav ',i0,' (',i0,' atoms), inputs ',i0,' (',i0,
     $   ' atoms)')
1003    format(' RESTART ERROR: time step changed: traj.sav dt=',
     $   g0,', INdyn line 5 dt=',g0)
1004    format(' RESTART ERROR: kout changed: traj.sav ',i0,
     $   ', INdyn line 4 ',i0)
1005    format(' RESTART ERROR: traj.sav holds ',i0,
     $   ' QTs, INdyn line 9 asks for ',i0)
1006    format(' RESTART ERROR: traj.sav ends inside QT ',i0,
     $   ' or cannot be read there (incomplete file)')
1007    format(' RESTART ERROR: traj.sav holds ',i0,
     $   ' QTs, INdyn line 9 asks for ',i0)
1008    format('   (traj.sav has no header.) If np has changed since',
     $   ' the run that'/'   wrote it, set np in INdyn to ',i0,
     $   '; otherwise traj.sav is'/'   incomplete or from another run.')
1009    format('   traj.sav is damaged: its header says ',i0,' QTs.')
        return
        end
!-------------------------------------------------------------
        subroutine savfail
!------- stop all ranks after a restart-file error (rank 0) -------
        include 'mpif.h'
        call flush(6)
        call mpi_abort(MPI_COMM_WORLD,1,ierr)
        end
