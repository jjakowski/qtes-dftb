C Output from Public domain Ratfor, version 1.0
      subroutine gettab_mpi(ntype)
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
      real*8 skhtab(10,maxtab,maxtyp,maxtyp) 
      real*8 skstab(10,maxtab,maxtyp,maxtyp)
      real*8 espin(maxtyp)
      real*8 skself(3,maxtyp),dr(maxtyp,maxtyp),qzeroh(3),uhubbh(3)
      real*8 coeff(6,maxint,maxtyp,maxtyp),xr(2,maxint,maxtyp,maxtyp)
      real*8 efkt(3,maxtyp,maxtyp),cutoff(maxtyp,maxtyp)
      integer numint(maxtyp,maxtyp)
      integer dim(maxtyp,maxtyp),ppp,lmax(maxtyp)
      character*64 skfile
      character*128 chdummy
      common /sktab/ skhtab,skstab,skself,dr,dim
      common /spin/ espin
      common /spltab/ coeff,xr,efkt,cutoff,numint
      common /mcharge/ qzero(maxtyp), uhubb(maxtyp)
      common /lmax/ lmax

c   parallel  trajectories execution 
      include 'mpif.h'
      include 'mpiio.inc'     !ninfile, noutfile
      
      DO i = 1,ntype 
        qzero(i) = 0.0d0
        uhubb(i) = 0.0d0

        Do j = 1,ntype 
          write(istdout,*)'enter file name for sk-data of pair ',i,j
          call flush(istdout)
          read(istdin,*)skfile
          open (3,file=skfile,status='unknown')
          rewind 3

          if(i .ne. j)then
            read (3,*) dr(i,j),dim(i,j)
          endif

          if(i .eq. j)then
C     read (3,*) dr(i,j),dim(i,j),lmax(i)
           read (3,*) dr(i,j),dim(i,j)
           read (3,*) (skself(l,i),l = 1,3),espin(i),
     &                (uhubbh(4-l), l=1,3), (qzeroh(4-l), l=1,3)
           uhubb(i) = uhubbh(1)

           if(uhubb(i).eq.0.4195007d0)then
           write(istdout,*)'+++modified version of gamma used of X-H+++'
           endif

           do k = 1,3 
             qzero(i) = qzero(i) + qzeroh(k)
           enddo

           do l=1,3 
             espin(i) = espin(i) + skself(l,i)*qzeroh(4-l)
           enddo
          endif

          do k = 1,dim(i,j) 
            read (3,*) (skhtab(l,k,i,j),l = 1,10), 
     &                 (skstab(l,k,i,j),l = 1,10)
          enddo

23016     if(.true.)then
            read(3,'(A)') chdummy
            if(chdummy.eq.'Spline')then
              goto 23017
            endif
            goto 23016
          endif
23017     continue

          read(3,*) numint(i,j),cutoff(i,j)

          if(numint(i,j).gt.maxint)then
            write(istdout,*)'Too many intervalls!'
            goto 99
          endif

          read(3,*) (efkt(ppp,i,j),ppp=1,3)

          do ppp=1,numint(i,j) 
            if(ppp.lt.numint(i,j))then
               read (3,*) xr(1,ppp,i,j),xr(2,ppp,i,j),
     &                    coeff(1,ppp,i,j),coeff(2,ppp,i,j),
     &                    coeff(3,ppp,i,j),coeff(4,ppp,i,j)
            else
               read (3,*) xr(1,ppp,i,j),xr(2,ppp,i,j),
     &                    coeff(1,ppp,i,j),coeff(2,ppp,i,j),
     &                    coeff(3,ppp,i,j),coeff(4,ppp,i,j),
     &                    coeff(5,ppp,i,j),coeff(6,ppp,i,j)
            endif
          enddo

          if(xr(2,numint(i,j),i,j).ne.cutoff(i,j))then
            write(istdout,*) 'Error in Datafile'
            goto 99
          endif

          write(istdout,*) 'skfile for pair ',i,j,' :',skfile
          close (3)

        ENDDo

      ENDDO
99    continue

c prepare buffer  for before  MPI transfer 

      maxt2 = maxtyp*maxtyp     ! currently: 36=6*6

c copy everything 
      intype   =  1                       ! Int
      ilmax    = intype    +     1        ! int
      idimm    = ilmax     +    maxtyp    ! Int
      inumint  = idimm     +    maxt2     ! Int
      idr      = inumint   +    maxt2
      iskself  = idr       +    maxt2 
      iespin   = iskself   +  3*maxtyp
      iuhubb   = iespin    +    maxtyp 
      iqzero   = ihubb     +    maxtyp
      iskhtab  = iqzero    +    maxtyp
      iskstab  = iskhtab   + 10*maxt2 *maxtab
! now splines
      icutoff  = iskstab   + 10*maxt2 *maxtab 
      iefkt    = icutoff   +    maxt2 
      ixr      = ieftk     +  3*maxt2
      icoeff   = ixr       +  2*maxt2 *maxint
      iend     = icoeff    +  6*maxt2 *maxint

      Mints    = idr  - intype                    ! how many  integers we store ??
      Mfloats  = iend - idr                       ! how many  floats

c  check if enough space to copy      
      call  tstcor(iend,mdv)

c now copy everything  to buffer
      call  cp_i2i(ntype,  V(intype),  1)
      call  cp_i2i(lmax,   V(ilmax),   maxtyp)
      call  cp_i2i(dim,    V(idimm),   maxtyp)
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


       
      call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
      call MPI_Bcast(V(1),iend,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     $               ierr)
!          call MPI_BCAST(outfile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
!          call MPI_BCAST(ntraj,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)



      
  
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
         stop
         endif

      end    
