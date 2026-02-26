        program potext
        implicit real*8 (a-h,o-z)
        include 'maxima.inc'
        real*8 C(3,nndim),Frc(3*nndim)
        integer IAt(nndim)
        !real*8 V(mdim*mdim)
        parameter (MDV=10000000)
        real*8 V(MDV)
        character infile*80, outfile*80
        character*70 blabel 
       

       include 'mpif.h'


       integer ierr, myproc, myid
       integer clock0,clock1, nclockrate
       character dum*80,line*255,name1*1,name2*2, name3*3,name4*4
       logical ex

c
c -parallel trajectories
c
       include 'mpiio.inc'    ! defines filehandlers:  istdio, istdout
       infile='run-new2.dftb'
       outfile='out-new2.dftb'


c
c  initialize MPI info
c
      call mpi_init(ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )



c
c  go to the working  directory     
c      
      !write(*,*)' MPI_COMM_WORLD =', MPI_COMM_WORLD
      !write(*,*) 'numprocs =', numprocs
11     format(I1)
12     format(I2)
13     format(I3)
14     format(I4)

c   working directory for  'myid'    (it is the same as myid)      
       if (myid.le.9) then
             write(name1(:1),11)myid
             call chdir(name1)
       elseif  (myid.ge.10.AND.myid.le.99)  then
             write(name2(:2),12)myid
             call chdir(name2)
       elseif  (myid.ge.100.AND.myid.le.999)  then
             write(name3(:3),13)myid
             call chdir(name3)
       elseif  (myid.ge.1000.AND.myid.le.9999) then
             write(name4(:),14)myid
             call chdir(name4)
       else
            call  MPI_ABORT()
       endif

       !call system("pwd")
       call getcwd(line)               ! get working directory
       write(*,*)'task:',myid,'--> ',trim(line)



        blabel = 'dupa:    '
        ii =0
! initialize dftb ,how many different atoms (see run-new.dftb read from dftb_init.f and gettab )
!        write(*,*) 'ntype = ?'
!        read(*,*) ntype
!        write(*,*) ntype


         !inquire(file="CHR.DAT",exist=ex)
         inquire(file=infile,exist=ex)
         if(ex)then
             open(istdin,file=infile,status='old')
           else
             write(*,*)myid,' cannot find file', infile
             write(*,*)' stop!'
             call mpi_abort()
           endif 
            
         open(istdout,file=outfile,status='unknown')

      if (numprocs.eq.1) then 
        call dftb_init(ntype,ibohm,nat,iat,C,Frc,Eelec,Emermin,blabel,V)
      else
        call dftb_init_mpi(ntype,ibohm,nat,iat,C,Frc,Eelec,Emermin,
     $          blabel,V,MDV)
        !----- temp -----------------------------------
        ! ignore first few lines in 'infile' if mpi task >0 to adjust offset
        if (myid.ne.0)  then        ! adjust offset in infile since alread set in  dftb_init_mpi
           read(istdin,*)           ! ntype
           read(istdin,*)           ! lmax
           do i=1,ntype*ntype       ! loop over  sklfiles 
              read(istdin,*)        ! names of skfiles   
              enddo
           endif     
        !----end of temp ------------------------------
      endif


100        write(istdout,*)'Nat = ?'
        !read(*,*) Nat
        read(istdin,*) Nat
        write(istdout,*)'Nat =',Nat 
        !read(*,*)
        read(istdin,*)
        do i =1,Nat
          !read(*,*)ii, iiat, x,y,z
          read(istdin,*)ii, iiat, x,y,z
          IAt(i) = iiat
          C(1,i) = x  
          C(2,i) = y  
          C(3,i) = z  
          enddo

           iii= iii+1
           write(blabel,101) iii
101       format('bohm: ', I6)
          write(istdout,*)'***** C :'
          do i=1,Nat 
             write(istdout,110)iat(i),C(1,i),C(2,i),C(3,i)
             enddo

          ibohm=1
cconvert  CC from Angstroms  to bohrs 
          bohr2A =  0.52917724
          A2bohr =1.0/bohr2A
          do i=1,Nat
             C(1,i) = C(1,i)*A2bohr
             C(2,i) = C(2,i)*A2bohr
             C(3,i) = C(3,i)*A2bohr
             enddo

c--          
!          call timedump()
!          open(istdin,file=infile,status='old')
!          open(istdout,file=outfile,status='unknown')
!          call pbdylcao(ibohm,nat,iat,C,Frc,Energy,Emermin,blabel,V)
          call derivs2(nat,iat,C,Energy,frc) 
!          call timedump()

          write(istdout,*)'------------done-----------'
          write(istdout,*)' Energy       = ', Energy

          write(istdout,*)'*** Forces'
          write(istdout,111) (Frc(i),i=1,3*nat)

          write(istdout,*)'============================================'
          write(istdout,*)' more ?  /1=yes/'
          read(istdin,*) more
          if (more.eq.1)  goto 100

110       format(I6,100(' ',E22.16))          
111       format(3('   ',F22.8))          

          close(istdin)
          close(istdout)

          call  mpi_finalize(ierr)

          end
c-----------------------------------------------------------
!       subroutine   derivs(nat,iat,C,Energy,frc)        
!       implicit real*8 (a-h,o-z)
!       real*8  C(3,*), frc(3,*)
!       integer iat(*)
!       character*70 blabel 
!       blabel=' '
!       ibohm=1
!        
!        call pbdylcao(ibohm,nat,iat,C,Frc,Eelec,Emermin,blabel,V)
!        energy =  Emermin 
!         return
!         end
