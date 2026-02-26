c-----------------------------------------------------------
        subroutine   derivs2(nat,iat,C,Energy,frc)        
        implicit real*8 (a-h,o-z)
        real*8  C(3,*), frc(3,*)
        integer iat(*),lmax(100)
        character*70 blabel 
        logical lopf,scfhelp
        blabel=' '
        ibohm=1
         
         inquire(101,opened=lopf)
         if (.NOT.LOPF) then 
           open(101,file='run.dftb',status='unknown')
         else   
!           rewind(101)
           endif
           
        ntype   = 2     ! max atoms type
!  lmax passed  from dftb_init via common         
!        lmax(1) = 2     ! carbon  
!        lmax(2) = 1     ! hydrogen
        qnel = 0.0 
        scfhelp = .true.
        scftol = 1e-6
!        telec = 1000.0
        telec = 6000.0



           call do_dftb(ibohm,nat,iat,C,Frc,Eelec,Emermin,blabel,
     $         ntype, qnel,scfhelp, scftol,telec, V)
         energy =  Emermin 


!         stop
          return
          end
