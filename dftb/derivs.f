c-----------------------------------------------------------
        subroutine   derivs(nat,iat,C,Energy,frc)        
        implicit real*8 (a-h,o-z)
        real*8  C(3,*), frc(3,*)
        integer iat(*)
        character*70 blabel 
        logical lopf
        blabel=' '
        ibohm=1
         
         inquire(101,opened=lopf)
         if (.NOT.LOPF) then 
           open(101,file='run.dftb',status='unknown')
         else   
           rewind(101)
           endif

         call pbdylcao(ibohm,nat,iat,C,Frc,Eelec,Emermin,blabel,V)
         energy =  Emermin 
!         stop
          return
          end
