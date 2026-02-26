        subroutine dumpm(A,Nx,Ny,Mx,My)
        implicit real*8 (a-h,o-z)

        real*8 A(Mx,My)
!        real*8 A(Mx,*)


        do i=1,Nx
          write(*,100)(A(i,j),j=1,Ny)
          enddo
100     format(100(' ',E22.14))
        end

c-------------------------------------------------------------------
        subroutine dumpms(A,Nx,Ny,Mx,My)
        implicit real*8 (a-h,o-z)

        real*8 A(Mx,My)
!        real*8 A(Mx,*)

        write(*,*)'Nx,Ny,Mx,My= ', Nx,Ny,Mx,My
        do i=1,Nx
          write(*,100)(A(i,j),j=1,i),(A(j,i),j=i+1,Ny)
          enddo
100     format(100(' ',E22.14))
        end

c-------------------------------------------------------------------
        subroutine dumpmd(A,Nx,Ny,Mx,My)
        implicit real*8 (a-h,o-z)

        real*8 A(Mx,My)
!        real*8 A(Mx,*)

        write(*,100)(A(i,i),i=1,Nx)
100     format(100(' ',E22.14))
        end

