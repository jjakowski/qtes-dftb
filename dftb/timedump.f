
       subroutine timedump()
       implicit real *8 (a-h,o-z)
       Character*(70) Date1
       include 'mpiio2.inc'

       
        write(date1,*) 'Current time:       '
        call  time(Date1(21:30))
        Date1(31:34) = ' on '
        Call Date(Date1(35:44))
 
        write(istdout,*)date1


        end

