PLATFORM=$(shell uname)
ifeq ($(PLATFORM),OSF1)
	RATFOR = ratfor
	FC = f90
	FOPT = -O5 -tune ev6
	LAPACK =
#	LAPACK = dsygvd.a
	BLAS = -ldxml
endif
ifeq ($(PLATFORM),AIX)
        RATFOR = ratfor
        FC = f77 
        FOPT = -O3 -qstrict  
        LAPACK = dsygvd.a
        BLAS = -lblas 
endif
ifeq ($(PLATFORM),HP-UX)
        RATFOR = ratfor
        FC = f77
        FOPT = -O3  +Olibcalls +Odataprefetch +Ofastaccess +Oentrysched \
          +Oregionsched +Onolimit  +Onoloop_transform 
#	LAPACK =
        LAPACK = dsygvd.a
        BLAS = +Ofastaccess -Wl,+FPD \
          -Wl,-x -Wl,-N -Wl,-n -lveclib -llapack 
endif
ifeq ($(PLATFORM),IRIX64)
        RATFOR = ratfor
        FC = f77
        FOPT = -c -O3 -mips4 -r10000 
        LAPACK = dsygvd.a
        BLAS = -lcomplib.sgimath 
endif
ifeq ($(PLATFORM),Linux)
#	FC =  pgf90  -Minfo=all  -D_IFMPI
#	FC =  pgf90  -Minfo=all 
#---pgi:
#	FC =  mpif90  -Minfo=all -Mmpi
#	FOPT = -mp -O2  -pgf90libs  
	#FOPT = -mp -fastsse -O2  -pgf90libs  
#	FOPT = -mp -O2  -pgf90libs  
#	LAPACK =  -llapack -pgf90libs
#--- gpu:
	FC =  mpif90  
	FOPT =  -O2  -fopenmp
	LAPACK = -mkl=sequential
        BLAS = -mkl=sequential
#	LAPACK =  -llapack 
#	BLAS = -lblas  -g

endif

