# makefile
objects = qtes.o dumpef.o read_input_mpi.o set_traj.o ./dftb/libdftb.a
FC = mpif90
#--- Intel + MKL (original):
#switch = -mkl=sequential
#LIB = -openmp
#--- gfortran + OpenBLAS:
switch = -O2 -fallow-argument-mismatch -std=legacy
LIB = -fopenmp
BLASLIB = -lopenblas
EXEC = qtes.x

$(EXEC): $(objects)
	$(FC) -o $(EXEC) $(LIB) $(switch) $(objects) $(BLASLIB)
%.o: %.f
	$(FC) -c $(LIB) $(switch) $<

allclean:
	rm -rf *.o qtes.x *.log core.* fort.* data/*
