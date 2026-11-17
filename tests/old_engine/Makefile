# Old built-in DFTB engine on one geometry (milestone 1, step 3, check c).
# Needs dftb/libdftb.a (cd ../../dftb && make libdftb). Use rmg-dev.

FC     = mpif90
FFLAGS = -O2 -fopenmp -fallow-argument-mismatch -std=legacy
LIBS   = ../../dftb/libdftb.a -lopenblas

old_engine: old_engine.f ../../dftb/libdftb.a
	$(FC) $(FFLAGS) -o $@ old_engine.f $(LIBS)

clean:
	rm -f old_engine *.o

.PHONY: clean
