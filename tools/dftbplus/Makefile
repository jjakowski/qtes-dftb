# Milestone-1 driver linked against the installed DFTB+ library.
# DFTB+ is found through pkg-config; see BUILD.md for building it.
# Use the rmg-dev environment (same gfortran as the DFTB+ build).

DFTBPLUS_PREFIX ?= /home/j2c/Development/DFTB/dftbplus-25.1/install-qtes
PKGCONF  = PKG_CONFIG_PATH=$(DFTBPLUS_PREFIX)/lib/pkgconfig pkg-config

FC      = gfortran
FFLAGS  = -O2 -g $(shell $(PKGCONF) --cflags dftbplus)
LIBS    = $(shell $(PKGCONF) --libs --static dftbplus)

dftbp_driver: dftbp_driver.f90
	$(FC) $(FFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f dftbp_driver *.o *.mod

.PHONY: clean
