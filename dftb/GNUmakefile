TARGET=dftb
MAXIMA=maxima
.SUFFIXES :
#.SUFFIXES : .o .r .f .inc .h
#.SUFFIXES : .o .f .inc
.SUFFIXES : .o .f .F .inc
#
include GNUmake.inc
#
#.h.inc : 
#	@$(RATFOR) < $*.h > $*.inc
#	@printf .
#.r.o :
#	@$(RATFOR) < $*.r > $*.f
#	@$(FC) $(FOPT) -c $*.f
#	@printf .
.f.o :
	@$(FC) $(FOPT) -c $*.f
	@printf .
#
#RSRC= gettab.f skpar.f my_mulliken.f \
#      slkode.f vnorm.f xnorm.f geometries.f repulsiv.f fermi.f \
#      my_gradient.f 
#
INC=  commonbox.inc
OBJ0=  pdftb.o eglcao.o skpar.o gettab.o xnorm.o geometries.o \
      vnorm.o slkode.o ewevge.o slktrafo.o self.o input.o \
      cgrad.o ranmar.o broyden.o my_mulliken.o my_gradient.o \
      short_range.o long_range.o shift.o repulsiv.o  \
      fermi.o output.o gamma.o gammamat.o atomen.o \
      dispersion_egr.o dispersionread.o  externalchgrad.o \
      externalshift.o  reada.o  timedump.o
OBJ1=  eglcao.o skpar.o gettab.o xnorm.o geometries.o \
      vnorm.o slkode.o ewevge.o slktrafo.o self.o input.o \
      cgrad.o ranmar.o broyden.o my_mulliken.o my_gradient.o \
      short_range.o long_range.o shift.o repulsiv.o  \
      fermi.o output.o gamma.o gammamat.o atomen.o \
      dispersion_egr.o dispersionread.o  externalchgrad.o \
      externalshift.o  reada.o   dumpm.o  \
      dftb_init.o  dftb_init_mpi.o do_dftb.o  derivs2.o  

#-- old:
#OBJ= $(OBJ1)  potext.o 
#-- new:
#OBJ= $(OBJ1)  potext2.o 
#--new+mpi:
OBJ= $(OBJ1)  potext2_mpi.o 

OBJLIB= $(OBJ1)
#
all: $(TARGET) libdftb
#
pdftb.o:  $(MAXIMA).inc
atomen.o: $(MAXIMA).inc
broyden.o: $(MAXIMA).inc
cgrad.o: $(MAXIMA).inc
dispersion_egr.o: $(MAXIMA).inc
dispersion_read.o: $(MAXIMA).inc
pdylcao.o: $(MAXIMA).inc
eglcao.o: $(MAXIMA).inc
ewevge.o: $(MAXIMA).inc
externalchgrad.o: $(MAXIMA).inc
externalshift.o: $(MAXIMA).inc
fermi.o: $(MAXIMA).inc
gammamat.o: $(MAXIMA).inc
gettab.o: $(MAXIMA).inc
input.o: $(MAXIMA).inc
my_gradient.o: $(MAXIMA).inc
my_mulliken.o: $(MAXIMA).inc
output.o: $(MAXIMA).inc
repulsiv.o: $(MAXIMA).inc
shift.o: $(MAXIMA).inc
skpar.o: $(MAXIMA).inc
slkode.o: $(MAXIMA).inc
derivs.o:  $(MAXIMA).inc	
#
#

$(TARGET): $(MAXIMA).inc $(INC) $(OBJ) 
	@$(FC) $(FOPT) -o $(TARGET) $(OBJ) $(LAPACK) $(BLAS)
	@strip dftb
	@echo 
	@echo "........  Make of $(TARGET) for $(PLATFORM) completed."
	@echo "........  Rquired datasize limit to run it: `size dftb | \
         awk '/data/{getline; a=$$4; print int((a/1024)/1024)+1}'`M"
#	@rm -f $(RSRC) $(INC) $(OBJ) $(MAXIMA).h $(MAXIMA).inc
libdftb:  $(MAXIMA).inc $(INC) $(OBJLIB)
	 @$(FC) $(FOPT) -c $(OBJLIB)
	 ar -r  libdftb.a  $(OBJLIB)
	
clean:
	@rm -f $(TARGET) $(RSRC) $(OBJ) 
	@echo "     ** "
	@echo "     ** $(TARGET) deleted..."
	@echo "     ** "
allclean: clean
	@rm -f dsygvd.a
	@echo "     ** "
	@echo "     **                        all clear..."
	@echo "     ** "
