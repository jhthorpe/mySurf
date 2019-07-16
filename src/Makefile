#local
export FC = gfortran
#export FFLAGS = -O3 
export FFLAGS = -g -fcheck=bounds -Wall -Wextra -Wshadow
export bindir = ../bin

#local
#export LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas

#HiPerGator
export MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
export LINAL = -L$(MKL) -Wl,-R$(MKL) -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm


objects = xsurf.o input.o MD.o V.o nco.o linal.o

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp xsurf.o $(bindir)/xsurf

input.o : input.f90
	$(FC) $(FFALGS) -c input.f90

linal.o : linal.f90
	$(FC) $(FFLAGS) $(LINAL) -c linal.f90 

nco.o : nco.f90 linal.o
	$(FC) $(FFLAGS) -c nco.f90 linal.o

V.o : V.f90
	$(FC) $(FFLAGS) -c V.f90

MD.o : MD.f90 V.o input.o nco.o 
	$(FC) $(FFLAGS) -c MD.f90 V.o input.o nco.o 

xsurf.o : surf.f90 input.o MD.o nco.o linal.o V.o
	$(FC) $(FFLAGS) $(LINAL)  -o xsurf.o surf.f90 input.o MD.o nco.o linal.o V.o

clean:
	rm *.o *.mod


