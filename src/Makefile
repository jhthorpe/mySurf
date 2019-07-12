#local
export FC = gfortran
#export FFLAGS = -O3 
export FFLAGS = -g -fcheck=bounds -Wall -Wextra -Wshadow
export bindir = ../bin

objects = xsurf.o input.o MD.o #potential.o 

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp xsurf.o $(bindir)/xsurf

input.o : input.f90
	$(FC) $(FFALGS) -c input.f90

MD.o : MD.f90
	$(FC) $(FFLAGS) -c MD.f90

potential.o : potential.f90
	$(FC) $(FFLAGS) -c potential.f90

xsurf.o : surf.f90 input.o MD.o 
	$(FC) $(FFLAGS) -o xsurf.o surf.f90 input.o MD.o 

clean:
	rm *.o *.mod


