#local
export FC = gfortran
#export FFLAGS = -O3 
export FFLAGS = -g -fcheck=bounds -Wall -Wextra -Wshadow
export bindir = ../bin

objects = xsurf.o input.o MD.o V.o 

all: $(objects)
	if [ ! -d $(bindir) ]; then mkdir -p $(bindir); fi
	cp xsurf.o $(bindir)/xsurf

input.o : input.f90
	$(FC) $(FFALGS) -c input.f90

V.o : V.f90
	$(FC) $(FFLAGS) -c V.f90

MD.o : MD.f90 V.o
	$(FC) $(FFLAGS) -c MD.f90 V.o

xsurf.o : surf.f90 input.o MD.o V.o
	$(FC) $(FFLAGS) -o xsurf.o surf.f90 input.o MD.o V.o 

clean:
	rm *.o *.mod


