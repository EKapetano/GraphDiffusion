#Compiler and Linker - Options
#Used Compiler:
FC = gfortran
#Compiler options:
FCOPTS = -std=f2003 -pedantic -Wall -O3
#Used Linker:
LN = $(FC)
#Linker options:
LNOPTS = 

OBJS = accuracy.o inout.o crwalk.o

crwalk: $(OBJS)
	$(FC) $(FCOPTS) -fopenmp -o crwalk $(OBJS)


crwalk.o: crwalk.f90
	$(FC) $(FCOPTS) -fopenmp -c crwalk.f90 

accuracy.o: accuracy.f90
	$(FC) $(FCOPTS) -c accuracy.f90

inout.o:
	$(FC) $(FCOPTS) -fopenmp -c inout.f90

#Cleaning-Options
.PHONY: clean realclean

clean:
	rm -f *.mod *.o *.dat

realclean: clean
	rm -f crwalk
