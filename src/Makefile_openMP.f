#
COMPILER = gfortran
OPTFLAGS = -O3
DBGFLAGS = -g
OBJ_1 = Modules.o Main.o PotentialProfile.o
OBJ_2 = fitpack.o MoveParticlePack.o CoulombCollisions.o

All: $(OBJ_1) $(OBJ_2)
	gfortran  $(OPTFLAGS) -fopenmp $(OBJ_1) $(OBJ_2) -o MPEX
	rm *.o *.mod

Modules.o: Modules.f90
	$(COMPILER) $(OPTFLAGS) -c Modules.f90

Main.o: Main.f90
	$(COMPILER) $(OPTFLAGS) -fopenmp -c Main.f90

PotentialProfile.o: PotentialProfile.f90
	$(COMPILER) $(OPTFLAGS) -c PotentialProfile.f90

fitpack.o: fitpack.f
	$(COMPILER) $(OPTFLAGS) -c -w fitpack.f

MoveParticlePack.o: MoveParticlePack.f90
	$(COMPILER) $(OPTFLAGS) -c MoveParticlePack.f90

CoulombCollisions.o: CoulombCollisions.f90
	$(COMPILER) $(OPTFLAGS) -fopenmp -c CoulombCollisions.f90

clean:
	rm *.o *.mod MPEX *.dat
