#
COMPILER = gfortran
OPTFLAGS = -O3
DBGFLAGS = -g

All: Modules.o Main.o PotentialProfile.o fitpack.o MoveParticlePack.o CoulombCollisions.o
	gfortran  $(OPTFLAGS) -fopenmp Modules.o Main.o PotentialProfile.o fitpack.o MoveParticlePack.o CoulombCollisions.o -o MPEX
	rm *.o *.mod

Modules.o: Modules.f90
	$(COMPILER) $(OPTFLAGS) -c Modules.f90

Main.o: Main.f90
	$(COMPILER) $(OPTFLAGS) -fopenmp -c Main.f90

PotentialProfile.o: PotentialProfile.f90
	$(COMPILER) $(OPTFLAGS) -c PotentialProfile.f90

fitpack.o: fitpack.f
	$(COMPILER) $(OPTFLAGS) -c fitpack.f

MoveParticlePack.o: MoveParticlePack.f90
	$(COMPILER) $(OPTFLAGS) -c MoveParticlePack.f90

CoulombCollisions.o: CoulombCollisions.f90
	$(COMPILER) $(OPTFLAGS) -fopenmp -c CoulombCollisions.f90

clean:
	rm *.o *.mod mpex
