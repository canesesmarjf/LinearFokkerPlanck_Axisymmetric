#
COMPILER = gfortran
OPTFLAGS = -O3
DBGFLAGS = -g
OBJ_1 = Modules.o linFP.o PotentialProfile.o
OBJ_2 = fitpack.o MoveParticlePack.o CoulombCollisions.o

All: $(OBJ_1) $(OBJ_2)
	$(COMPILER) $(OPTFLAGS) -fopenmp $(OBJ_1) $(OBJ_2) -o linFP
	rm *.o *.mod

Modules.o: Modules.f90
	$(COMPILER) $(OPTFLAGS) -c Modules.f90

linFP.o: linFP.f90
	$(COMPILER) $(OPTFLAGS) -fopenmp -c linFP.f90

PotentialProfile.o: PotentialProfile.f90
	$(COMPILER) $(OPTFLAGS) -c PotentialProfile.f90

fitpack.o: fitpack.f
	$(COMPILER) $(OPTFLAGS) -c -w fitpack.f

MoveParticlePack.o: MoveParticlePack.f90
	$(COMPILER) $(OPTFLAGS) -c -fopenmp MoveParticlePack.f90

CoulombCollisions.o: CoulombCollisions.f90
	$(COMPILER) $(OPTFLAGS) -fopenmp -c CoulombCollisions.f90

clean:
	rm *.o *.mod linFP *.dat
