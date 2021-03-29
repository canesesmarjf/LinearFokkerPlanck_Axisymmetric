#
COMPILER = gfortran
OPTFLAGS = -O3 -fopenmp
DBGFLAGS = -g
OBJ = T1.o T2.o T3.o T4.o T5.o T6.o
PROG = linFP

All: $(OBJ)
	$(COMPILER) $(OPTFLAGS) $(OBJ) -o $(PROG)
	rm *.o *.mod

T1.o: Modules.f90
	$(COMPILER) $(OPTFLAGS) -c Modules.f90 -o T1.o

T2.o: linFP.f90
	$(COMPILER) $(OPTFLAGS) -c linFP.f90 -o T2.o

T3.o: MoveParticlePack.f90
	$(COMPILER) $(OPTFLAGS) -c MoveParticlePack.f90 -o T3.o

T4.o: CoulombCollisions.f90
	$(COMPILER) $(OPTFLAGS) -c CoulombCollisions.f90 -o T4.o

T5.o: PIC.f90
	$(COMPILER) $(OPTFLAGS) -c PIC.f90 -o T5.o

T6.o: Fields.f90 MoveParticlePack.f90
	$(COMPILER) $(OPTFLAGS) -c Fields.f90 -o T6.o

clean:
	rm *.o *.mod linFP *.dat
