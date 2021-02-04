#
COMPILER=gfortran
DBGFLAGS=-g
OPTFLAGS=-O3
OBJ_1 = Modules.o test_spline.o
OBJ_2 = fitpack.o MoveParticlePack.o

All: $(OBJ_1) $(OBJ_2)
	$(COMPILER) $(OPTFLAGS) -fopenmp $(OBJ_1) $(OBJ_2) -o test_spline
	rm *.o *.mod

Modules.o: Modules.f90
	$(COMPILER) $(OPTFLAGS) -c Modules.f90

test_spline.o: test_spline.f90
	$(COMPILER) $(OPTFLAGS) -fopenmp -c test_spline.f90

fitpack.o: fitpack.f
	$(COMPILER) $(OPTFLAGS) -c -w fitpack.f

MoveParticlePack.o: MoveParticlePack.f90
	$(COMPILER) $(OPTFLAGS) -c MoveParticlePack.f90

clean:
	rm *.o *.mod test_spline *.dat
