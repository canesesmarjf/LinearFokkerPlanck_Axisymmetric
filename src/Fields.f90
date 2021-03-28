SUBROUTINE AdvanceEfield(plasma,mesh,params)
USE LOCAL
USE dataTYP
USE PhysicalConstants

IMPLICIT NONE	

! Declare interface variables:
TYPE(plasmaTYP), INTENT(IN)    :: plasma
TYPE(meshTYP),   INTENT(INOUT) :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Declare local variables:
REAL(r8), DIMENSION(mesh%NZmesh + 4) :: E, n
REAL(r8), DIMENSION(mesh%NZmesh) :: dn
REAL(r8) :: Te, dz
INTEGER(i4) :: i, NZ, NZg, frame
INTEGER(i4) , DIMENSION(mesh%NZmesh) :: rng

! Input:
n   = mesh%n
NZ  = mesh%NZmesh
NZg = NZ + 4
dz  = mesh%dzm
Te  = params%Te0

! Plasma density gradient:
rng = (\ (i, i = 3, (NZ + 2)) \)
dn = (n(rng+1) - n(rng-1))/(2*dz) 

! Electric field:
E = 0.
E(rng) = -Te*dn/n(rng)

! Apply smoothing:
frame = 9
CALL MovingMean(E(rng),NZ,frame)

! Output: 
mesh%E = E

RETURN
END SUBROUTINE AdvanceEfield
