SUBROUTINE AdvanceEfield(mesh,params)
USE LOCAL
USE dataTYP
USE PhysicalConstants

IMPLICIT NONE	

! Declare interface variables:
TYPE(meshTYP),   INTENT(INOUT) :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Declare local variables:
REAL(r8), DIMENSION(mesh%NZmesh + 4) :: E, n
REAL(r8), DIMENSION(mesh%NZmesh) :: dn
REAL(r8) :: Te, dz
INTEGER(i4) :: i, NZ, NZg, frame
INTEGER(i4), DIMENSION(mesh%NZmesh) :: rng

! Input:
n   = mesh%n
NZ  = mesh%NZmesh
NZg = NZ + 4
dz  = mesh%dzm
Te  = params%Te0

! Plasma density gradient:
rng = (/ (i, i = 3, (NZ + 2)) /)
dn = (n(rng+1) - n(rng-1))/(2*dz) 

! Electric field:
E = 0.
E(rng) = -Te*dn/n(rng)

! Apply smoothing:
frame = 9
CALL MovingMean_F(E,NZg,frame)

! Output: 
mesh%E = E

RETURN
END SUBROUTINE AdvanceEfield

! =======================================================================================================
SUBROUTINE MovingMean_F(y,NX,k)
! =======================================================================================================
USE LOCAL
USE dataTYP

IMPLICIT NONE

! Define interface variables:
REAL(r8)   , DIMENSION(NX), INTENT(INOUT) :: y
INTEGER(i4), INTENT(IN) :: NX
INTEGER(i4), INTENT(INOUT) :: k

! Define local variables:
REAL(r8), DIMENSION(NX) :: ym
REAL(r8) :: yd
INTEGER(i4) :: s, ii, jj, istart, iend, N

! Check frame:
IF (MOD(k,2) .EQ. 0) THEN
        k = k + 1
END IF

! Half frame size:
s = (k-1)/2

DO ii = 1,NX
        ! Start and end of frame:
        istart = ii-s
        iend   = ii+s

        ! Correct frame at edges:
        IF (istart .LE. 0 ) istart = 1
        IF (iend   .GT. NX) iend   = NX

        ! Effective frame size:
        N = iend - istart + 1

        ! Calculate mean:
        yd = 0.
        DO jj = istart,iend
                yd = yd + y(jj)/N
        END DO
        ym(ii) = yd
END DO

y = ym

END SUBROUTINE MovingMean_F

