! =======================================================================================================
SUBROUTINE AdvanceParticles(plasma,mesh,fieldspline,params)
! =======================================================================================================
USE local
USE dataTYP
USE OMP_LIB

IMPLICIT NONE

! Declare interface variables:
TYPE(plasmaTYP)     , INTENT(INOUT) :: plasma
TYPE(meshTYP)       , INTENT(IN)    :: mesh
TYPE(paramsTYP)     , INTENT(IN)    :: params
TYPE(fieldSplineTYP), INTENT(IN)    :: fieldspline

! Declare local variables:
REAL(r8) :: dresNum, resNum0, resNum1
INTEGER(i4) :: i

! Initialize resonance number:
dresNum = 0.
resNum0 = 0.
resNum1 = 0.

! Advance particle positions and velocities:
!$OMP PARALLEL DO FIRSTPRIVATE(dresNum, resNum0, resNum1) SCHEDULE(STATIC)
DO i = 1,params%NC
	! 2.1 - Reset flags:
	CALL ResetFlags(i,plasma)

	! 2.2- Compute resonance number:
	IF (params%iHeat) THEN
		CALL CyclotronResonanceNumber(i,plasma,fieldspline,params,resNum0)
        END IF

	! 2.3- Advance zp, kep, xip:
	IF (params%iPush) THEN
		! CALL  MoveParticle(i,plasma,fieldspline,params)
		CALL MoveParticle_RK4(i,plasma,mesh,params)
  END IF

	! 2.4- Check boundaries: set f1 and f2 flags
	IF (params%iPush) THEN
		CALL CheckBoundary(i,plasma,params)
	END IF

	! 2.5- Compute resonance number:
	IF (params%iHeat) THEN
		CALL CyclotronResonanceNumber(i,plasma,fieldspline,params,resNum0)
        END IF

	! 2.6- Check resonance: set f3 flag
	IF (params%iHeat) THEN
		dresNum = dsign(1.d0,resNum0*resNum1)
		IF (dresNum .LT. 0 .AND. plasma%zp(i) .GT. params%zRes1 .AND. plasma%zp(i) .LT. params%zRes2) THEN
        		plasma%f3(i) = 1
			CALL RFoperatorTerms(i,plasma,fieldspline,params)
 		END IF
        END IF
END DO
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE AdvanceParticles

! ======================================================================================================
SUBROUTINE InterpolateElectromagneticFields(plasma,mesh,params)
! ======================================================================================================
USE LOCAL
USE dataTYP

IMPLICIT NONE

! Declare interface variables:
TYPE(plasmaTYP), INTENT(INOUT) :: plasma
TYPE(meshTYP)  , INTENT(IN)    :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Declare local variables:
INTEGER(i4) :: i, ix
REAL(r8), DIMENSION(3) ::  w, E, B, dB, ddB

!$OMP PARALLEL DO PRIVATE(ix, w, E, B, dB, ddB)
DO i = 1,params%NC
	IF (plasma%f1(i) .EQ. 0 .AND. plasma%f2(i) .EQ. 0) THEN
                ! Get nearest grid point:
                ix = plasma%m(i) + 2

                ! Assignment function:
                w(1) = plasma%wL(i)
                w(2) = plasma%wC(i)
                w(3) = plasma%wR(i)

                ! Electric field:
                E(1) = mesh%E(ix - 1)
                E(2) = mesh%E(ix)
                E(3) = mesh%E(ix + 1)
                plasma%Ep(i) = DOT_PRODUCT(w,E)

                ! Magnetic field:
                B(1) = mesh%B(ix - 1)
                B(2) = mesh%B(ix)
                B(3) = mesh%B(ix + 1)
                plasma%Bp(i) = DOT_PRODUCT(w,B)

                ! Magnetic field: 1st derivative
                dB(1) = mesh%dB(ix - 1)
                dB(2) = mesh%dB(ix)
                dB(3) = mesh%dB(ix + 1)
                plasma%dBp(i) =  DOT_PRODUCT(w,dB)

                ! Magnetic field: 2nd derivative
		!IF (params%iHeat) THEN
		IF (.TRUE.) THEN
	                ddB(1) = mesh%ddB(ix - 1)
        	        ddB(2) = mesh%ddB(ix)
	        	      ddB(3) = mesh%ddB(ix + 1)
                	plasma%ddBp(i) =  DOT_PRODUCT(w,ddB)
		END IF
	END IF
END DO
!$OMP END PARALLEL DO

END SUBROUTINE InterpolateElectromagneticFields

! =======================================================================================================
SUBROUTINE InterpolateMoments(plasma,mesh,params)
! =======================================================================================================
USE LOCAL
USE dataTYP

IMPLICIT NONE

! Declare interface variables:
TYPE(plasmaTYP), INTENT(INOUT) :: plasma
TYPE(meshTYP)  , INTENT(INOUT) :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Declare local variables:
INTEGER(i4) :: i

END SUBROUTINE InterpolateMoments

! =======================================================================================================
SUBROUTINE CheckBoundary(i,plasma,params)
! =======================================================================================================
USE LOCAL
USE dataTYP

IMPLICIT NONE

! Declare interface variables:
INTEGER(i4),     INTENT(IN)    :: i
TYPE(plasmaTYP), INTENT(INOUT) :: plasma
TYPE(paramsTYP), INTENT(IN)    :: params

IF (plasma%zp(i) .GE. params%zmax) THEN
	plasma%f2(i)  = 1
	plasma%dE2(i) = plasma%kep(i)
ELSE IF (plasma%zp(i) .LE. params%zmin) THEN
	plasma%f1(i)  = 1
	plasma%dE1(i) = plasma%kep(i)
END IF

RETURN
END SUBROUTINE CheckBoundary


! =======================================================================================================
SUBROUTINE MoveParticle(i,plasma,fieldspline,params)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE

! Define type for interface arguments
INTEGER(i4)         , INTENT(IN)    :: i
TYPE(plasmaTYP)     , INTENT(INOUT) :: plasma
TYPE(paramsTYP)     , INTENT(IN)    :: params
TYPE(fieldSplineTYP), INTENT(IN)    :: fieldspline

! Local variables:
REAL(r8) :: zpnew, Xipnew, uparnew, upernew, munew  ! Position, kinetic energy and pitch of the ith particle
REAL(r8) :: Ma, dt

! Storage for the MoverParticle and RHS subroutines:
REAL(r8) :: xip0, xip1
REAL(r8) :: kep0, kep1
REAL(r8) :: K1, K2, K3, K4
REAL(r8) :: L1, L2, L3, L4
REAL(r8) :: zp0  , zp1  , zp2  , zp3
REAL(r8) :: upar0, upar1, upar2, upar3
REAL(r8) :: mu0  , mu1  , mu2  , mu3
REAL(r8) :: M1   , M2   , M3   , M4
REAL(r8) :: u2
REAL(r8) :: B

! Input data:
zp0  = plasma%zp(i)
kep0 = plasma%kep(i)
xip0 = plasma%xip(i)

! Time step:
dt = params%dt

! Test particle mass:
Ma = params%Ma

! Calculate initial parallel particle speed:
upar0 = xip0*sqrt(2.*e_c*kep0/Ma)

! Calculate initial particle speed squared:
u2 = 2.*e_c*kep0/Ma

! Calculate initial magnetic moment:
CALL Interp1(zp0,B,fieldspline%B)

mu0 = 0.5*Ma*u2*(1 - xip0*xip0)/B

! Begin assembling RK4 solution:
call RightHandSide(zp0,upar0,mu0,K1,L1,M1,params,fieldspline)             ! Update values of fK, fL and fM
zp1   = zp0   + (K1*dt/2.)
upar1 = upar0 + (L1*dt/2.)
mu1   = mu0   + (M1*dt/2.)

call RightHandSide(zp1,upar1,mu1,K2,L2,M2,params,fieldspline)             ! Update values of fK, fL and fM
zp2   = zp0   + (K2*dt/2.)
upar2 = upar0 + (L2*dt/2.)
mu2   = mu0   + (M2*dt/2.)

call RightHandSide(zp2,upar2,mu2,K3,L3,M3,params,fieldspline)             ! Update values of fK, fL and fM
zp3   = zp0   + K3*dt
upar3 = upar0 + L3*dt
mu3   = mu0   + M3*dt

call RightHandSide(zp3,upar3,mu3,K4,L4,M4,params,fieldspline)             ! Update values of fK, fL and fM
zpnew   = zp0   + ( (K1 + (2.*K2) + (2.*K3) + K4)/6. )*dt
uparnew = upar0 + ( (L1 + (2.*L2) + (2.*L3) + L4)/6. )*dt
munew   = mu0 +   ( (M1 + (2.*M2) + (2.*M3) + M4)/6. )*dt

! Calculate the magnetic field at zpnew:
CALL Interp1(zpnew,B,fieldspline%B)

! Based on new B and new mu, calculate new uper:
upernew = sqrt(2.*munew*B/Ma)

! New u due to new upar and new uper:
u2 = uparnew**2. + upernew**2.

! New pitch angle due to new upar and new u:
Xipnew = uparnew/sqrt(u2)

! Output data:
plasma%zp(i)  = zpnew
plasma%xip(i) = Xipnew
plasma%kep(i) = 0.5*Ma*u2/e_c

RETURN
END SUBROUTINE MoveParticle

! =======================================================================================================
SUBROUTINE RightHandSide(zp0,upar0,mu0,K,L,M,params,fieldspline)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE

! Define type of interface arguments:
REAL(r8)            , INTENT(IN) :: zp0, upar0, mu0
REAL(r8)            , INTENT(OUT) :: K, L, M
TYPE(paramsTYP)     , INTENT(IN) :: params
TYPE(fieldSplineTYP), INTENT(IN) :: fieldspline

! Local variables:
REAL(r8) :: dB, E
REAL(r8) :: Ma, qa

! Test particle mass:
Ma = params%Ma

! Test particle charge:
qa = params%qa

! Calculate the magnetic field gradient at zp0:
CALL Interp1(zp0,dB,fieldspline%dB)

! Calculate electric field at zp0:
CALL Interp1(zp0,E,fieldspline%E)

! Assign values to output variables:
K = upar0
L = (-mu0*dB + qa*E)/Ma
M = 0.

RETURN
END SUBROUTINE RightHandSide

! =======================================================================================================
SUBROUTINE AssignCell(plasma,mesh,params)
! =======================================================================================================
USE LOCAL
USE dataTYP
USE OMP_LIB

IMPLICIT NONE

! Define interface variables:
TYPE(plasmaTYP), INTENT(INOUT) :: plasma
TYPE(meshTYP)  , INTENT(IN)    :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Declare local variables:
INTEGER(i4) :: i
REAL(r8)    :: Z, dz, zi, zoffset

! Mesh grid size:
dz      = mesh%dzm
! Mesh offset:
zoffset = mesh%zmin

!$OMP PARALLEL DO PRIVATE(zi,Z)
DO i = 1,params%NC
	IF (plasma%f1(i) .EQ. 0 .AND. plasma%f2(i) .EQ. 0) THEN

	! Using zm and zp, find nearest grid point (NGP):
	zi = plasma%zp(i)
	plasma%m(i) = NINT(0.5 + (zi - zoffset)/dz)

	! Compute wL(i), wC(i), wR(i)
	Z = mesh%zm(plasma%m(i)) - zi
	plasma%wC(i) = 0.75 - (Z/dz)**2.
	plasma%wL(i) = 0.5*(1.5 + ((Z - dz)/dz) )**2.
	plasma%wR(i) = 0.5*(1.5 - ((Z + dz)/dz) )**2.

	END IF
END DO
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE AssignCell

! =======================================================================================================
SUBROUTINE ExtrapolateMomentsToMesh(plasma,mesh,params)
! =======================================================================================================
USE LOCAL
USE PhysicalConstants
USE dataTYP
USE OMP_LIB

IMPLICIT NONE

! Define interface variables:
TYPE(plasmaTYP), INTENT(IN)    :: plasma
TYPE(meshTYP)  , INTENT(INOUT) :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Define local variables:
INTEGER(i4) :: i, ix, frame, NZ
INTEGER(i4), DIMENSION(3) :: ixLeft, ixRight
REAL(r8) :: vpar, Ma, Ep, alpha, xip, vper, v, a, phi, vx
REAL(r8), DIMENSION(mesh%NZmesh + 4) :: n, nU, unU, P11, P22, nUE, cf

! Initialize local mesh quantities:
n = 0.
nU = 0.
unU = 0.
P11 = 0.
P22 = 0.
nUE = 0.
! Species "a" mass:
Ma = params%Ma
! Scaling factor:
alpha = plasma%alpha
! Size of mesh quantities:
NZ = mesh%NZmesh
! Set the edge cell number:
ixLeft  = (/1,2,3/)
ixRight = (/2,3,4/) + NZ

! Calculate moments and extrapolate to mesh points
!$OMP PARALLEL DO PRIVATE(a,Ep,v,xip,vpar,vper,ix) REDUCTION(+:n, nU, unU, P11, P22, nUE)
DO i = 1,params%NC
	IF (plasma%f1(i) .EQ. 0 .AND. plasma%f2(i) .EQ. 0) THEN

		! Derived quantities:
		a    = plasma%a(i)
		Ep   = e_c*plasma%kep(i)
		v    = sqrt(2.*Ep/Ma)
		xip  = plasma%xip(i)
		vpar = v*xip
		vper = v*sqrt(1 - xip**2.)
		call random_number(phi)
		phi = phi*2*pi
		vx = vper*cos(phi)

		! Mesh point with ghost cells:
		ix = plasma%m(i) + 2

		! Particle density: [NSP m^-3]
		n(ix-1)  = n(ix-1)    + plasma%wL(i)*a
		n(ix)    = n(ix)      + plasma%wC(i)*a
		n(ix+1)  = n(ix+1)    + plasma%wR(i)*a

		! Particle flux density: [NSP m^-2 s^-1]
		nU(ix-1)  = nU(ix-1)  + plasma%wL(i)*a*vpar
		nU(ix)    = nU(ix)    + plasma%wC(i)*a*vpar
		nU(ix+1)  = nU(ix+1)  + plasma%wR(i)*a*vpar

		! Computational particle flux density: [NC m^-2 s^-1]
		unU(ix-1) = unU(ix-1) + plasma%wL(i)*vpar
		unU(ix)   = unU(ix)   + plasma%wC(i)*vpar
		unU(ix+1) = unU(ix+1) + plasma%wR(i)*vpar

		! Stress tensor 11:
		P11(ix-1)  = P11(ix-1)  + plasma%wL(i)*a*Ma*vpar**2.
		P11(ix)    = P11(ix)    + plasma%wC(i)*a*Ma*vpar**2.
		P11(ix+1)  = P11(ix+1)  + plasma%wR(i)*a*Ma*vpar**2.

		! Stress tensor 22:
		P22(ix-1)  = P22(ix-1)  + plasma%wL(i)*a*Ma*vx**2.
		P22(ix)    = P22(ix)    + plasma%wC(i)*a*Ma*vx**2.
		P22(ix+1)  = P22(ix+1)  + plasma%wR(i)*a*Ma*vx**2.

		! Energy flux density:
		nUE(ix-1)  = nUE(ix-1)  + plasma%wL(i)*a*vpar*Ep
		nUE(ix)    = nUE(ix)    + plasma%wC(i)*a*vpar*Ep
		nUE(ix+1)  = nUE(ix+1)  + plasma%wR(i)*a*vpar*Ep

	END IF
END DO
!$OMP END PARALLEL DO

! Correct edge cells of particle density:
n(ixLeft)  = n(4)
n(ixRight) = n(NZ+1)

! Apply magnetic compression:
cf = (mesh%B/mesh%B0)/params%Area0
mesh%n   = n*cf
mesh%nU  = nU*cf
mesh%unU = unU
mesh%P11 = P11*cf
mesh%P22 = P22*cf
mesh%nUE = nUE*cf

! Apply scaling factor:
mesh%n   = alpha*  mesh%n/mesh%dzm
mesh%nU  = alpha* mesh%nU/mesh%dzm
mesh%unU = alpha*mesh%unU/mesh%dzm
mesh%P11 = alpha*mesh%P11/mesh%dzm
mesh%P22 = alpha*mesh%P22/mesh%dzm
mesh%nUE = alpha*mesh%nUE/mesh%dzm

! Apply smoothing:
frame = 9
!$OMP PARALLEL
IF (OMP_GET_THREAD_NUM() .EQ. 0) CALL MovingMean(mesh%n(3:(NZ+2))  ,NZ,frame)
IF (OMP_GET_THREAD_NUM() .EQ. 1) CALL MovingMean(mesh%nU(3:(NZ+2)) ,NZ,frame)
IF (OMP_GET_THREAD_NUM() .EQ. 2) CALL MovingMean(mesh%unU(3:(NZ+2)),NZ,frame)
IF (OMP_GET_THREAD_NUM() .EQ. 3) CALL MovingMean(mesh%P11(3:(NZ+2)),NZ,frame)
IF (OMP_GET_THREAD_NUM() .EQ. 4) CALL MovingMean(mesh%P22(3:(NZ+2)),NZ,frame)
IF (OMP_GET_THREAD_NUM() .EQ. 5) CALL MovingMean(mesh%nUE(3:(NZ+2)),NZ,frame)
!$OMP END PARALLEL

! Calculate U:
mesh%U    = mesh%nU/mesh%n

! Calculate Ppar and Pper:
mesh%Ppar = mesh%P11 - Ma*mesh%nU*mesh%U
mesh%Pper = mesh%P22

! Calculate Tpar and Tper:
mesh%Tpar = mesh%Ppar/(e_c*mesh%n)
mesh%Tper = mesh%Pper/(e_c*mesh%n)

END SUBROUTINE ExtrapolateMomentsToMesh

! =======================================================================================================
SUBROUTINE MovingMean(y,NX,k)
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

END SUBROUTINE MovingMean

! =======================================================================================================
SUBROUTINE ReinjectParticles(i,plasma,params)
! =======================================================================================================
USE local
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define interface variables:
INTEGER(i4)    , INTENT(IN)    :: i
TYPE(plasmaTYP), INTENT(INOUT) :: plasma
TYPE(paramsTYP), INTENT(IN)    :: params

! Define local variables:
REAL(r8) :: zp0, kep0, xip0
REAL(r8), DIMENSION(6) :: Rm6
REAL(r8) :: Ma, T0, T, vT, sigma_v, E, U, Ux, Uy, Uz
REAL(r8) :: R_1, R_3, t_2, t_4
REAL(r8) :: wx, wy, wz, vx, vy, vz, v

! Input variables:
zp0  = plasma%zp(i)
kep0 = plasma%kep(i)
xip0 = plasma%xip(i)

! Choose particle BC:
! 1: Isotropic plasma source
! 2: NBI
! 3: Periodic boundary
IF (params%BC_Type .EQ. 1 .OR. params%BC_Type .EQ. 2) THEN

  IF (params%BC_Type .EQ. 1) THEN
    T = params%Te0
    E = 0
  ELSE IF (params%BC_Type .EQ. 2) THEN
    T = params%BC_Tp
    E = params%BC_Ep
  END IF

  ! Velocity distribution:
  Ma = params%Ma
  vT = sqrt(2.*e_c*T/Ma)
  U  = sqrt(2.*e_c*E/Ma)
  Ux = 0
  Uy = U*sqrt(1. - params%BC_xip**2.)
  Uz = U*params%BC_xip
  sigma_v = vT/sqrt(2.)

  ! Box-muller terms:
  call random_number(Rm6)
  R_1 = sigma_v*sqrt(-2.*log(Rm6(1)))
  t_2 = 2.*pi*Rm6(2)
  R_3 = sigma_v*sqrt(-2.*log(Rm6(3)))
  t_4 = 2.*pi*Rm6(4)

  wx = R_1*cos(t_2)
  wy = R_1*sin(t_2)
  wz = R_3*cos(t_4)

  vx = Ux + wx
  vy = Uy + wy
  vz = Uz + wz
  v = sqrt( vx**2. + vy**2. + vz**2. )

  ! Populate output variables:
  plasma%kep(i) = 0.5*(Ma/e_c)*v**2
  plasma%xip(i) = vz/v

  ! Position distribution:
  plasma%zp(i) = params%BC_zp_std*sqrt(-2.*log(Rm6(5)))*cos(2.*pi*Rm6(6)) + params%BC_zp_mean

ELSE if (params%BC_Type .EQ. 3) THEN
  IF (zp0 .GE. params%zmax) THEN
    plasma%zp(i) = params%zmin
  ELSE
    plasma%zp(i) = params%zmax
  END IF
END IF

RETURN
END SUBROUTINE ReinjectParticles

! =======================================================================================================
SUBROUTINE CyclotronResonanceNumber(i,plasma,fieldspline,params,resNum0)
! =======================================================================================================

USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define interface variables:
INTEGER(i4)            , INTENT(IN)    :: i
TYPE(plasmaTYP)        , INTENT(IN)    :: plasma
REAL(r8)               , INTENT(INOUT) :: resNum0
TYPE(paramsTYP)        , INTENT(IN)    :: params
TYPE(fieldSplineTYP)   , INTENT(IN)    :: fieldspline

! Define local variables
REAL(r8) :: zp0, kep0, xip0
REAL(r8) :: upar, Bf, Omega, Omega_RF
REAL(r8) :: Ma, qa

! Input:
zp0  = plasma%zp(i)
kep0 = plasma%kep(i)
xip0 = plasma%xip(i)
! Test particle mass:
Ma = params%Ma
! Test particle charge:
qa = params%qa
! Parallel velocity of test particle:
upar = sqrt(2.*e_c*kep0/Ma)*xip0
! Magnetic field at location zp0 of test particle:
CALL Interp1(zp0,Bf,fieldspline%B)
! Cyclotron frequency of test particle:
Omega = abs(qa)*Bf/Ma
! RF frequency in rad/s:
Omega_RF = 2*pi*params%f_RF
! Cyclotron resonance number:
resNum0 = Omega_RF - params%kpar*upar - params%n_harmonic*Omega

RETURN
END SUBROUTINE CyclotronResonanceNumber

! =======================================================================================================
SUBROUTINE RFoperatorTerms(i,plasma,fieldspline,params)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define interface variables:
INTEGER(i4)            , INTENT(IN)    :: i
TYPE(plasmaTYP)        , INTENT(INOUT) :: plasma
TYPE(paramsTYP)        , INTENT(IN)    :: params
TYPE(fieldSplineTYP)   , INTENT(IN)    :: fieldspline

! Define local variables:
REAL(r8) :: zp0, kep0, xip0
REAL(r8) :: u0, upar0, uper0
REAL(r8) :: kep_par0, kep_per0
REAL(r8) :: dB, ddB, E
REAL(r8) :: Bf, Omega, dOmega, ddOmega
REAL(r8) :: Omega_dot, Omega_ddot, tau_rf
REAL(r8) :: rl, flr, besselterm
REAL(r8) :: mean_dkep_per
REAL(r8) :: Ma, qa

! Input variables:
zp0  = plasma%zp(i)
kep0 = plasma%kep(i)
xip0 = plasma%xip(i)

! Test particle mass:
Ma = params%Ma

! Test particle charge:
qa = params%qa

! Calculate derived quantities:
u0       = sqrt(2.*e_c*kep0/Ma)
upar0    = u0*xip0
uper0    = u0*(1. - xip0**2)**0.5
kep_par0 = kep0*xip0**2.
kep_per0 = kep0*(1. - xip0**2.)

! Gradients:
CALL Interp1(zp0,dB ,fieldspline%dB )
CALL Interp1(zp0,ddB,fieldspline%ddB)
CALL Interp1(zp0,E  ,fieldspline%E  )

! Spatial derivatives of the magnetic field:
CALL Interp1(zp0,Bf,fieldspline%B)
Omega     = params%n_harmonic*e_c*Bf/Ma
dOmega    = params%n_harmonic*e_c*dB/Ma
ddOmega   = params%n_harmonic*e_c*ddB/Ma

! Calculate the first and second time derivative of Omega:
Omega_dot = upar0*dOmega
Omega_ddot = (upar0**2.)*ddOmega  - (uper0**2.)*dOmega*dOmega/(2.*Omega) + qa*E*dOmega/Ma

! Calculate the interaction time (tau_RF):
IF ( (Omega_ddot**2.) .GT. 4.8175*ABS(Omega_dot**3.) )  then
        ! Approximate Ai(x) ~ 0.3833
        tau_rf = (2.*pi)*(ABS(2./Omega_ddot)**(1/3.))*0.3833
ELSE
        tau_rf = sqrt(2.*pi/ABS(Omega_dot))
END IF

! Calculate Bessel term:
rl       = uper0/(abs(qa)*Bf/Ma)
flr      = params%kper*rl
besselterm = BESSEL_JN(params%n_harmonic-1,flr)

!Calculate the mean RF kick per unit electric field squared:
mean_dkep_per = 0.5*(e_c/Ma)*(besselterm*tau_rf)**2.

! Populate plasma structure:
plasma%udErf(i)   = mean_dkep_per
plasma%doppler(i) = params%kpar*upar0/Omega
plasma%udE3(i)    = plasma%udErf(i)*(1. + plasma%doppler(i))

RETURN
END SUBROUTINE RFoperatorTerms

! =======================================================================================================
SUBROUTINE RFOperator(i,plasma,fieldspline,params,uE3)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP
USE OMP_LIB

IMPLICIT NONE
! Define interface variables:
REAL(r8)               , INTENT(IN)    :: uE3
INTEGER(i4)            , INTENT(IN)    :: i
TYPE(plasmaTYP)        , INTENT(INOUT) :: plasma
TYPE(paramsTYP)        , INTENT(IN)    :: params
TYPE(fieldSplineTYP)   , INTENT(IN)    :: fieldspline

! Define local variables:
REAL(r8) :: kep0, xip0, Ew2
REAL(r8) :: kep_par0, kep_per0
REAL(r8) :: mean_dkep_per, dkep_per, Rm1
REAL(r8) :: dkep_par, dkep, kep1
REAL(r8) :: kep_per1, kep_par1
REAL(r8) :: upar1, u1
REAL(r8) :: Ma

! Input variables:
kep0 = plasma%kep(i)
xip0 = plasma%xip(i)

! Test particle mass:
Ma = params%Ma

! Initial energy state:
kep_par0 = kep0*xip0**2.
kep_per0 = kep0*(1. - xip0**2.)

! Calculate the mean RF kick:
Ew2 = params%Prf/uE3
mean_dkep_per = plasma%udErf(i)*Ew2

IF (OMP_GET_THREAD_NUM() .EQ. 0) THEN
!   WRITE(*,*) 'Ew', sqrt(Ew2)
END IF


! Calculate the change in perp, parallel and total energy:
CALL RANDOM_NUMBER(Rm1)
Rm1 = 2.*Rm1 - 1.
dkep_per = mean_dkep_per + Rm1*sqrt(2.*kep_per0*mean_dkep_per)

! Given the perp kick in energy, apply the parallel energy kick
! This arises from the non-linear effect of the perturbed magnetic field
! See Stix section 10.3 "Trapped electromagnetic modes"
dkep_par = (plasma%doppler(i))*dkep_per

! total change in energy kick
dkep = dkep_par + dkep_per

! Calculate changes in energy
! Perp degree of freedon:
kep_per1 = kep_per0 + dkep_per
! Parallel degree of freedom:
kep_par1 = kep_par0 + dkep_par
! Total final energy:
kep1     = kep_per1 + kep_par1

if (kep1 .LT. 0) then
  WRITE(*,*) 'final kep < 0'
  print *, 'kep0', kep0
  print *, 'dkep', dkep
  print *, 'kep1', kep1
end if

! Assign value of the new energy:
plasma%kep(i) = kep1

! Calculate the new pitch angle:
upar1 = sqrt( (2.*e_c/Ma)*abs(kep_par1) )*dsign(1.d0,xip0)*dsign(1.d0,kep_par1)
u1    = sqrt( (2.*e_c/Ma)*kep1 )
plasma%xip(i) = upar1/u1

! Record resonance event:
plasma%f3(i) = 1
! Record energy kick:
plasma%dE3(i) = dkep

RETURN
END SUBROUTINE RFOperator

! =======================================================================================================
SUBROUTINE loadParticles(i,plasma,params)
! =======================================================================================================
  USE local
  use PhysicalConstants
  USE dataTYP

  IMPLICIT NONE
  ! Define interface variables:
  INTEGER(i4)    , INTENT(IN)    :: i
  TYPE(plasmaTYP), INTENT(INOUT) :: plasma
  TYPE(paramsTYP), INTENT(IN)    :: params

  ! Declare internal variables:
  REAL(r8) :: zp0, kep0, xip0
  REAL(r8) :: X1, X2, X3, X4
  REAL(r8) :: R_1, R_2, R_3, R_4, t_2, t_4
  REAL(r8) :: wx, wy, wz, vx, vy, vz, v
  REAL(r8) :: zmin, zmax, sigma_v, Ma
  REAL(r8) :: Ux, Uy, Uz, vT, U, T, E

  ! Input variables:
  zp0  = plasma%zp(i)
  kep0 = plasma%kep(i)
  xip0 = plasma%xip(i)

  ! Particle position:
  zmin = params%zmin
  zmax = params%zmax
  if (params%IC_Type .EQ. 1) then
      ! Uniform load
      call random_number(X1)
      plasma%zp(i) = zmin + (zmax - zmin)*X1
  elseif (params%IC_Type .EQ. 2) then
      ! Gaussian load
      call random_number(X1)
      call random_number(X2)
      plasma%zp(i) = params%IC_zp_std*sqrt(-2.*log(X1))*cos(2.*pi*X2)  +  params%IC_zp_mean
  end if

  ! Particle kinetic energy and pitch angle:
  ! Generate random numbers:
  call random_number(X1)
  call random_number(X2)
  call random_number(X3)
  call random_number(X4)

  ! Derived quantities:
  Ma = params%Ma
  T = params%IC_Tp
  E = params%IC_Ep
  vT = sqrt(2.*e_c*T/Ma)
  U  = sqrt(2.*e_c*E/Ma)
  Ux = 0
  Uy = U*sqrt(1. - params%IC_xip**2.)
  Uz = U*params%IC_xip
  sigma_v = vT/sqrt(2.)

  ! Box-muller:
  R_1 = sigma_v*sqrt(-2.*log(X1))
  t_2 = 2.*pi*X2
  R_3 = sigma_v*sqrt(-2.*log(X3))
  t_4 = 2.*pi*X4

  wx = R_1*cos(t_2)
  wy = R_1*cos(t_2)
  wz = R_3*cos(t_4)

  vx = Ux + wx
  vy = Uy + wy
  vz = Uz + wz
  v = sqrt( vx**2. + vy**2. + vz**2. )

  ! Populate output variables:
  plasma%kep(i) = 0.5*(Ma/e_c)*v**2
  plasma%xip(i) = vz/v

RETURN
END SUBROUTINE loadParticles
