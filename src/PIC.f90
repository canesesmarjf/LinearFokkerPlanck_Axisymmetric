! =======================================================================================================
SUBROUTINE MoveParticle_RK4(i,plasma,mesh,params)
! =======================================================================================================
USE LOCAL
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE

! Define interface variables:
INTEGER(i4)    , INTENT(IN)    :: i
TYPE(plasmaTYP), INTENT(INOUT) :: plasma
TYPE(meshTYP)  , INTENT(IN)    :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Define local variables:
REAL(r8), DIMENSION(3) :: Z0, ZN, EM
REAL(r8), DIMENSION(3) :: F, dZ1, dZ2, dZ3, dZ4
REAL(r8) :: v, vz, vp
REAL(r8) :: z, K, X
REAL(r8) :: E, B, dB
REAL(r8) :: qa, Ma, dt

! Input data:
z  = plasma.zp(i)
K  = plasma.kep(i)
X  = plasma.xip(i)
B  = plasma.Bp(i)
E  = plasma.Ep(i)
dB = plasma.dBp(i)

! Parameters:
Ma = params.Ma
qa = params.qa
dt = params.dt

! Derived quantities:
v  = sqrt(2.*e_c*K/Ma)
vz = v*X
vp = v*sqrt( 1. - (X**2.) )

! Assemble RK4 inputs:
Z0 = (/z, vz, vp/)
EM = (/E,  B, dB/)

! Start RK4 solution:
! ==============================================================================
! Equations of motion are expressed as follows:
!
! dZ/dt = F, thus Z1 = Z0 + dZ
!
! where
!
! dZ = F(Z0)*dt
! Z = [z, vz, vp]
! F(1) = +vz
! F(2) = -0.5*vp*vp*dB/B + (q/Ma)*E
! F(3) = +0.5*vp*vz*dB/B

! Step 1:
ZN = Z0
CALL CalculateF(ZN,EM,params,F)
dZ1 = F*dt

! Step 2:
ZN = Z0 + dZ1/2.
CALL InterpEM(ZN,EM,mesh,params)
CALL CalculateF(ZN,EM,params,F)
dZ2 = F*dt

! Step 3:
ZN = Z0 + dZ2/2.
CALL InterpEM(ZN,EM,mesh,params)
CALL CalculateF(ZN,EM,params,F)
dZ3 = F*dt

! Step 4:
ZN = Z0 + dZ3
CALL InterpEM(ZN,EM,mesh,params)
CALL CalculateF(ZN,EM,params,F)
dZ4 = F*dt

! Assemble RK4 solution:
ZN = Z0 + (dZ1 + 2*dZ2 + 3*dZ3 + dZ4)/6.

! End of RK solution:
! ==============================================================================

! Derived quantities:
z  = ZN(1)
vz = ZN(2)
vp = ZN(3)
v = sqrt(vz**2. + vp**2.)
X = vz/v
K = (0.5*Ma/e_c)*v**2.

! Output data:
plasma.zp(i)  = z
plasma.kep(i) = K
plasma.xip(i) = X

! Notes:
! Electromagnetic fields needs to be interpolated at this new postion.
! This means that plasma%Ep, plasma%Bp, plasma%dBp and plasma%ddBp need update

RETURN
END SUBROUTINE MoveParticle_RK4

! =======================================================================================================
SUBROUTINE InterpEM(Z,EM,mesh,params)
! =======================================================================================================
USE LOCAL
USE dataTYP

IMPLICIT NONE

! Define interface variables:
REAL(r8), DIMENSION(3), INTENT(IN)    :: Z
REAL(r8), DIMENSION(3), INTENT(INOUT) :: EM
TYPE(meshTYP)         , INTENT(IN)    :: mesh
TYPE(paramsTYP)       , INTENT(IN)    :: params

! Define local variables:
REAL(r8) :: zp, zL, dz, X
INTEGER(i4) :: m, ix
REAL(r8), DIMENSION(3) :: f, W

! Assign cell:
zp = Z(1)
zL = params%zmin
dz = mesh%dzm
m = NINT(0.5 + (zp - zL)/dz)
X = mesh%zm(m) -zp

! Assignment function:
! Left:
W(1) = 0.5*(1.5 + ((X - dz)/dz) )**2.
! Center:
W(2) = 0.75 - (X/dz)**2.
! Right:
W(3) = 0.5*(1.5 - ((X + dz)/dz) )**2.

! Get nearest grid point:
ix = m + 2

! Interpolate:
! E:
f(1) = mesh%E(ix - 1)
f(2) = mesh%E(ix)
f(3) = mesh%E(ix + 1)
EM(1) = DOT_PRODUCT(W,f)

! B:
f(1) = mesh%B(ix - 1)
f(2) = mesh%B(ix)
f(3) = mesh%B(ix + 1)
EM(2) = DOT_PRODUCT(W,f)

! dB:
f(1) = mesh%dB(ix - 1)
f(2) = mesh%dB(ix)
f(3) = mesh%dB(ix + 1)
EM(3) = DOT_PRODUCT(W,f)

RETURN
END SUBROUTINE InterpEM

! =======================================================================================================
SUBROUTINE CalculateF(Z,EM,params,F)
! =======================================================================================================
USE LOCAL
USE dataTYP
USE PhysicalConstants

IMPLICIT NONE

! Define interface variables:
REAL(r8), DIMENSION(3), INTENT(IN)    :: Z
REAL(r8), DIMENSION(3), INTENT(IN)    :: EM
TYPE(paramsTYP)       , INTENT(IN)    :: params
REAL(r8), DIMENSION(3), INTENT(INOUT) :: F

! Define local variables:
REAL(r8) :: vz, vp
REAL(r8) :: E, B, dB
REAL(r8) :: Ma, qa

! Input variables:
vz = Z(2)
vp = Z(3)
E  = EM(1)
B  = EM(2)
dB = EM(3)

! Parameters:
Ma = params%Ma
qa = params%qa

! Output:
F(1) = vz
F(2) = -0.5*vp*vp*dB/B + qa*E/Ma
F(3) = +0.5*vp*vz*dB/B

RETURN
END SUBROUTINE CalculateF
