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
REAL(r8) :: dB, dV
REAL(r8) :: Ma, qa

! Test particle mass:
Ma = params%Ma

! Test particle charge:
qa = params%qa

! Calculate the magnetic field gradient at zp0:
CALL Interp1(zp0,dB,fieldspline%dB)

! Calculate electric potential gradient at zp0:
CALL Interp1(zp0,dV,fieldspline%dV)

! Assign values to output variables:
K = upar0
L = -(1/Ma)*(mu0*dB + qa*dV)
M = 0

RETURN
END SUBROUTINE RightHandSide

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

! Record event:
!ecnt = ecnt + kep0
!pcnt = pcnt + 1

! Choose particle BC:
! 1: Isotropic plasma source
! 2: NBI
! 3: Periodic boundary
if (params%BC_Type .EQ. 1 .OR. params%BC_Type .EQ. 2) then

  if (params%BC_Type .EQ. 1) then
    T = params%Te0
    E = 0
  else if (params%BC_Type .EQ. 2) then
    T = params%BC_Tp
    E = params%BC_Ep
  end if

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

else if (params%BC_Type .EQ. 3) then
  if (zp0 .GE. params%zmax) then
    plasma%zp(i) = params%zmin
  else
    plasma%zp(i) = params%zmax
  end if
end if

RETURN
END SUBROUTINE ReinjectParticles

! =======================================================================================================
SUBROUTINE CyclotronResonanceNumber(i,plasma,resNum0,fieldspline,params)
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
REAL(r8) :: dB, ddB, dV
REAL(r8) :: Bf, Omega, dOmega, ddOmega
REAL(r8) :: Omega_dot, Omega_ddot, tau_rf
REAL(r8) :: rl, flr, besselterm
REAL(r8) :: mean_dkep_per
!REAL(r8) :: dkep_par, dkep, kep1
!REAL(r8) :: kep_per1, kep_par1
!REAL(r8) :: upar1, u1
REAL(r8) :: Ma, qa

! Input variables:
zp0  = plasma%zp(i)
kep0 = plasma%kep(i)
xip0 = plasma%xip(i)

! Test particle mass:
Ma = params%Ma

! Test particle charge:
qa = params%qa

! Calculate derived quantities
u0       = sqrt(2.*e_c*kep0/Ma)
upar0    = u0*xip0
uper0    = u0*(1. - xip0**2)**0.5
kep_par0 = kep0*xip0**2.
kep_per0 = kep0*(1. - xip0**2.)

! Gradients:
CALL Interp1(zp0,dB ,fieldspline%dB)
CALL Interp1(zp0,ddB,fieldspline%ddB)
CALL Interp1(zp0,dV ,fieldspline%dV )

! Spatial derivatives of the magnetic field:
CALL Interp1(zp0,Bf,fieldspline%B)
Omega     = params%n_harmonic*e_c*Bf/Ma
dOmega    = params%n_harmonic*e_c*dB/Ma
ddOmega   = params%n_harmonic*e_c*ddB/Ma

! Calculate the first and second time derivative of Omega:
Omega_dot = upar0*dOmega
Omega_ddot = (upar0**2.)*ddOmega  - (uper0**2.)*dOmega*dOmega/(2.*Omega) - qa*dV*dOmega/Ma

! Calculate the interaction time (tau_RF):
if ( (Omega_ddot**2.) .GT. 4.8175*ABS(Omega_dot**3.) )  then
        ! Approximate Ai(x) ~ 0.3833
        tau_rf = (2.*pi)*(ABS(2./Omega_ddot)**(1/3.))*0.3833
else
        tau_rf = sqrt(2.*pi/ABS(Omega_dot))
end if

! Calculate Bessel term:
rl       = uper0/(abs(qa)*Bf/Ma)
flr      = params%kper*rl
besselterm = BESSEL_JN(params%n_harmonic-1,flr)

!Calculate the mean RF kick per unit electric field:
mean_dkep_per = 0.5*(e_c/Ma)*(besselterm*tau_rf)**2.

! Populate plasma structure:
plasma%Erf_hat(i) = mean_dkep_per
plasma%doppler(i) = params%kpar*(upar0)/Omega

RETURN
END SUBROUTINE RFoperatorTerms

! =======================================================================================================
SUBROUTINE RFOperator(i,Edot3_hat,plasma,fieldspline,params)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define interface variables:
REAL(r8)               , INTENT(IN)    :: Edot3_hat
INTEGER(i4)            , INTENT(IN)    :: i
TYPE(plasmaTYP)        , INTENT(INOUT) :: plasma
TYPE(paramsTYP)        , INTENT(IN)    :: params
TYPE(fieldSplineTYP)   , INTENT(IN)    :: fieldspline

! Define local variables:
REAL(r8) :: kep0, xip0, Ew2
!REAL(r8) :: u0, upar0, uper0
REAL(r8) :: kep_par0, kep_per0
!REAL(r8) :: dB, ddB, dV
!REAL(r8) :: Bf, Omega, dOmega, ddOmega
!REAL(r8) :: Omega_dot, Omega_ddot, tau_rf
!REAL(r8) :: rl, flr, besselterm
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
Ew2 = params%Ew**2.
mean_dkep_per = plasma%Erf_hat(i)*Ew2

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
plasma%E3(i) = dkep

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
  TYPE(paramsTYP)    , INTENT(IN)    :: params

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
  zmin = params%zmin !+ .01*(in0%zmax-in0%zmin)
  zmax = params%zmax !- .01*(in0%zmax-in0%zmin)
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
