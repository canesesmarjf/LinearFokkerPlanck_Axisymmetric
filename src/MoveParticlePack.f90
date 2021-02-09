! =======================================================================================================
SUBROUTINE MoveParticle(zp0,kep0,xip0,in0,spline_B,spline_dB,spline_dV)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE

! Define type for interface arguments 
REAL(r8), INTENT(INOUT) :: zp0, kep0, xip0
TYPE(inTYP), INTENT(IN)  :: in0
TYPE(splTYP), INTENT(IN) :: spline_B, spline_dB, spline_dV

! Local variables:
REAL(r8) :: zpnew, Xipnew, uparnew, upernew, munew  ! Position, kinetic energy and pitch of the ith particle
REAL(r8) :: Ma, dt

! Storage for the MoverParticle and RHS subroutines:
REAL(r8) :: zp1, zp2, zp3
REAL(r8) :: K1, K2, K3, K4
REAL(r8) :: upar0, upar1, upar2, upar3
REAL(r8) :: L1, L2, L3, L4
REAL(r8) :: mu0, mu1, mu2, mu3
REAL(r8) :: M1, M2, M3, M4
REAL(r8) :: u2
REAL(r8) :: B
!REAL(r8) :: curv2

! Time step:
dt = in0%dt

! Test particle mass:
Ma = in0%Ma

! Calculate initial parallel particle speed:
upar0 = xip0*sqrt(2.*e_c*kep0/Ma)

! Calculate initial particle speed squared:
u2 = 2.*e_c*kep0/Ma

! Calculate initial magnetic moment:
CALL Interp1(zp0,B,spline_B)
!B = curv2(zp0,spline0%n,spline0%x,spline0%y,spline0%yp,spline0%sigma)

mu0 = 0.5*Ma*u2*(1 - xip0*xip0)/B

! Begin assembling RK4 solution:
call RightHandSide(zp0,upar0,mu0,K1,L1,M1,in0,spline_dB,spline_dV)             ! Update values of fK, fL and fM
zp1   = zp0   + (K1*dt/2.)
upar1 = upar0 + (L1*dt/2.)
mu1   = mu0   + (M1*dt/2.)

call RightHandSide(zp1,upar1,mu1,K2,L2,M2,in0,spline_dB,spline_dV)             ! Update values of fK, fL and fM
zp2   = zp0   + (K2*dt/2.)
upar2 = upar0 + (L2*dt/2.)
mu2   = mu0   + (M2*dt/2.)

call RightHandSide(zp2,upar2,mu2,K3,L3,M3,in0,spline_dB,spline_dV)             ! Update values of fK, fL and fM
zp3   = zp0   + K3*dt
upar3 = upar0 + L3*dt
mu3   = mu0   + M3*dt

call RightHandSide(zp3,upar3,mu3,K4,L4,M4,in0,spline_dB,spline_dV)             ! Update values of fK, fL and fM
zpnew   = zp0   + ( (K1 + (2.*K2) + (2.*K3) + K4)/6. )*dt
uparnew = upar0 + ( (L1 + (2.*L2) + (2.*L3) + L4)/6. )*dt
munew   = mu0 +   ( (M1 + (2.*M2) + (2.*M3) + M4)/6. )*dt

! Calculate the magnetic field at zpnew:
CALL Interp1(zpnew,B,spline_B)
!B = curv2(zpnew,spline0%n,spline0%x,spline0%y,spline0%yp,spline0%sigma)

! Based on new B and new mu, calculate new uper:
upernew = sqrt(2.*munew*B/Ma)

! New u due to new upar and new uper:
u2 = uparnew**2. + upernew**2.

! New pitch angle due to new upar and new u:
Xipnew = uparnew/sqrt(u2)

! Output data:
zp0  = zpnew
xip0 = Xipnew
kep0 = 0.5*Ma*u2/e_c

RETURN
END SUBROUTINE MoveParticle

! =======================================================================================================
SUBROUTINE RightHandSide(zp0,upar0,mu0,K,L,M,in0,spline_dB,spline_dV)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE

! Define type of interface arguments:
REAL(r8), INTENT(IN) :: zp0, upar0, mu0
REAL(r8), INTENT(OUT) :: K, L, M
TYPE(inTYP), INTENT(IN) :: in0
TYPE(splTYP), INTENT(IN) :: spline_dB, spline_dV

! Local variables:
REAL(r8) :: dB, dV
REAL(r8) :: Ma, qa
!REAL(r8) :: curv2

! Test particle mass:
Ma = in0%Ma

! Test particle charge:
qa = in0%qa

! Calculate the magnetic field gradient at zp0:
CALL Interp1(zp0,dB,spline_dB)
!dB = curv2(zp0,spline_dB%n,spline_dB%x,spline_dB%y,spline_dB%y2,1)

! Calculate electric potential gradient at zp0:
CALL Interp1(zp0,dV,spline_dV)
!dV = curvd(zp0,spline_dV%n,spline_dV%x,spline_dV%y,spline_dV%y2,1)

! Assign values to output variables:
K = upar0
L = -(1/Ma)*(mu0*dB + qa*dV)
M = 0

RETURN
END SUBROUTINE RightHandSide

! =======================================================================================================
SUBROUTINE ReinjectParticles(zp0,kep0,xip0,in0,ecnt,pcnt)
! =======================================================================================================
USE local
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables:
TYPE(inTYP)  :: in0
REAL(r8) :: zp0, kep0, xip0, ecnt, pcnt
!REAL(r8) :: uper, upar, u, sigma_u0
REAL(r8), DIMENSION(6) :: Rm6
REAL(r8) :: Ma, T0, T, vT, sigma_v, E, U, Ux, Uy, Uz
REAL(r8) :: R_1, R_3, t_2, t_4
REAL(r8) :: wx, wy, wz, vx, vy, vz, v

! Record event:
ecnt = ecnt + kep0
pcnt = pcnt + 1

! Choose particle BC:
! 1: Isotropic plasma source
! 2: NBI
! 3: Periodic boundary
if (in0%particleBC .EQ. 1 .OR. in0%particleBC .EQ. 2) then

  if (in0%particleBC .EQ. 1) then
    T = in0%Te0
    E = 0
  else if (in0%particleBC .EQ. 2) then
    T = in0%Tp_init
    E = in0%Ep_init
  end if

  ! Velocity distribution:
  Ma = in0%Ma
  vT = sqrt(2.*e_c*T/Ma)
  U  = sqrt(2.*e_c*E/Ma)
  Ux = 0
  Uy = U*sqrt(1. - in0%xip_init**2.)
  Uz = U*in0%xip_init
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
  kep0 = 0.5*(Ma/e_c)*v**2
  xip0 = vz/v

  ! Position distribution:
  zp0 = in0%zp_init_std*sqrt(-2.*log(Rm6(5)))*cos(2.*pi*Rm6(6)) + in0%zp_init

else if (in0%particleBC .EQ. 3) then
  if (zp0 .GE. in0%zmax) then
    zp0 = in0%zmin
  else
    zp0 = in0%zmax
  end if
end if

return
END SUBROUTINE ReinjectParticles

! =======================================================================================================
SUBROUTINE CyclotronResonanceNumber(zp0,kep0,xip0,f0,in0,spline_B)
! =======================================================================================================

USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables
REAL(r8) :: zp0, kep0, xip0, f0     ! Input variables
REAL(r8) :: upar, Bf, Omega, Omega_RF
REAL(r8) :: Ma, qa
TYPE(inTYP)  :: in0
TYPE(splTYP) :: spline_B

! Test particle mass:
Ma = in0%Ma
! Test particle charge:
qa = in0%qa
! Parallel velocity of test particle:
upar = sqrt(2.*e_c*kep0/Ma)*xip0
! Magnetic field at location zp0 of test particle:
CALL Interp1(zp0,Bf,spline_B)
!Bf = Interp1(zp0,spline_B)

! Cyclotron frequency of test particle:
Omega = abs(qa)*Bf/Ma
! RF frequency in rad/s:
Omega_RF = 2*pi*in0%f_RF
! Cyclotron resonance number:
f0 = Omega_RF - in0%kpar*upar - in0%n_harmonic*Omega

return
END SUBROUTINE CyclotronResonanceNumber

! =======================================================================================================
SUBROUTINE RFHeatingOperator(zp0,kep0,xip0,ecnt,pcnt,in0,spline_B,spline_dB,spline_ddB,spline_dV)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables:
TYPE(inTYP)  :: in0
TYPE(splTYP) :: spline_B, spline_dB, spline_ddB, spline_dV
REAL(r8) :: zp0, kep0, xip0, ecnt, pcnt
REAL(r8) :: u0, upar0, uper0
REAL(r8) :: kep_par0, kep_per0
REAL(r8) :: dB, ddB, dV
REAL(r8) :: Bf, Omega, dOmega, ddOmega
REAL(r8) :: Omega_dot, Omega_ddot, tau_rf
REAL(r8) :: rl, flr, besselterm
REAL(r8) :: mean_dkep_per, dkep_per, Rm1
REAL(r8) :: dkep_par, dkep, kep1
REAL(r8) :: kep_per1, kep_par1
REAL(r8) :: upar1, u1
REAL(r8) :: Ma, qa

! Test particle mass:
Ma = in0%Ma

! Test particle charge:
qa = in0%qa

! Calculate derived quantities
u0       = sqrt(2.*e_c*kep0/Ma)
upar0    = u0*xip0
uper0    = u0*(1. - xip0**2)**0.5
kep_par0 = kep0*xip0**2.
kep_per0 = kep0*(1. - xip0**2.)

! Gradients:
CALL Interp1(zp0,dB ,spline_dB )
CALL Interp1(zp0,ddB,spline_ddB)
CALL Interp1(zp0,dV ,spline_dV )

! Spatial derivatives of the magnetic field:
CALL Interp1(zp0,Bf,spline_B)
Omega     = in0%n_harmonic*e_c*Bf/Ma
dOmega    = in0%n_harmonic*e_c*dB/Ma
ddOmega   = in0%n_harmonic*e_c*ddB/Ma

! Calculate the first and second time derivative of Omega:
Omega_dot = upar0*dOmega
Omega_ddot = (upar0**2.)*ddOmega  - (uper0**2.)*dOmega*dOmega/(2.*Omega) - qa*dV*dOmega/Ma

! Calculate the interaction time (tau_RF):
if ( (Omega_ddot**2.) .GT. 4.8175*ABS(Omega_dot**3.) )  then
        ! Approximate Ai(x) ~ 0.3833
        tau_rf = (2.*pi)*(ABS(2./Omega_ddot)**(1/3.))*0.3833
        !WRITE(*,*) 'Airy function, tau_rf', tau_rf
        !WRITE(*,*) 'zp at res', zp0
else
        tau_rf = sqrt(2.*pi/ABS(Omega_dot))
end if

! Calculate Bessel term:
rl       = uper0/(abs(qa)*Bf/Ma)
flr      = in0%kper*rl
besselterm = BESSEL_JN(in0%n_harmonic-1,flr)

! Calculate the cyclotron interaction:
! Using method based on VS. Chan PoP 9,2 (2002)
! Consistent with J. Carlsson'd PhD thesis (1998)
mean_dkep_per = 0.5*(e_c/Ma)*(in0%Ew*besselterm*tau_rf)**2.

! Calculate the change in perp, parallel and total energy:
CALL RANDOM_NUMBER(Rm1)
Rm1 = 2.*Rm1 - 1.
dkep_per = mean_dkep_per + Rm1*sqrt(2.*kep_per0*mean_dkep_per)

!WRITE(*,*) "dkep_per", dkep_per

! Given the perp kick in energy, apply the parallel energy kick
! This arises from the non-linear effect of the perturbed magnetic field
! See Stix section 10.3 "Trapped electromagnetic modes"
!dkep_par = (in0%kpar*abs(upar0)/Omega)*dkep_per
dkep_par = (in0%kpar*(upar0)/Omega)*dkep_per

! total change in energy kick
dkep = dkep_par + dkep_per

! Calculate changes in energy
! Perp degree of freedon:
kep_per1 = kep_per0 + dkep_per
! Parallel degree of freedom:
kep_par1 = kep_par0 + dkep_par

!if (kep_par1 .LT. 0) then
  !WRITE(*,*) "kep_par1", kep_par1
!end if

! Total final energy:
kep1     = kep_per1 + kep_par1

if (kep1 .LT. 0) then
  print *, 'kep0', kep0
  print *, 'dkep', dkep
  print *, 'kep1', kep1
end if

! Assign value of the new energy:
kep0 = kep1

! Calculate the new pitch angle:
upar1 = sqrt( (2.*e_c/Ma)*abs(kep_par1) )*dsign(1.d0,xip0)*dsign(1.d0,kep_par1)
u1    = sqrt( (2.*e_c/Ma)*kep1 )
!WRITE(*,*) "zp0", zp0
!WRITE(*,*) "xip0", xip0
xip0 = upar1/u1
!WRITE(*,*) "xip1", xip0

! Record resonance event:
pcnt = pcnt + 1
! Record energy kick:
ecnt = ecnt + dkep

RETURN
END SUBROUTINE RFHeatingOperator

! =======================================================================================================
SUBROUTINE loadParticles(zp0,kep0,xip0,in0)
! =======================================================================================================
  USE local
  use PhysicalConstants
  USE dataTYP

  IMPLICIT NONE
  ! Declare internal variables:
  TYPE(inTYP) :: in0
  REAL(r8) :: zp0, kep0, xip0
  REAL(r8) :: X1, X2, X3, X4
  REAL(r8) :: R_1, R_2, R_3, R_4, t_2, t_4
  REAL(r8) :: wx, wy, wz, vx, vy, vz, v
  REAL(r8) :: zmin, zmax, sigma_v, Ma
  REAL(r8) :: Ux, Uy, Uz, vT, U, T, E

  ! Particle position:
  zmin = in0%zmin + .01*(in0%zmax-in0%zmin)
  zmax = in0%zmax - .01*(in0%zmax-in0%zmin)
  if (in0%zp_InitType .EQ. 1) then
      ! Uniform load
      call random_number(X1)
      zp0 = zmin + (zmax - zmin)*X1
  elseif (in0%zp_InitType .EQ. 2) then
      ! Gaussian load
      call random_number(X1)
      call random_number(X2)
      zp0 = in0%zp_init_std*sqrt(-2.*log(X1))*cos(2.*pi*X2)  +  in0%zp_init
  end if

  ! Particle kinetic energy and pitch angle:
  ! Generate random numbers:
  call random_number(X1)
  call random_number(X2)
  call random_number(X3)
  call random_number(X4)

  ! Derived quantities:
  Ma = in0%Ma
  T = in0%Tp_init
  E = in0%Ep_init
  vT = sqrt(2.*e_c*T/Ma)
  U  = sqrt(2.*e_c*E/Ma)
  Ux = 0
  Uy = U*sqrt(1. - in0%xip_init**2.)
  Uz = U*in0%xip_init
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
  kep0 = 0.5*(Ma/e_c)*v**2
  xip0 = vz/v

RETURN
END SUBROUTINE loadParticles
