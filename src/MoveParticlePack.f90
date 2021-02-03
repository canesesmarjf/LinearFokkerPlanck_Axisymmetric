! =======================================================================================================
SUBROUTINE MoveParticle(zp0,kep0,xip0,in0,spline0,spline1)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables
TYPE(splTYP) :: spline0, spline1
TYPE(inTYP)  :: in0
REAL(r8) :: zp0, kep0, xip0                  ! Position, kinetic energy and pitch of the ith particle
REAL(r8) :: zpnew, Xipnew, uparnew, upernew, munew  ! Position, kinetic energy and pitch of the ith particle
REAL(r8) :: m_t, dt

! Storage for the MoverParticle and RHS subroutines
REAL(r8) :: zp1, zp2, zp3
REAL(r8) :: K1, K2, K3, K4
REAL(r8) :: upar0, upar1, upar2, upar3
REAL(r8) :: L1, L2, L3, L4
REAL(r8) :: mu0, mu1, mu2, mu3
REAL(r8) :: M1, M2, M3, M4
REAL(r8) :: u2
REAL(r8) :: B, Phi
REAL(r8) :: Interp1

! Time step:
dt = in0%dt

! Test particle mass:
m_t = in0%m_t

! Calculate initial parallel particle speed:
upar0 = xip0*sqrt(2.*e_c*kep0/m_t)

! Calculate initial particle speed squared:
u2 = 2.*e_c*kep0/m_t

! Calculate initial magnetic moment:
B = Interp1(zp0,spline0)
mu0 = 0.5*m_t*u2*(1 - xip0*xip0)/B

! Begin assembling RK4 solution:
call RightHandSide(zp0,upar0,mu0,K1,L1,M1,in0,spline0,spline1)             ! Update values of fK, fL and fM
zp1   = zp0   + (K1*dt/2.)
upar1 = upar0 + (L1*dt/2.)
mu1   = mu0   + (M1*dt/2.)

call RightHandSide(zp1,upar1,mu1,K2,L2,M2,in0,spline0,spline1)             ! Update values of fK, fL and fM
zp2   = zp0   + (K2*dt/2.)
upar2 = upar0 + (L2*dt/2.)
mu2   = mu0   + (M2*dt/2.)

call RightHandSide(zp2,upar2,mu2,K3,L3,M3,in0,spline0,spline1)             ! Update values of fK, fL and fM
zp3   = zp0   + K3*dt
upar3 = upar0 + L3*dt
mu3   = mu0   + M3*dt

call RightHandSide(zp3,upar3,mu3,K4,L4,M4,in0,spline0,spline1)             ! Update values of fK, fL and fM
zpnew   = zp0   + ( (K1 + (2.*K2) + (2.*K3) + K4)/6. )*dt
uparnew = upar0 + ( (L1 + (2.*L2) + (2.*L3) + L4)/6. )*dt
munew   = mu0 +   ( (M1 + (2.*M2) + (2.*M3) + M4)/6. )*dt

! Calculate the magnetic field at zpnew:
B = Interp1(zpnew,spline0)

! Based on new B and new mu, calculate new uper:
upernew = sqrt(2.*munew*B/m_t)

! New u due to new upar and new uper:
u2 = uparnew**2. + upernew**2.

! New pitch angle due to new upar and new u:
Xipnew = uparnew/sqrt(u2)

zp0  = zpnew
xip0 = Xipnew
kep0 = 0.5*m_t*u2/e_c

return
END SUBROUTINE MoveParticle

! =======================================================================================================
SUBROUTINE RightHandSide(zp0,upar0,mu0,K,L,M,in0,spline0,spline1)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables:
TYPE(splTYP) :: spline0, spline1
TYPE(inTYP) :: in0
REAL(r8) :: zp0, upar0, mu0, K, L, M
REAL(r8) :: dB, dPhi
REAL(r8) :: m_t, q_t
REAL(r8) :: diff1

! Test particle mass:
m_t = in0%m_t

! Test particle charge:
q_t = in0%q_t

! Calculate the magnetic field gradient at zp0:
dB = diff1(zp0,spline0)

! Calculate electric potential gradient at zp0:
dPhi = diff1(zp0,spline1)

! Assign values to output variables:
K = upar0
L = -(1/m_t)*(mu0*dB + q_t*dPhi)
M = 0

return
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
REAL(r8) :: m_t, T0, T, vT, sigma_v, E, U, Ux, Uy, Uz
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
  m_t = in0%m_t
  vT = sqrt(2.*e_c*T/m_t)
  U  = sqrt(2.*e_c*E/m_t)
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
  kep0 = 0.5*(m_t/e_c)*v**2
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
SUBROUTINE CyclotronResonanceNumber(zp0,kep0,xip0,f0,in0,spline0)
! =======================================================================================================

USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables
REAL(r8) :: zp0, kep0, xip0, f0     ! Input variables
REAL(r8) :: upar, Bf, Omega, Omega_RF
REAL(r8) :: m_t, q_t
TYPE(inTYP)  :: in0
TYPE(splTYP) :: spline0
REAL(r8) :: Interp1

! Test particle mass:
m_t = in0%m_t
! Test particle charge:
q_t = in0%q_t
! Parallel velocity of test particle:
upar = sqrt(2.*e_c*kep0/m_t)*xip0
! Magnetic field at location zp0 of test particle:
Bf = Interp1(zp0,spline0)
! Cyclotron frequency of test particle:
Omega = abs(q_t)*Bf/m_t
! RF frequency in rad/s:
Omega_RF = 2*pi*in0%f_RF
! Cyclotron resonance number:
f0 = Omega_RF - in0%kpar*upar - in0%n_harmonic*Omega

return
END SUBROUTINE CyclotronResonanceNumber

! =======================================================================================================
SUBROUTINE RFHeatingOperator(zp0,kep0,xip0,ecnt,pcnt,in0,spline0,spline1,spline2)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables:
TYPE(inTYP)  :: in0
TYPE(splTYP) :: spline0, spline1, spline2
REAL(r8) :: zp0, kep0, xip0, ecnt, pcnt
REAL(r8) :: u0, upar0, uper0
REAL(r8) :: kep_par0, kep_per0
REAL(r8) :: dB, ddB, dPhi
REAL(r8) :: Bf, Omega, dOmega, ddOmega
REAL(r8) :: Omega_dot, Omega_ddot, tau_rf
REAL(r8) :: rl, flr, besselterm
REAL(r8) :: mean_dkep_per, dkep_per, Rm1
REAL(r8) :: dkep_par, dkep, kep1
REAL(r8) :: kep_per1, kep_par1
REAL(r8) :: upar1, u1
REAL(r8) :: m_t, q_t
REAL(r8) :: Interp1, diff1

! Test particle mass:
m_t = in0%m_t

! Test particle charge:
q_t = in0%q_t

! Calculate derived quantities
u0       = sqrt(2.*e_c*kep0/m_t)
upar0    = u0*xip0
uper0    = u0*(1. - xip0**2)**0.5
kep_par0 = kep0*xip0**2.
kep_per0 = kep0*(1. - xip0**2.)

! Gradients:
dB = diff1(zp0,spline0)
ddB = Interp1(zp0,spline1)
dPhi = diff1(zp0,spline2)

! Spatial derivatives of the magnetic field:
Bf        = Interp1(zp0,spline0)
Omega     = in0%n_harmonic*e_c*Bf/m_t
dOmega    = in0%n_harmonic*e_c*dB/m_t
ddOmega   = in0%n_harmonic*e_c*ddB/m_t

! Calculate the first and second time derivative of Omega:
Omega_dot = upar0*dOmega
Omega_ddot = (upar0**2.)*ddOmega  - (uper0**2.)*dOmega*dOmega/(2.*Omega) - q_t*dPhi*dOmega/m_t

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
rl       = uper0/(abs(q_t)*Bf/m_t)
flr      = in0%kper*rl
besselterm = BESSEL_JN(in0%n_harmonic-1,flr)

! Calculate the cyclotron interaction:
! Using method based on VS. Chan PoP 9,2 (2002)
! Consistent with J. Carlsson'd PhD thesis (1998)
mean_dkep_per = 0.5*(e_c/m_t)*(in0%Ew*besselterm*tau_rf)**2.

! Calculate the change in perp, parallel and total energy:
call random_number(Rm1)
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
upar1 = sqrt( (2.*e_c/m_t)*abs(kep_par1) )*dsign(1.d0,xip0)*dsign(1.d0,kep_par1)
u1    = sqrt( (2.*e_c/m_t)*kep1 )
!WRITE(*,*) "zp0", zp0
!WRITE(*,*) "xip0", xip0
xip0 = upar1/u1
!WRITE(*,*) "xip1", xip0

! Record resonance event:
pcnt = pcnt + 1
! Record energy kick:
ecnt = ecnt + dkep

return
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
  REAL(r8) :: zmin, zmax, sigma_v, m_t
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
  m_t = in0%m_t
  T = in0%Tp_init
  E = in0%Ep_init
  vT = sqrt(2.*e_c*T/m_t)
  U  = sqrt(2.*e_c*E/m_t)
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
  kep0 = 0.5*(m_t/e_c)*v**2
  xip0 = vz/v

return
END SUBROUTINE loadParticles

! =======================================================================================================
FUNCTION Interp1(xq, spline0)
! =======================================================================================================
  USE local
  USE spline_fits
  USE dataTYP

  IMPLICIT NONE
  TYPE(splTYP) :: spline0
  REAL(r8) :: xq, Interp1, curv2

  Interp1 = curv2(xq,spline0%n,spline0%x,spline0%y,spline0%yp,spline0%sigma)

END FUNCTION Interp1

! =======================================================================================================
FUNCTION diff1(xq, spline0)
! =======================================================================================================
  USE local
  USE spline_fits
  USE dataTYP

  IMPLICIT NONE
  TYPE(splTYP) :: spline0
  REAL(r8) :: xq, diff1, curvd

  diff1 = curvd(xq,spline0%n,spline0%x,spline0%y,spline0%yp,spline0%sigma)

END FUNCTION diff1
