! =======================================================================================================
SUBROUTINE MoveParticle(zp0,kep0,xip0,in0,der0,spline0,spline1)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables
TYPE(splTYP) :: spline0, spline1
TYPE(derTYP) :: der0
TYPE(inTYP)  :: in0
REAL(r8) :: zp0, kep0, xip0                  ! Position, kinetic energy and pitch of the ith particle
REAL(r8) :: zpnew, Xipnew, uparnew, upernew, munew  ! Position, kinetic energy and pitch of the ith particle
REAL(r8) :: m_t, delta_t

! Storage for the MoverParticle and RHS subroutines
REAL(r8) :: zp1, zp2, zp3
REAL(r8) :: K1, K2, K3, K4
REAL(r8) :: upar0, upar1, upar2, upar3
REAL(r8) :: L1, L2, L3, L4
REAL(r8) :: mu0, mu1, mu2, mu3
REAL(r8) :: M1, M2, M3, M4
REAL(r8) :: u2
REAL(r8) :: B, Phi                            ! Variables to hold potential field

! Time step:
delta_t = in0%dt

! Test particle mass:
m_t = der0%m_t

! Calculate initial parallel particle speed:
upar0 = xip0*sqrt(2.*e_c*kep0/m_t)

! Calculate initial particle speed squared:
u2 = 2.*e_c*kep0/m_t

! Calculate initial magnetic moment:
B = Interp1(zp0,spline0)
mu0 = 0.5*m_t*u2*(1 - xip0*xip0)/B

! Begin assembling RK4 solution:
call RightHandSide(zp0,upar0,mu0,K1,L1,M1,der0,spline0,spline1)             ! Update values of fK, fL and fM
zp1   = zp0   + (K1*delta_t/2.)
upar1 = upar0 + (L1*delta_t/2.)
mu1   = mu0   + (M1*delta_t/2.)

call RightHandSide(zp1,upar1,mu1,K2,L2,M2,der0,spline0,spline1)             ! Update values of fK, fL and fM
zp2   = zp0   + (K2*delta_t/2.)
upar2 = upar0 + (L2*delta_t/2.)
mu2   = mu0   + (M2*delta_t/2.)

call RightHandSide(zp2,upar2,mu2,K3,L3,M3,der0,spline0,spline1)             ! Update values of fK, fL and fM
zp3   = zp0   + K3*delta_t
upar3 = upar0 + L3*delta_t
mu3   = mu0   + M3*delta_t

call RightHandSide(zp3,upar3,mu3,K4,L4,M4,der0,spline0,spline1)             ! Update values of fK, fL and fM
zpnew   = zp0   + ( (K1 + (2.*K2) + (2.*K3) + K4)/6. )*delta_t
uparnew = upar0 + ( (L1 + (2.*L2) + (2.*L3) + L4)/6. )*delta_t
munew   = mu0 +   ( (M1 + (2.*M2) + (2.*M3) + M4)/6. )*delta_t

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
SUBROUTINE RightHandSide(zp0,upar0,mu0,K,L,M,der0,spline0,spline1)
! =======================================================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables:
TYPE(splTYP) :: spline0, spline1
TYPE(derTYP) :: der0
REAL(r8) :: zp0, upar0, mu0, K, L, M
REAL(r8) :: dB, dPhi
REAL(r8) :: m_t, q_t

! Test particle mass:
m_t = der0%m_t

! Test particle charge:
q_t = der0%q

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
SUBROUTINE ReinjectParticles(zp0,kep0,xip0,in0,der0,ecnt,pcnt)
! =======================================================================================================
USE local
!USE ParticlePusher
USE PhysicalConstants
!USE plasma_params
!USE collision_data
USE dataTYP

IMPLICIT NONE
! Define local variables:
TYPE(inTYP)  :: in0
TYPE(derTYP) :: der0
REAL(r8) :: zp0, kep0, xip0, ecnt, pcnt     ! Input variables
REAL(r8) :: uper, upar, u, sigma_u0
REAL(r8), DIMENSION(6) :: Rm6               ! Variable for storing 6 random numbers
REAL(r8) :: m_t, T0

! Record event:
ecnt = ecnt + kep0
pcnt = pcnt + 1

! Test particle mass:
m_t = der0%m_t

! Background particle temperature:
T0 = in0%Te0

! Particle velocity standard deviation:
sigma_u0     = sqrt(e_c*T0/m_t)

! Re-inject particle at source with new zp, kep, xip
call random_number(Rm6)
zp0 = in0%zp_init_std*sqrt(-2.*log(Rm6(1)))*cos(2.*pi*Rm6(2)) + in0%zp_init
uper = sigma_u0*sqrt(-2.*log(Rm6(3)))
upar = sigma_u0*sqrt(-2.*log(Rm6(4)))*cos(2.*pi*Rm6(5))
u    = sqrt( uper**2 + upar**2 )
kep0 = (m_t*u**2.)/(2.*e_c)
xip0 = upar/u

return
END SUBROUTINE ReinjectParticles

! =======================================================================================================
SUBROUTINE CyclotronResonanceNumber(zp0,kep0,xip0,f0,in0,der0,spline0)
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
TYPE(derTYP) :: der0

! Test particle mass:
m_t = der0%m_t
! Test particle charge:
q_t = der0%q
! Parallel velocity of test particle:
upar = sqrt(2.*e_c*kep0/m_t)*xip0
! Magnetic field at location zp0 of test particle:
Bf = Interp1(zp0,spline0)
! Cyclotron frequenc of test particle:
Omega = abs(q_t)*Bf/m_t
! RF frequency in rad/s:
Omega_RF = 2*pi*in0%f_RF
! Cyclotron resonance number:
f0 = Omega_RF - in0%kpar*upar - in0%n_harmonic*Omega

return
END SUBROUTINE CyclotronResonanceNumber

! =======================================================================================================
SUBROUTINE RFHeatingOperator(zp0,kep0,xip0,ecnt,pcnt,in0,der0,spline0,spline1,spline2)
! =======================================================================================================
USE local
USE spline_fits
!USE ParticlePusher
USE PhysicalConstants
!USE plasma_params
!use rf_heating_data
USE dataTYP

IMPLICIT NONE
! Define local variables:
TYPE(inTYP)  :: in0
TYPE(derTYP) :: der0
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

! Test particle mass:
m_t = der0%m_t

! Test particle charge:
q_t = der0%q

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
SUBROUTINE loadParticles(zp,kep,xip,in0,der0)
! =======================================================================================================
  USE local
  !USE ParticlePusher
  use PhysicalConstants
  USE dataTYP
  IMPLICIT NONE
  ! Declare internal variables:
  TYPE(inTYP)  :: in0
  TYPE(derTYP) :: der0
  REAL(r8), DIMENSION(in0%Nparts) :: zp, kep, xip
  REAL(r8), DIMENSION(in0%Nparts) :: RmArray1, RmArray2, RmArray3
  REAL(r8), DIMENSION(in0%Nparts) :: uperArray, uparArray, uArray
  REAL(r8) :: zmin, zmax, sigma_u_init, m_t

  WRITE(*,*) 'm_t', der0%m_t

  ! Particle position:
  zmin = in0%zmin + .01*(in0%zmax-in0%zmin)
  zmax = in0%zmax - .01*(in0%zmax-in0%zmin)
  if (in0%zp_InitType .EQ. 1) then
      ! Uniform load
      call random_number(RmArray1)
      zp = zmin + (zmax - zmin)*RmArray1
  elseif (in0%zp_InitType .EQ. 2) then
      ! Gaussian load
      call random_number(RmArray1)
      call random_number(RmArray2)
      zp = in0%zp_init_std*sqrt(-2.*log(RmArray1))*cos(2.*pi*RmArray2)  +  in0%zp_init
  end if

  ! Particle kinetic energy:
  call random_number(RmArray1)
  call random_number(RmArray2)
  call random_number(RmArray3)
  m_t = der0%m_t
  sigma_u_init = sqrt(e_c*in0%kep_init/m_t)
  uperArray = sigma_u_init*sqrt(-2.*log(RmArray1))
  uparArray = sigma_u_init*sqrt(-2.*log(RmArray2))*cos(2.*pi*RmArray3)
  uArray    = sqrt( uperArray**2 + uparArray**2 )

  if (in0%kep_InitType .EQ. 1) then
      ! Maxwellian EEDF:
      kep = (m_t*uArray**2.)/(2.*e_c)
  elseif (in0%kep_InitType .EQ. 2) then
      ! Beam EEDF:
      kep = in0%kep_init
  end if

  ! Particle pitch angle:
  if (in0%xip_InitType .EQ. 1) then
      ! Isotropic pitch angle
      xip = uparArray/uArray
  elseif (in0%xip_InitType .EQ. 2) then
      ! Beam pitch angle
      xip = in0%xip_init
  end if
return
END SUBROUTINE loadParticles
