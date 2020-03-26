! =======================================================================================================
SUBROUTINE MoveParticle(zp0,kep0,xip0,spline0,spline1)
! =======================================================================================================
USE local
USE spline_fits
USE ParticlePusher
USE PhysicalConstants

IMPLICIT NONE
! Define local variables
TYPE(splTYP) :: spline0, spline1
REAL(r8) :: zp0, kep0, xip0                  ! Position, kinetic energy and pitch of the ith particle
REAL(r8) :: zpnew, Xipnew, uparnew, upernew, munew  ! Position, kinetic energy and pitch of the ith particle

! Storage for the MoverParticle and RHS subroutines
REAL(r8) :: zp1, zp2, zp3
REAL(r8) :: K1, K2, K3, K4
REAL(r8) :: upar0, upar1, upar2, upar3
REAL(r8) :: L1, L2, L3, L4
REAL(r8) :: mu0, mu1, mu2, mu3
REAL(r8) :: M1, M2, M3, M4
REAL(r8) :: u2
REAL(r8) :: B, Phi                            ! Variables to hold potential field
!REAL(r8) :: curv2, curvd

! Calculate initial parallel particle speed
upar0 = xip0*sqrt(2.*e_c*kep0/m_t)

! Calculate initial particle speed squared
u2 = 2.*e_c*kep0/m_t

! Calculate initial magnetic moment
!B    = curv2(zp0,nz,z_Ref,B_Ref,b_spl,sigma)
B = Interp1(zp0,spline0)
mu0 = 0.5*m_t*u2*(1 - xip0*xip0)/B

! Begin assembling RK4 solution:
call RightHandSide(zp0,upar0,mu0,K1,L1,M1,spline0,spline1)             ! Update values of fK, fL and fM
zp1   = zp0   + (K1*dt/2.)
upar1 = upar0 + (L1*dt/2.)
mu1   = mu0   + (M1*dt/2.)

call RightHandSide(zp1,upar1,mu1,K2,L2,M2,spline0,spline1)             ! Update values of fK, fL and fM
zp2   = zp0   + (K2*dt/2.)
upar2 = upar0 + (L2*dt/2.)
mu2   = mu0   + (M2*dt/2.)

call RightHandSide(zp2,upar2,mu2,K3,L3,M3,spline0,spline1)             ! Update values of fK, fL and fM
zp3   = zp0   + K3*dt
upar3 = upar0 + L3*dt
mu3   = mu0   + M3*dt

call RightHandSide(zp3,upar3,mu3,K4,L4,M4,spline0,spline1)             ! Update values of fK, fL and fM
zpnew   = zp0   + ( (K1 + (2.*K2) + (2.*K3) + K4)/6. )*dt
uparnew = upar0 + ( (L1 + (2.*L2) + (2.*L3) + L4)/6. )*dt
munew   = mu0 +   ( (M1 + (2.*M2) + (2.*M3) + M4)/6. )*dt

! Calculate the magnetic field at zpnew
B = Interp1(zpnew,spline0)

! Based on new B and new mu, calculate new uper
upernew = sqrt(2.*munew*B/m_t)

! New u due to new upar and new uper
u2 = uparnew**2. + upernew**2.

! New pitch angle due to new upar and new u
Xipnew = uparnew/sqrt(u2)

zp0  = zpnew
xip0 = Xipnew
kep0 = 0.5*m_t*u2/e_c

return
END SUBROUTINE MoveParticle

! =======================================================================================================
SUBROUTINE RightHandSide(zp0,upar0,mu0,K,L,M,spline0,spline1)
! =======================================================================================================
USE local
USE spline_fits
USE ParticlePusher
USE PhysicalConstants
USE plasma_params
USE collision_data

IMPLICIT NONE
! Define local variables
TYPE(splTYP) :: spline0, spline1
REAL(r8) :: zp0, upar0, mu0, K, L, M        ! Input variables
REAL(r8) :: dB, dPhi          ! Variables to hold magnetic field, gradient of magnetic and potential field
!REAL(r8) :: curvd

! Calculate the magnetic field and electric potential at zp0
!dB   = curvd(zp0,nz,z_Ref,B_Ref,b_spl,sigma)
dB = diff1(zp0,spline0)
!dPhi = curvd(zp0,nz,z_Ref,Phi_Ref,phi_spl,sigma)
dPhi = diff1(zp0,spline1)

! Assign values to output variables
K = upar0
L = -(1/m_t)*(mu0*dB + q*dPhi)
M = 0

return
END SUBROUTINE RightHandSide

! =======================================================================================================
SUBROUTINE ReinjectParticles(zp0,kep0,xip0,ecnt,pcnt)
! =======================================================================================================
USE local
USE ParticlePusher
USE PhysicalConstants
USE plasma_params
USE collision_data

IMPLICIT NONE
! Define local variables
REAL(r8) :: zp0, kep0, xip0, ecnt, pcnt     ! Input variables
REAL(r8) :: uper, upar, u, sigma_u0
REAL(r8), DIMENSION(6) :: Rm6               ! Variable for storing 6 random numbers

! Record event:
ecnt = ecnt + kep0
pcnt = pcnt + 1

! Particle velocity standard deviation:
sigma_u0     = sqrt(e_c*Ti0/m_t)

! Re-inject particle at source with new zp, kep, xip
call random_number(Rm6)
zp0 = zp_init_std*sqrt(-2.*log(Rm6(1)))*cos(2.*pi*Rm6(2)) + zp_init
uper = sigma_u0*sqrt(-2.*log(Rm6(3)))
upar = sigma_u0*sqrt(-2.*log(Rm6(4)))*cos(2.*pi*Rm6(5))
u    = sqrt( uper**2 + upar**2 )
kep0 = (m_t*u**2.)/(2.*e_c)
xip0 = upar/u

return
END SUBROUTINE ReinjectParticles

! =======================================================================================================
SUBROUTINE CyclotronResonanceNumber(zp0,kep0,xip0,f0,in0,spline0)
! =======================================================================================================

USE local
USE spline_fits
USE ParticlePusher
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE
! Define local variables
REAL(r8) :: zp0, kep0, xip0, f0     ! Input variables
REAL(r8) :: upar, Bf, Omega, Omega_RF
REAL(r8) :: curv2
TYPE(inTYP) :: in0
TYPE(splTYP) :: spline0

upar = sqrt(2.*e_c*kep0/m_t)*xip0
Bf = Interp1(zp0,spline0)
Omega = abs(q)*Bf/m_t
Omega_RF = 2*pi*in0%f_RF
f0 = Omega_RF - in0%kpar*upar - in0%n_harmonic*Omega

return
END SUBROUTINE CyclotronResonanceNumber

! =======================================================================================================
SUBROUTINE RFHeatingOperator(zp0,kep0,xip0,ecnt,pcnt,spline0,spline1,spline2)
! =======================================================================================================
USE local
USE spline_fits
USE ParticlePusher
USE PhysicalConstants
USE plasma_params
use rf_heating_data

IMPLICIT NONE
! Define local variables:
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
!Bf        = curv2(zp0,nz,z_Ref,B_Ref,b_spl,sigma)
Bf        = Interp1(zp0,spline0)
Omega     = n_harmonic*e_c*Bf/m_t
dOmega    = n_harmonic*e_c*dB/m_t
ddOmega   = n_harmonic*e_c*ddB/m_t

! Calculate the first and second time derivative of Omega:
Omega_dot = upar0*dOmega
Omega_ddot = (upar0**2.)*ddOmega  - (uper0**2.)*dOmega*dOmega/(2.*Omega) - q*dPhi*dOmega/m_t

! Calculate the interaction time (tau_RF):
if ( (Omega_ddot**2.) .GT. 4.8175*ABS(Omega_dot**3.) )  then
        ! Approximate Ai(x) ~ 0.3833
        tau_rf = (2.*pi)*(ABS(2./Omega_ddot)**(1/3.))*0.3833
else
        tau_rf = sqrt(2.*pi/ABS(Omega_dot))
end if

! Calculate Bessel term:
rl       = uper0/(abs(q)*Bf/m_t)
flr      = kper*rl
besselterm = BESSEL_JN(n_harmonic-1,flr)

! Calculate the cyclotron interaction:
! Using method based on VS. Chan PoP 9,2 (2002)
! Consistent with J. Carlsson'd PhD thesis (1998)
mean_dkep_per = 0.5*(e_c/m_t)*(Ew*besselterm*tau_rf)**2.

! Calculate the change in perp, parallel and total energy:
call random_number(Rm1)
Rm1 = 2.*Rm1 - 1.
dkep_per = mean_dkep_per + Rm1*sqrt(2.*kep_per0*mean_dkep_per)

! Given the perp kick in energy, apply the parallel energy kick
! This arises from the non-linear effect of the perturbed magnetic field
! See Stix section 10.3 "Trapped electromagnetic modes"
dkep_par = (kpar*abs(upar0)/Omega)*dkep_per

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
kep0 = kep1

! Calculate the new pitch angle:
upar1 = sqrt( (2.*e_c/m_t)*kep_par1 )
u1    = sqrt( (2.*e_c/m_t)*kep1 )
xip0 = upar1/u1

! Record resonance event:
pcnt = pcnt + 1
! Record energy kick:
ecnt = ecnt + dkep

return
END SUBROUTINE RFHeatingOperator

! =======================================================================================================
SUBROUTINE loadParticles(in0,out0,der0)
! =======================================================================================================
  USE local
  USE ParticlePusher
  use PhysicalConstants
  USE dataTYP
  IMPLICIT NONE
  ! Declare internal variables:
  TYPE(inTYP)  :: in0
  TYPE(outTYP) :: out0
  TYPE(derTYP) :: der0
  REAL(r8)     :: zmin, zmax, sigma_u_init, m_test
  REAL(r8), DIMENSION(in0%Nparts) :: RmArray1, RmArray2, RmArray3
  REAL(r8), DIMENSION(in0%Nparts) :: uperArray, uparArray, uArray

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
  m_test = der0%m_t
  sigma_u_init = sqrt(e_c*in0%kep_init/m_test)
  uperArray = sigma_u_init*sqrt(-2.*log(RmArray1))
  uparArray = sigma_u_init*sqrt(-2.*log(RmArray2))*cos(2.*pi*RmArray3)
  uArray    = sqrt( uperArray**2 + uparArray**2 )

  if (in0%kep_InitType .EQ. 1) then
      ! Maxwellian EEDF:
      kep = (m_test*uArray**2.)/(2.*e_c)
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
