! =======================================================================================================
SUBROUTINE collisionOperator(zp0,kep0,xip0,ecnt,pcnt)
! =======================================================================================================

USE local
USE plasma_params
USE collision_data
USE ParticlePusher
USE PhysicalConstants

IMPLICIT NONE
REAL(r8) :: xip0, kep0, zp0, ecnt, pcnt
REAL(r8) :: xip1, kep1
REAL(r8) :: u, xb
REAL(r8) :: ln_lambda, Tb
REAL(r8) :: phi, phip, phi2p, Gb
REAL(r8) :: nuab0, nu_s, nu_D, nu_E, nu_prll
REAL(r8) :: E_nuE_d_nu_E_dE
REAL(r8) :: xsi1, xsi2, ran_num
REAL(r8) :: d1, d2, S1, S2, E0, E1, E2, S3, mass_term

REAL(r8) :: upar, uper
REAL(r8) :: vpar, vper, v, Cs, vthb, vthb3
REAL(r8) :: xip_pf_0, kep_pf_0
REAL(r8) :: xip_pf_1, kep_pf_1
REAL(r8) :: dE_pf, dE_lf

if(species_b .EQ. 1) then       ! Field electron
	Mb = m_e
	Tb = Te0
	Za = 1.; Zb = 1.
else if(species_b .EQ. 2) then  ! Field-ion
	Mb = Aion*m_p
	Tb = Ti0
	Za = Zion; Zb = Zion
end if

! Particle velocities in lab frame:
! --------------------------------
! Input: kep(i), xip(i)
u = sqrt(2.*e_c*kep0/m_t)
upar = xip0*u
uper = sqrt(1 - (xip0**2.) )*u

! Plasma drift in the lab frame:
! -----------------------------
if (iDrag) then
	if (zp0 .GE. zp_init) then
		Cs = +1.*sqrt(e_c*(Te0 + Ti0)/(Aion*m_p))
	else
		Cs = -1.*sqrt(e_c*(Te0 + Ti0)/(Aion*m_p))
	end if
 else
	Cs = 0.
 end if

! Plasma frame:
! ----------------------------------------------------------------------------------------------------------------
! Convert quantities to plasma frame:
vper = uper
vpar = upar - Cs
v = sqrt( (vpar**2.) + (vper**2.) )
xip_pf_0 = vpar/v
kep_pf_0 = (0.5*m_t/e_c)*v**2.

! Define initial pitch and energy in plasma frame:
! -----------------------------------------------
!xip0 = xip_pf
!kep0 = kep_pf

! Diagnostic
if (xip_pf_0**2 .GT. 1) then
	!print *,'Start Sub: xip > 1', xip(i)
end if

! Test particle to background thermal velocity ratio:
! --------------------------------------------------
vthb = sqrt(2.*e_c*Tb/Mb)
xb = v/vthb

! Define Chandrasekhar functions:
! -------------------------------
phi = erf(xb)
phip = (2./sqrt(pi))*exp(-xb*xb)
phi2p = -(4.*xb/sqrt(pi))*exp(-xb*xb)
Gb = (phi - xb*phip)/(2.*xb*xb)

! Fundamental collision rate:
! --------------------------
vthb3 = vthb**3.
ln_lambda = 30. - log(sqrt(ne0)*Te0**(-3/2))
nuab0 = ne0*(e_c**4.)*((Za*Zb)**2.)*ln_lambda/(2.*pi*m_t*m_t*e_0*e_0*vthb3)

! Perpendicular deflection rate:
! -----------------------------
nu_D = nuab0*(phi - Gb)/(xb**3)

! Slowing down rate:
! ------------------
nu_s = nuab0*(1. + (m_t/Mb))*Gb/xb

! Parallel diffusion rate:
! ------------------------
nu_prll = nuab0*Gb/(xb**3)

! Energy loss rate:
! -----------------
if (.false.) then
	! From Hinton 1983 EQ 92 and L. Chen 1988 EQ 50
	nu_E = nuab0*( (2*(m_t/Mb)*Gb/xb) - (phip/(xb**2)) )
else
	! From L. Chen 1983 EQ 57 commonly used for NBI
	nu_E = nuab0*(2.*(m_t/Mb))*(Gb/xb)
end if

! Energy loss derivative:
! -----------------------
E_nuE_d_nu_E_dE = 0.5*(  ( 3*xb*phip - 3*phi - (xb**2)*phi2p )/( phi - (xb*phip) )  )

! Special case for electron-ion-impurity:need only pitch angle scatt.:
! -------------------------------------------------------------------
if(species_a .eq. 1 .and. species_b .eq. 2) then
	nu_D = nuab0*Zeff*(phi - Gb)/(xb**3)
end if

! Generate random numbers:
! ------------------------
call random_number(ran_num)
xsi1 = sign(1._r8, ran_num - 0.5_r8)
call random_number(ran_num)
xsi2 = sign(1._r8, ran_num - 0.5_r8)

! Apply pitch angle scattering in plasma frame:
! ---------------------------------------------
d1 = nu_D*dt
S1 = xip_pf_0*(1. - d1)
S2 = ( 1. - xip_pf_0*xip_pf_0 )*d1
S3 = xsi1*sqrt(S2)

xip_pf_1 = S1 + S3

! Diagnotics
if (isnan(xip_pf_1)) xip_pf_1 = xip_pf_0
if  (xip_pf_1**2 .GT. 1) then
	call random_number(ran_num)
	xip_pf_1 = 2.*ran_num - 1.
end if

! Apply energy scattering in plasma frame:
! ----------------------------------------
! Select mass term
if (CollOperType .EQ. 1) mass_term = 1.              ! Boozer-Only term
if (CollOperType .EQ. 2) mass_term = 1. + (Mb/m_t)   ! Boozer-Kim term

d2 = nu_E*dt/mass_term
E0 = (1.5 + E_nuE_d_nu_E_dE)*Tb
E1 = d2*(kep_pf_0 - E0)
E2 = sqrt(Tb*kep_pf_0*d2)

dE_pf = 2.*E1

kep_pf_1 = kep_pf_0 - dE_pf + (2.*xsi2*E2)

! Diagnotics
if (kep_pf_1 .le. 0.) kep_pf_1 = kep_pf_0

! Record energy loss due to collisions in the plasma frame:
! ------------------------------------------------------
if (kep_pf_1 .GT. elevel*Te0 ) then
	! Record slowing down energy during time step dt
	ecnt = ecnt + dE_pf
	! Count how many particles are involved in the slowing down power calculation
	if(species_a .EQ. species_b) then
		pcnt = pcnt + 1
	end if
end if

! End of plasma frame calculation
! ----------------------------------------------------------------------------------------------------------------

! Convert back to lab frame:
! --------------------------
dE_lf = dE_pf
v = sqrt(2.*e_c*kep_pf_1/m_t)
vpar = xip_pf_1*v
vper = sqrt( 1 - (xip_pf_1**2.) )*v
upar = vpar + Cs
uper = vper
u = sqrt( (upar**2.) + (uper**2.) )

! New pitch and energy in the lab frame:
! -------------------------------------
xip0 = upar/u
kep0 = 0.5*(m_t/e_c)*u**2.

return
END SUBROUTINE collisionOperator
