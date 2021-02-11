! =======================================================================================================
SUBROUTINE collisionOperator(zp0,kep0,xip0,ecnt,pcnt,in0)
! =======================================================================================================

USE local
USE PhysicalConstants
USE dataTYP

IMPLICIT NONE

! Define interface variables:
TYPE(inTYP), INTENT(IN) :: in0
REAL(r8), INTENT(IN) :: zp0
REAL(r8), INTENT(INOUT) :: xip0, kep0, ecnt, pcnt

! Define local variables:
INTEGER(i4) :: species_a, species_b, ss
REAL(r8) :: Za, Ma, qa
REAL(r8) :: Zb, Mb, Tb, nb
REAL(r8) :: u, w, wpar, wper, v, vpar, vper
REAL(r8) :: xip_pf_0, kep_pf_0
REAL(r8) :: xip_pf_1, kep_pf_1
REAL(r8) :: xab, wTb
REAL(r8) :: nu_ab0, nu_s, nu_D, nu_E, nu_prll
REAL(r8) :: xsi1, xsi2, ran_num
REAL(r8) :: d1, d2, S1, S2, E0, E1, E2, S3, mass_term
REAL(r8) :: dE_pf, dE_lf

! Define functions:
REAL(r8) :: lnA, phi, phip, phi2p, Gb, E_nuE_d_nu_E_dE

! Test particle properties:
! -----------------------------------------------------------------------------------------------------------------
Ma = in0%Ma
qa = in0%qa
species_a = in0%species_a
if (species_a .EQ. 1) then      ! Test electrons
   Za = -1.
else if (species_a .EQ. 1) then ! Test ions
   Za = in0%Zion
end if

species_b_loop: do ss = 1,2
! Background species properties:
! -----------------------------------------------------------------------------------------------------------------
if (ss .EQ. 1) species_b = 1
if (ss .EQ. 2) species_b = 2

if(species_b .EQ. 1) then       ! Background electron
  Mb = m_e
  Tb = in0%Te0
  Zb = -1.
  nb = in0%ne0
else if(species_b .EQ. 2) then  ! Background ion
  Mb = in0%Aion*m_p
  Tb = in0%Ti0
  Zb = in0%Zion
  nb = in0%ne0
end if

! Coordinate systems:
! -----------------------------------------------------------------------------------------------------------------
! We make use of two frames: (1) laboratory frame and (2) plasma frame.
! In the non-relativistic limit, the particle velocity in the lab frame is the following:
! v_vec = u_vec + w_vec
! Where u_vec is the mean drift velocit vector and w_vec is the thermal motion relative
! to the mean drift

! Lab frame: 
! -----------------------------------------------------------------------------------------------------------------
v = sqrt(2.*e_c*kep0/Ma)
vpar = xip0*v
vper = sqrt(1 - (xip0**2.) )*v

! Plasma drift in the lab frame:
if (in0%iDrag) then
   if (zp0 .GE. in0%IC_zp_mean) then
      u = +1.*sqrt( e_c*(in0%Te0 + in0%Ti0)/(in0%Aion*m_p) )
   else
      u = -1.*sqrt( e_c*(in0%Te0 + in0%Ti0)/(in0%Aion*m_p) )
   end if
else
   u = 0.
end if

! Plasma frame:
! ----------------------------------------------------------------------------------------------------------------
! Thermal velocities in the plasma frame:
wper = vper
wpar = vpar - u
w = sqrt( (wpar**2.) + (wper**2.) )

! Initial energy and pitch angle in the plasma frame:
xip_pf_0 = wpar/w
kep_pf_0 = (0.5*Ma/e_c)*w**2.

! Diagnostic:
if (xip_pf_0**2. .GT. 1.) then
	print *,'xip_pf_0^2 > 1', xip_pf_0
end if

! Test particle to background thermal velocity ratio:
wTb = sqrt(2.*e_c*Tb/Mb)
xab = w/wTb

! Generate random numbers:
call random_number(ran_num)
xsi1 = sign(1._r8, ran_num - 0.5_r8)
call random_number(ran_num)
xsi2 = sign(1._r8, ran_num - 0.5_r8)

! Apply pitch angle scattering in plasma frame:
! ---------------------------------------------
d1 = nu_D(xab,nb,Tb,Mb,Zb,Za,Ma)*in0%dt
S1 = xip_pf_0*(1. - d1)
S2 = ( 1. - xip_pf_0*xip_pf_0 )*d1
S3 = xsi1*sqrt(S2)
xip_pf_1 = S1 + S3

! Apply energy scattering in plasma frame:
! ----------------------------------------
! Select mass term
if (in0%CollOperType .EQ. 1) mass_term = 1.                 ! Boozer-Only term
if (in0%CollOperType .EQ. 2) mass_term = 1. + (Mb/Ma)   ! Boozer-Kim term

d2 = nu_E(xab,nb,Tb,Mb,Zb,Za,Ma)*in0%dt/mass_term
E0 = (1.5 + E_nuE_d_nu_E_dE(xab))*Tb
E1 = d2*(kep_pf_0 - E0)
E2 = sqrt(Tb*kep_pf_0*d2)
dE_pf = 2.*E1
kep_pf_1 = kep_pf_0 - dE_pf + (2.*xsi2*E2)

! Diagnotics
if (isnan(xip_pf_1)) xip_pf_1 = xip_pf_0
if  (xip_pf_1**2. .GT. 1.) then
	call random_number(ran_num)
	xip_pf_1 = 2.*ran_num - 1.
end if
if (kep_pf_1 .le. 0.) kep_pf_1 = kep_pf_0

! Record energy loss due to collisions in the plasma frame:
if (kep_pf_1 .GT. in0%elevel*Tb) then
	! Record slowing down energy during time step dt
	ecnt = ecnt + dE_pf
	! Count how many particles are involved in the slowing down power calculation
	if(species_a .EQ. species_b) then
		pcnt = pcnt + 1
	end if
end if

! Based on the modifed xip_pf and kep_pf, we get:
w    = sqrt(2.*e_c*kep_pf_1/Ma)
wpar = xip_pf_1*w
wper = sqrt(1. - (xip_pf_1**2.) )*w

! Lab frame:
! ----------------------------------------------------------------------------------------------------------------
dE_lf = dE_pf
vpar  = u + wpar
vper  = wper
v     = sqrt( (vpar**2.) + (vper**2.) )

! New pitch and energy in the lab frame:
xip0 = vpar/v
kep0 = 0.5*(Ma/e_c)*v**2.

end do species_b_loop 

RETURN
END SUBROUTINE collisionOperator


! ===============================================================================================================
REAL(r8) FUNCTION phi(x)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x

phi = erf(x)

END FUNCTION phi

! ===============================================================================================================
REAL(r8) FUNCTION phip(x)
USE local
USE PhysicalConstants

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x

phip = (2./sqrt(pi))*exp(-x*x)

END FUNCTION phip

! ===============================================================================================================
REAL(r8) FUNCTION phi2p(x)
USE local
USE PhysicalConstants

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x

phi2p = -(4.*x/sqrt(pi))*exp(-x*x)

END FUNCTION phi2p

! ===============================================================================================================
REAL(r8) FUNCTION Gb(x)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x
REAL(r8) :: phi, phip

Gb = (phi(x) - x*phip(x))/(2.*x*x)

END FUNCTION Gb

! ===============================================================================================================
REAL(r8) FUNCTION lnA(nb,Tb)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: nb,Tb

lnA = 30. - log( sqrt( nb*Tb**(-3./2.) ) )

END FUNCTION lnA

! ===============================================================================================================
REAL(r8) FUNCTION nu_ab0(nb,Tb,Mb,Zb,Za,Ma)
! Fundamental collision rate:
USE local
USE PhysicalConstants

IMPLICIT NONE
REAL(r8), INTENT(IN) :: nb,Tb,Mb,Zb,Za,Ma
REAL(r8) :: lnA
REAL(r8) :: wTb

wTb = sqrt(2.*e_c*Tb/Mb)

nu_ab0 = nb*(e_c**4.)*((Za*Zb)**2.)*lnA(nb,Tb)/(2.*pi*Ma*Ma*e_0*e_0*wTb**3.)

END FUNCTION nu_ab0

! ===============================================================================================================
REAL(r8) FUNCTION E_nuE_d_nu_E_dE(x)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x
REAL(r8) :: phi, phip, phi2p

E_nuE_d_nu_E_dE = 0.5*(  ( 3.*x*phip(x) - 3.*phi(x) - (x**2.)*phi2p(x) )/( phi(x) - (x*phip(x)) )  )

END FUNCTION E_nuE_d_nu_E_dE

! ===============================================================================================================
REAL(r8) FUNCTION nu_D(x,nb,Tb,Mb,Zb,Za,Ma)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x
REAL(r8), INTENT(IN) :: nb,Tb,Mb,Zb,Za,Ma
REAL(r8) :: nu_ab0, phi, Gb

nu_D = nu_ab0(nb,Tb,Mb,Zb,Za,Ma)*(phi(x) - Gb(x))/(x**3.)

END FUNCTION nu_D

! ===============================================================================================================
REAL(r8) FUNCTION nu_s(x,nb,Tb,Mb,Zb,Za,Ma)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x
REAL(r8), INTENT(IN) :: nb,Tb,Mb,Zb,Za,Ma
REAL(r8) :: nu_ab0, Gb

nu_s = nu_ab0(nb,Tb,Mb,Zb,Za,Ma)*(1. + (Ma/Mb))*Gb(x)/x

END FUNCTION nu_s

! ===============================================================================================================
REAL(r8) FUNCTION nu_prll(x,nb,Tb,Mb,Zb,Za,Ma)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x
REAL(r8), INTENT(IN) :: nb,Tb,Mb,Zb,Za,Ma
REAL(r8) :: nu_ab0, Gb

nu_prll = nu_ab0(nb,Tb,Mb,Zb,Za,Ma)*Gb(x)/(x**3.)

END FUNCTION nu_prll

! ===============================================================================================================
REAL(r8) FUNCTION nu_E(x,nb,Tb,Mb,Zb,Za,Ma)
USE local

IMPLICIT NONE
REAL(r8), INTENT(IN) :: x
REAL(r8), INTENT(IN) :: nb,Tb,Mb,Zb,Za,Ma
REAL(r8) :: nu_ab0, Gb, phip

! Energy loss rate:
if (.false.) then
	! From Hinton 1983 EQ 92 and L. Chen 1988 EQ 50
	nu_E = nu_ab0(nb,Tb,Mb,Zb,Za,Ma)*( (2.*(Ma/Mb)*Gb(x)/x) - (phip(x)/(x**2.)) )
else
	! From L. Chen 1983 EQ 57 commonly used for NBI
	nu_E = nu_ab0(nb,Tb,Mb,Zb,Za,Ma)*(2.*(Ma/Mb))*(Gb(x)/x)
end if

END FUNCTION nu_E


