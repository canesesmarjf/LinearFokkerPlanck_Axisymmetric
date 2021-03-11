! =======================================================================================================
SUBROUTINE ApplyCollisionOperator(plasma,mesh,params)
! =======================================================================================================
USE LOCAL
USE PhysicalConstants
USE dataTYP
USE OMP_LIB

IMPLICIT NONE

! Define interface variables:
TYPE(plasmaTYP), INTENT(INOUT) :: plasma
TYPE(meshTYP)  , INTENT(IN)    :: mesh
TYPE(paramsTYP), INTENT(IN)    :: params

! Define local variables:
INTEGER(i4) :: i, ix
REAL(r8), DIMENSION (3) :: w, n, U, Tpar, Tper

! Interpolate n ,U ,Tpar and Tper to particle positions:
! =========================================
!$OMP PARALLEL DO PRIVATE(ix, w, n, U, Tpar, Tper)
DO i = 1,params%NC
	IF (plasma%f1(i) .EQ. 0 .AND. plasma%f2(i) .EQ. 0) THEN
		! Get nearest grid point:
		ix = plasma%m(i) + 2
		
		! Assignment function:
		w(1) = plasma%wL(i)
		w(2) = plasma%wC(i)
		w(3) = plasma%wR(i)
		
		! Plasma density:
		n(1) = mesh%n(ix - 1)
		n(2) = mesh%n(ix)
		n(3) = mesh%n(ix + 1)
		plasma%np(i) = w(1)*n(1) + w(2)*n(2) + w(3)*n(3)  

		! Parallel drift velocity:
		U(1) = mesh%U(ix - 1)
		U(2) = mesh%U(ix)
		U(3) = mesh%U(ix + 1)
		plasma%Up(i) = w(1)*U(1) + w(2)*U(2) + w(3)*U(3)  

		! Parallel Temperature:
		Tpar(1) = mesh%Tpar(ix - 1)
		Tpar(2) = mesh%Tpar(ix)
		Tpar(3) = mesh%Tpar(ix + 1)
		plasma%Tparp(i) = w(1)*Tpar(1) + w(2)*Tpar(2) + w(3)*Tpar(3) 

 		! Perpendicular Temperature:
		Tper(1) = mesh%Tper(ix - 1)
		Tper(2) = mesh%Tper(ix)
		Tper(3) = mesh%Tper(ix + 1)
		plasma%Tperp(i) = w(1)*Tper(1) + w(2)*Tper(2) + w(3)*Tper(3) 
	END IF
END DO
!$OMP END PARALLEL DO

! Apply collision operator:
! =========================
!$OMP PARALLEL DO
 DO i = 1,params%NC
	 CALL CollisionOperator(i,plasma,params)
 END DO
!$OMP END PARALLEL DO

END SUBROUTINE ApplyCollisionOperator

! =======================================================================================================
SUBROUTINE collisionOperator(i,plasma,params)
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
REAL(r8) :: zp0, xip0, kep0
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

! Input variables:
zp0  = plasma%zp(i)
kep0 = plasma%kep(i)
xip0 = plasma%xip(i)

! Test particle properties:
! -----------------------------------------------------------------------------------------------------------------
Ma = params%Ma
qa = params%qa
species_a = params%species_a
IF (species_a .EQ. 1) THEN      ! Test electrons
   Za = -1.
ELSE IF (species_a .EQ. 1) THEN ! Test ions
   Za = params%Zion
END IF

species_b_loop: DO ss = 1,2
! Background species properties:
! -----------------------------------------------------------------------------------------------------------------
IF (ss .EQ. 1) species_b = 1
IF (ss .EQ. 2) species_b = 2

IF (species_b .EQ. 1) THEN       ! Background electron
  Mb = m_e
  Tb = params%Te0
  Zb = -1.
  nb = params%ne0
ELSE IF (species_b .EQ. 2) THEN  ! Background ion
  Mb = params%Aion*m_p
  Tb = params%Ti0
  Zb = params%Zion
  nb = params%ne0
END IF

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
IF (params%iDrag) THEN
   IF (zp0 .GE. params%BC_zp_mean) THEN
      u = +1.*sqrt( e_c*(params%Te0 + params%Ti0)/(params%Aion*m_p) )
   ELSE
      u = -1.*sqrt( e_c*(params%Te0 + params%Ti0)/(params%Aion*m_p) )
   END IF
ELSE
   u = 0.
END IF

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
IF (xip_pf_0**2. .GT. 1.) THEN
	print *,'xip_pf_0^2 > 1', xip_pf_0
END IF

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
d1 = nu_D(xab,nb,Tb,Mb,Zb,Za,Ma)*params%dt
S1 = xip_pf_0*(1. - d1)
S2 = ( 1. - xip_pf_0*xip_pf_0 )*d1
S3 = xsi1*sqrt(S2)
xip_pf_1 = S1 + S3

! Apply energy scattering in plasma frame:
! ----------------------------------------
! Select mass term
IF (params%CollOperType .EQ. 1) mass_term = 1.                 ! Boozer-Only term
IF (params%CollOperType .EQ. 2) mass_term = 1. + (Mb/Ma)   ! Boozer-Kim term

d2 = nu_E(xab,nb,Tb,Mb,Zb,Za,Ma)*params%dt/mass_term
E0 = (1.5 + E_nuE_d_nu_E_dE(xab))*Tb
E1 = d2*(kep_pf_0 - E0)
E2 = sqrt(Tb*kep_pf_0*d2)
dE_pf = 2.*E1
kep_pf_1 = kep_pf_0 - dE_pf + (2.*xsi2*E2)

! Diagnotics
IF (isnan(xip_pf_1)) xip_pf_1 = xip_pf_0
IF  (xip_pf_1**2. .GT. 1.) THEN
	call random_number(ran_num)
	xip_pf_1 = 2.*ran_num - 1.
END IF
if (kep_pf_1 .le. 0.) kep_pf_1 = kep_pf_0

! Record energy loss due to collisions in the plasma frame:
IF (kep_pf_1 .GT. params%elevel*Tb) THEN
	! Record slowing down energy during time step dt
	!ecnt = ecnt + dE_pf
        plasma%dE4(i) = plasma%dE4(i) +  dE_pf
	! Count how many particles are involved in the slowing down power calculation
	IF (species_a .EQ. species_b) THEN
		!pcnt = pcnt + 1
                plasma%f4(i) = 1
	END IF
END IF

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

END DO species_b_loop

! Output data from calculation:
plasma%xip(i) = vpar/v
plasma%kep(i) = 0.5*(Ma/e_c)*v**2.

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
