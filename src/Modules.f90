! MODULE local
! =============================================================================
! Module containing definitions of real and integer kinds
MODULE local

IMPLICIT NONE
INTEGER, PARAMETER ::i4=SELECTED_INT_KIND(9),i8=SELECTED_INT_KIND(18),  &
r4=SELECTED_REAL_KIND(6,37),r8=SELECTED_REAL_KIND(13,307)

END MODULE local

! MODULE PhysicalConstants
! =============================================================================
MODULE PhysicalConstants
USE local
IMPLICIT NONE

REAL(r8), PARAMETER :: e_0   = 8.854e-12
REAL(r8), PARAMETER :: pi    = 3.1415926
REAL(r8), PARAMETER :: e_c   = 1.602e-19
REAL(r8), PARAMETER :: m_e   = 9.109e-31
REAL(r8), PARAMETER :: m_p   = 1.672e-27
REAL(r8), PARAMETER :: c     = 299792458 ! Speed of light [m/s]

END MODULE PhysicalConstants

! MODULE spline_fits
! =============================================================================
! Module containing arrays used in spline fits
MODULE spline_fits
USE local

IMPLICIT NONE
!----------------------------------------------------------------------------
TYPE splTYP
  INTEGER(i4) :: n
  REAL(r8), DIMENSION(:), ALLOCATABLE :: x, y, y2
  REAL(r8) :: slp1, slpn
END TYPE splTYP

CONTAINS
! ----------------------------------------------------------------------------
  SUBROUTINE AllocateSpline(spline0,n,slp1,slpn)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0
    INTEGER(i4) :: n
    REAL(r8) :: slp1, slpn
    ALLOCATE(spline0%x(n))
    ALLOCATE(spline0%y(n))
    ALLOCATE(spline0%y2(n))
    spline0%n = n
    spline0%slp1 = slp1
    spline0%slpn = slpn
  END SUBROUTINE AllocateSpline

! ----------------------------------------------------------------------------
  SUBROUTINE ReadSpline(spline0,fileName)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0
    CHARACTER*250 :: fileName
    INTEGER(i4) :: i
    open(unit=8,file=fileName,status="old")
    do i=1,(spline0%n)
        read(8,*) spline0%x(i), spline0%y(i)
    end do
    close(unit=8)
  END SUBROUTINE ReadSpline

! ----------------------------------------------------------------------------
  SUBROUTINE DiffSpline(spline0,spline1)
    ! DiffSpline takes in spline0, replicates it on spline1 AND populates the
    ! y field with the derivative of spline0
    IMPLICIT NONE
    TYPE(splTYP), INTENT(IN)  :: spline0
    TYPE(splTYP), INTENT(OUT) :: spline1

    spline1 = spline0
    CALL diff(spline0%x,spline0%y,spline0%n,spline1%y)
  END SUBROUTINE DiffSpline

! ----------------------------------------------------------------------------
  SUBROUTINE ComputeSpline(spline0)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0
    CALL spline(spline0%x,spline0%y,spline0%n,spline0%slp1,spline0%slpn,spline0%y2)
  END SUBROUTINE ComputeSpline

! ----------------------------------------------------------------------------
 SUBROUTINE Interp1(xq,yq,spline0)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0
    REAL(r8) :: xq, yq
    CALL splint(spline0%x,spline0%y,spline0%y2,spline0%n,xq,yq)
 END SUBROUTINE Interp1

! ----------------------------------------------------------------------------
 SUBROUTINE diff(x,y,n,dy)
 USE local
 IMPLICIT NONE
 REAL(r8), DIMENSION(n) :: x, y, dy
 INTEGER(i4) :: n, i

 do i=1,(n-1)
  dy(i) = (y(i+1) - y(i))/(x(i+1) - x(i))
 end do

 dy(n) = dy(n-1)

 RETURN
 END SUBROUTINE diff

! ----------------------------------------------------------------------------
   SUBROUTINE spline(x, y, n, yp1, ypn, y2)
! source: http://web.gps.caltech.edu/~cdp/cloudmodel/Current/util/
!   use nrtype
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for a
! natural spline with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipiated value of n
! (adopted from Numerical Recipes in FORTRAN 77)
!
   INTEGER, PARAMETER :: DP=KIND(1.0D0)
   INTEGER:: n
   INTEGER, PARAMETER:: nmax=500
   REAL(DP):: yp1, ypn, x(n), y(n), y2(n)
   INTEGER:: i, k
   REAL(DP):: p, qn, sig, un, u(nmax)

     if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
     else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
     endif

     do i=2, n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/&
             & (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
     enddo

     if (ypn.gt..99e30) then
        qn=0.
        un=0.
     else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
     endif

     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

     do k=n-1, 1, -1
        y2(k)=y2(k)*y2(k+1)+u(k)
     enddo

     return
     END SUBROUTINE spline

! ----------------------------------------------------------------------------
   SUBROUTINE splint(xa, ya, y2a, n, x, y)
! source: http://web.gps.caltech.edu/~cdp/cloudmodel/Current/util/
!   USE nrtype
!
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function
! (with the xa(i) in order), and given the array y2a(1:n), which is the output
! from the subroutine spline, and given a value of x, this routine returns a
! cubic spline interpolated value y.
! (adopted from Numerical Recipes in FORTRAN 77)
!
   INTEGER, PARAMETER :: DP = KIND(1.0D0)
   INTEGER:: n
   REAL(DP):: x, y, xa(n), y2a(n), ya(n)
   INTEGER:: k, khi, klo
   REAL(DP):: a, b, h

     klo=1
     khi=n
1   if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if (xa(k).gt.x) then
           khi=k
        else
           klo=k
        endif
        goto 1
     endif

     h=xa(khi)-xa(klo)
     if (h.eq.0.) pause 'bad xa input in splint'

     a=(xa(khi)-x)/h
     b=(x-xa(klo))/h
     y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

     return
     END SUBROUTINE splint

END MODULE spline_fits

! MODULE dataTYP
! =============================================================================
! The following types are prototyped:
! paramsTYP: Contain simulation basic and derived parameters
! plasmaTYP: Contain arrays of simulation data
! fieldSplineTYP: Contain splines for fields
MODULE dataTYP
USE local
USE spline_fits

IMPLICIT NONE
! -----------------------------------------------------------------------------
TYPE paramsTYP
  ! Simulation name:
  CHARACTER*150 :: fileDescriptor
  ! Magnetic field input data:
  CHARACTER*150 :: repoDir, BFieldFile, BFieldFileDir
  INTEGER(i4)   :: nz
  ! Simulation conditions:
  INTEGER(i4) :: NC, NS
  REAL(r8)    :: dt
  ! Time steps to record:
  INTEGER(i4) :: jstart, jend, jincr
  ! Domain geometry:
  REAL(r8) :: zmin, zmax, dtheta, r1, r2, Area0
  ! Physics:
  LOGICAL :: iDrag, iPotential, iSave, iPush, iHeat, iColl
  ! Collision operator conditions:
  REAL(r8)    :: Ti0, Te0, ne0, elevel
  REAL(r8)    :: Aion, Zion, Ma, qa
  INTEGER(i4) :: species_a, species_b, CollOperType
  ! Particle boundary conditions:
  INTEGER(i4) :: BC_Type
  REAL(r8)    :: BC_zp_mean, BC_Ep, BC_xip, BC_zp_std, BC_Tp
  ! Initial conditions:
  INTEGER(i4) :: IC_Type
  REAL(r8)    :: IC_zp_mean, IC_Ep, IC_xip, IC_zp_std, IC_Tp
  ! RF heating operator conditions:
  REAL(r8)    :: f_RF, kpar, kper, zRes1, zRes2, Prf
  INTEGER(i4) :: n_harmonic
  ! Electric potential conditions:
  REAL(r8) :: s1, s2, s3, phi1, phi2, phi3
  ! Compute time data:
  INTEGER(i4) :: threads_given
  REAL(r8)    :: tComputeTime, tSimTime

END TYPE paramsTYP

! -----------------------------------------------------------------------------
TYPE plasmaTYP
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: zp, kep, xip, a
 INTEGER(i4), DIMENSION(:), ALLOCATABLE :: f1, f2, f3, f4
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: E1, E2, E3, E4, E3_hat
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: Erf_hat, doppler
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: NR , NSP
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: Eplus, Eminus
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: Ndot1, Ndot2, Ndot3, Ndot4
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: Edot1, Edot2, Edot3, Edot4, Edot3_hat
END TYPE plasmaTYP

! -----------------------------------------------------------------------------
TYPE fieldSplineTYP
 TYPE(splTYP) :: B, dB, ddB, V, dV, U
END TYPE fieldSplineTYP

! -----------------------------------------------------------------------------
TYPE outputTYP
 REAL(r8) , DIMENSION(:,:), ALLOCATABLE :: zp, kep, xip, a
 REAL(r8) , DIMENSION(:)  , ALLOCATABLE :: tp, jrng
 INTEGER(i4) :: jsize
END TYPE outputTYP

CONTAINS
! ----------------------------------------------------------------------------
SUBROUTINE AllocatePlasma(plasma,params)
   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(plasmaTYP), INTENT(INOUT) :: plasma
   TYPE(paramsTYP), INTENT(IN)    :: params

   ! Declare local variables:
   INTEGER(i4) :: NC, NS

   ! NC: Number of computational particles
   NC = params%NC
   ! NS: Number of time steps
   NS = params%NS

   ! Allocate memory: For all computational particles
   ALLOCATE(plasma%zp(NC) ,plasma%kep(NC) ,plasma%xip(NC) ,plasma%a(NC))
   ALLOCATE(plasma%f1(NC) ,plasma%f2(NC)  ,plasma%f3(NC)  ,plasma%f4(NC))
   ALLOCATE(plasma%E1(NC) ,plasma%E2(NC)  ,plasma%E3(NC)  ,plasma%E4(NC))
   ALLOCATE(plasma%Erf_hat(NC))
   ALLOCATE(plasma%doppler(NC))
   ALLOCATE(plasma%E3_hat(NC))

   ! Allocate memory: For all time steps
   ALLOCATE(plasma%NR(NS)   ,plasma%NSP(NS))
   ALLOCATE(plasma%Eplus(NS),plasma%Eminus(NS))
   ALLOCATE(plasma%Ndot1(NS),plasma%Ndot2(NS) ,plasma%Ndot3(NS),plasma%Ndot4(NS))
   ALLOCATE(plasma%Edot1(NS),plasma%Edot2(NS) ,plasma%Edot3(NS),plasma%Edot4(NS))
   ALLOCATE(plasma%Edot3_hat(NS))

END SUBROUTINE AllocatePlasma

! --------------------------------------------------------------------------
SUBROUTINE InitializePlasma(plasma,params)
   USE PhysicalConstants
   USE OMP_LIB
   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(plasmaTYP), INTENT(INOUT) :: plasma
   TYPE(paramsTYP), INTENT(INOUT) :: params
   INTEGER(i4) :: i, j

   ! Derived parameters: Reference cross sectional area
    params%Area0 = 0.5*params%dtheta*( params%r2**2. - params%r1**2.)

   ! Derived parameters: Select test species
   if (params%species_a .EQ. 1) then
       params%qa = -e_c
       params%Ma = m_e
   else
       params%qa = +params%Zion*e_c
       params%Ma = params%Aion*m_p
   end if

   !$OMP PARALLEL
   !$OMP DO
   ! Initialize plasma: time dependent quantities:
   DO j = 1,params%NS
     ! Initialize plasma: NR and NC
     ! When source contrained fueling is used, these quantities
     ! can be time depedent
     plasma%NR(j)         = params%ne0*params%Area0*(params%zmax - params%zmin)
     plasma%NSP(j)        = params%NC
     plasma%Eplus(j)      = 0.
     plasma%Eminus(j)     = 0.
     plasma%Ndot1(j)      = 0.
     plasma%Ndot2(j)      = 0.
     plasma%Ndot3(j)      = 0.
     plasma%Ndot4(j)      = 0.
     plasma%Edot1(j)      = 0.
     plasma%Edot2(j)      = 0.
     plasma%Edot3(j)      = 0.
     plasma%Edot4(j)      = 0.
     plasma%Edot3_hat(j)  = 0.
   END DO
   !$OMP DO

   !$OMP DO
   ! Initialize plasma: For all computational particles
   DO i = 1,params%NC
     plasma%a(i)  = 1.
     CALL loadParticles(i,plasma,params)
   END DO
   !$OMP END DO
   !$OMP END PARALLEL

END SUBROUTINE InitializePlasma

! --------------------------------------------------------------------------
SUBROUTINE ResetFlags(i,plasma)
   USE LOCAL
   IMPLICIT NONE
   ! Declare interface variables:
   INTEGER(i4)    , INTENT(IN)    :: i
   TYPE(plasmaTYP), INTENT(INOUT) :: plasma

   plasma%f1(i) = 0. ; plasma%f2(i) = 0. ; plasma%f3(i) = 0. ; plasma%f4(i) = 0.
   plasma%E1(i) = 0. ; plasma%E2(i) = 0. ; plasma%E3(i) = 0. ; plasma%E4(i) = 0.
   plasma%Erf_hat(i) = 0.
   plasma%doppler(i) = 0.;

END SUBROUTINE ResetFlags

! --------------------------------------------------------------------------
SUBROUTINE AllocateFieldSpline(fieldspline,params)
   USE spline_fits

   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(fieldSplineTYP), INTENT(INOUT) :: fieldspline
   TYPE(paramsTYP)     , INTENT(IN)    :: params

   ! Declare local variables:
   INTEGER(i4) :: NZ

   ! Size of spline data:
   NZ = params%nz

   ! Allocate memory:
   CALL AllocateSpline(fieldspline%B  ,NZ,0._8,0._8)
   CALL AllocateSpline(fieldspline%dB ,NZ,0._8,0._8)
   CALL AllocateSpline(fieldspline%ddB,NZ,0._8,0._8)
   CALL AllocateSpline(fieldspline%V  ,NZ,0._8,0._8)
   CALL AllocateSpline(fieldspline%dV ,NZ,0._8,0._8)
   CALL AllocateSpline(fieldspline%U  ,NZ,0._8,0._8)

END SUBROUTINE AllocateFieldSpline

SUBROUTINE ComputeFieldSpline(fieldspline)
   USE spline_fits

   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(fieldSplineTYP), INTENT(INOUT) :: fieldspline

   ! Complete setting up the spline data:
   CALL ComputeSpline(fieldspline%B)
   CALL ComputeSpline(fieldspline%dB)
   CALL ComputeSpline(fieldspline%ddB)
   CALL ComputeSpline(fieldspline%V)
   CALL ComputeSpline(fieldspline%dV)
   CALL ComputeSpline(fieldspline%U)

END SUBROUTINE ComputeFieldSpline

!------------------------------------------------------------------------
SUBROUTINE AllocateOutput(output,params)
   IMPLICIT none
   ! Declare interface variables:
   TYPE(outputTYP), INTENT(INOUT) :: output
   TYPE(paramsTYP), INTENT(IN)    :: params
   INTEGER(i4) :: j

   ! Determine size of temporal snapshots to record:
   output%jsize = (params%jend-params%jstart+1)/params%jincr

   ! Allocate memory:
   ALLOCATE(output%jrng(output%jsize))
   ALLOCATE(output%zp(params%NC ,output%jsize))
   ALLOCATE(output%kep(params%NC,output%jsize))
   ALLOCATE(output%xip(params%NC,output%jsize))
   ALLOCATE(output%a(params%NC  ,output%jsize))
   ALLOCATE(output%tp(output%jsize))

   ! Create array with the indices of the time steps to save:
   output%jrng = (/ (j, j=params%jstart, params%jend, params%jincr) /)

END SUBROUTINE AllocateOutput

END MODULE dataTYP
