!  **************************************************************
! MODULE local
! module containing definitions of real and integer kinds
!  **************************************************************
MODULE local

IMPLICIT NONE
INTEGER, PARAMETER ::i4=SELECTED_INT_KIND(9),i8=SELECTED_INT_KIND(18),  &
r4=SELECTED_REAL_KIND(6,37),r8=SELECTED_REAL_KIND(13,307)

END MODULE local

!  **************************************************************
! MODULE dataTYP
! module containing definition of an object to contain all data used in
! computation
!  **************************************************************
MODULE dataTYP
USE local

IMPLICIT NONE
TYPE inTYP
  CHARACTER*150 :: fileDescriptor, rootDir
  CHARACTER*150 :: BFieldFile, BFieldFileDir
  REAL(r8) :: Ti0, Te0, ne0, dt
  REAL(r8) :: Aion, Zeff, Zion
  INTEGER(i4) :: Nparts, Nsteps, nz, species_a, species_b
  INTEGER(i4) :: jstart, jend, jincr
  INTEGER(i4) :: threads_request, threads_given, particleBC
  LOGICAL:: iDrag, iPotential, iSave, iPush, iHeat, iColl
  INTEGER(i4) :: zp_InitType,kep_InitType, xip_InitType
  REAL(r8) :: zp_init, kep_init, xip_init, zp_init_std
  REAL(r8) :: zmin, zmax
  INTEGER(i4) :: CollOperType
  REAL(r8) :: elevel
  REAL(r8) :: s1, s2, s3, phi1, phi2, phi3
  REAL(r8) :: f_RF, kpar, kper, Ew, zRes1, zRes2
  INTEGER(i4) :: n_harmonic
  REAL(r8) :: tComputeTime, tSimTime                              ! Variables to hold cpu time at start and end of computation
  REAL(r8) :: m_t, q_t
END TYPE inTYP

END MODULE dataTYP

!  **************************************************************
! MODULE spline_fits
! module containing arrays used in spline fits
!  **************************************************************
MODULE spline_fits
USE local

IMPLICIT NONE
TYPE splTYP
  INTEGER(i4) :: n, islpsw, ierr
  REAL(r8), DIMENSION(:), ALLOCATABLE :: x, y, yp, temp
  REAL(r8) :: slp1, slpn, sigma
END TYPE splTYP

TYPE spltestTYP
  INTEGER(i4) :: n
  REAL(r8), DIMENSION(:), ALLOCATABLE :: x, y1, y2, y3, y4, y5, y6
END TYPE spltestTYP

CONTAINS
  SUBROUTINE InitSpline(spline0,n,slp1,slpn,islpsw,sigma)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0
    INTEGER(i4) :: n, islpsw
    REAL(r8) :: slp1, slpn, sigma
    ALLOCATE(spline0%x(n))
    ALLOCATE(spline0%y(n))
    ALLOCATE(spline0%yp(n))
    ALLOCATE(spline0%temp(n))
    spline0%n = n
    spline0%slp1 = slp1
    spline0%slpn = slpn
    spline0%islpsw = islpsw
    spline0%sigma = sigma
  END SUBROUTINE InitSpline

  SUBROUTINE InitSplineTest(spline0,n)
    IMPLICIT NONE
    TYPE(spltestTYP) :: spline0
    INTEGER(i4) :: n
    ALLOCATE(spline0%x(n))
    ALLOCATE(spline0%y1(n))
    ALLOCATE(spline0%y2(n))
    ALLOCATE(spline0%y3(n))
    ALLOCATE(spline0%y4(n))
    ALLOCATE(spline0%y5(n))
    ALLOCATE(spline0%y6(n))
  END SUBROUTINE InitSplineTest

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

  SUBROUTINE ComputeSpline(spline0)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0

    call curv1(spline0%n,spline0%x,spline0%y,spline0%slp1,spline0%slpn, &
    spline0%islpsw,spline0%yp,spline0%temp,spline0%sigma,spline0%ierr)

  END SUBROUTINE ComputeSpline

  FUNCTION Interp1(xq, spline0)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0
    REAL(r8) :: xq, Interp1, curv2

    Interp1 = curv2(xq,spline0%n,spline0%x,spline0%y,spline0%yp,spline0%sigma)

  END FUNCTION Interp1

  FUNCTION diff1(xq, spline0)
    IMPLICIT NONE
    TYPE(splTYP) :: spline0
    REAL(r8) :: xq, diff1, curvd

    diff1 = curvd(xq,spline0%n,spline0%x,spline0%y,spline0%yp,spline0%sigma)

  END FUNCTION diff1

END MODULE spline_fits

!  **************************************************************
! MODULE PhysicalConstants
!  **************************************************************
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
