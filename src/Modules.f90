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
  REAL(r8) :: Ti0, Te0, ne0, ni0, dt
  REAL(r8) :: Aion, Zeff, Zion
  INTEGER(i4) :: Nparts, Nsteps, nz, species_a
  INTEGER(i4) :: jstart, jend, jincr
  INTEGER(i4) :: threads_request, threads_given
  LOGICAL:: iDrag, iPotential, iSave, iPush, iHeat, iColl
  INTEGER(i4) :: zp_InitType,kep_InitType, xip_InitType
  REAL(r8) :: zp_init, kep_init, xip_init, zp_init_std
  REAL(r8) :: zmin, zmax
  INTEGER(i4) :: CollOperType
  REAL(r8) :: elevel
  REAL(r8) :: s1, s2, s3, phi1, phi2, phi3
  REAL(r8) :: f_RF, kpar, kper, Ew, zRes1, zRes2
  INTEGER(i4) :: n_harmonic
END TYPE inTYP

TYPE outTYP
  REAL(r8), DIMENSION(:), ALLOCATABLE :: kep, xip, zp
  REAL(r8), DIMENSION(:), ALLOCATABLE :: pcount1, pcount2, pcount3, pcount4
  REAL(r8), DIMENSION(:), ALLOCATABLE :: ecount1, ecount2, ecount3, ecount4
END TYPE outTYP

TYPE derTYP
  REAL(r8) :: q
  REAL(r8) :: m_t
  REAL(r8) :: species_b
  REAL(r8),DIMENSION(:), ALLOCATABLE :: fcurr, fnew
  ! Need to add frame variable
END TYPE derTYP

CONTAINS
  SUBROUTINE InitOut(out0,in0)
    IMPLICIT NONE
    TYPE(inTYP)  :: in0
    TYPE(outTYP) :: out0
    INTEGER(i4)  :: n1, n2

    ! Allocate memory:
    n1 = in0%Nparts
    n2 = in0%Nsteps
    WRITE(*,*) "Nparts", in0%Nparts
    WRITE(*,*) "Nsteps", in0%Nsteps
    ALLOCATE(out0%kep(n1), out0%xip(n1), out0%zp(n1))
    ALLOCATE(out0%pcount1(n2), out0%pcount2(n2), out0%pcount3(n2), out0%pcount4(n2))
    ALLOCATE(out0%ecount1(n2), out0%ecount2(n2), out0%ecount3(n2), out0%ecount4(n2))

    ! Initialize variables
    out0%kep = 0.; out0%xip = 0.; out0%zp = 0.
    out0%pcount1 = 0.; out0%pcount2 = 0.; out0%pcount3 = 0.; out0%pcount4 = 0.
    out0%ecount1 = 0.; out0%ecount2 = 0.; out0%ecount3 = 0.; out0%ecount4 = 0.
  END SUBROUTINE InitOut

  SUBROUTINE InitDer(der0,in0)
    IMPLICIT NONE
    TYPE(inTYP)  :: in0
    TYPE(derTYP) :: der0
    INTEGER(i4)  :: n1

    ! Allocate memory:
    n1 = in0%Nparts
    ALLOCATE(der0%fcurr(n1), der0%fnew(n1))

    ! Initialize variables
    der0%fcurr = 0.; der0%fnew = 0.
  END SUBROUTINE InitDer

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

REAL(r8), DIMENSION(:), ALLOCATABLE :: z_Ref, B_Ref, Phi_Ref  ! Variables to hold reference B and Phi data
REAL(r8), DIMENSION(:), ALLOCATABLE :: ddB_Ref                ! Variables to hold 2nd spatial derivative of the reference B
REAL(r8) :: s1, s2, s3, phi1, phi2, phi3                      ! Needed for define shape of Phi
INTEGER(i4) :: nz                                             ! Variable to hold number of points in reference B and Phi data
INTEGER(i4) :: islpsw, ierr1, ierr2                           ! Variables for curv1 and curv2 from fitpack
REAL(r8) :: slp1,slpn,sigma                                   ! Variables for curv1 and curv2 from fitpack
REAL(r8), DIMENSION(:), ALLOCATABLE :: b_spl, b_temp          ! B field reference spline data
REAL(r8), DIMENSION(:), ALLOCATABLE :: ddb_spl, ddb_temp      ! B field reference spline data
REAL(r8), DIMENSION(:), ALLOCATABLE :: phi_spl, phi_temp      ! Phi field reference spline data
REAL(r8), DIMENSION(:), ALLOCATABLE :: zz, b1, ddb1           ! Variables to hold data for interpolation test
REAL(r8), DIMENSION(501) :: j0_spl, j0_temp                   ! Bessel J1 reference spline data
REAL(r8), DIMENSION(501) :: j1_spl, j1_temp                   ! Bessel J2 reference spline data
REAL(R8), DIMENSION(501) :: x_j_ref, j0_ref, j1_ref            ! Save bessel function order 1 and 2

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
! MODULE ParticlePusher
!  **************************************************************
MODULE ParticlePusher
USE local

IMPLICIT NONE
! Total number of particles
INTEGER(i4) :: Nparts
! Number of time steps
INTEGER(i4) :: Nsteps
! time step [s]
REAL(r8) :: dt
! Charge of the test particles
REAL(r8) :: q
! Mass of the test particles
REAL(r8) :: m_t

! Main simulation variables
REAL(r8), DIMENSION(:), ALLOCATABLE :: xip, zp, kep                 ! Particle position (zp), kinetic energy (KEp), pitch angle (Xip)
REAl(r8) :: tp                                                      ! Hold simulation time

! Variables to hold selected time steps
REAL(r8), DIMENSION(:,:), ALLOCATABLE :: zp_hist, kep_hist, xip_hist
REAl(r8), DIMENSION(:), ALLOCATABLE :: t_hist

END MODULE ParticlePusher




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




!  **************************************************************
!  MODULE plasma_params
!  Module containing arrays used for plasma variables
!  **************************************************************
MODULE plasma_params
USE local

IMPLICIT NONE
REAL(r8) :: Ti0, Te0, ne0, ni0     !constant density/temperature for now
REAL(r8) :: kep_init, xip_Init       ! Test particle temperature during initialization setp
REAL(r8), DIMENSION(:), ALLOCATABLE :: Ti, Te, ne, ni
INTEGER(i4) :: num_threads, id, threads_request
!INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
!INTEGER, EXTERNAL :: OMP_SET_NUM_THREADS, OMP_GET_WTIME

END MODULE plasma_params




!  **************************************************************
!  MODULE collision_data
!  **************************************************************
MODULE collision_data
USE local

IMPLICIT NONE
REAL(r8) :: Za, Zb, Ma, Mb, Aion, Zeff, Zion
INTEGER(i4) :: species_a, species_b, CollOperType
REAL(r8) :: elevel

! Switch to enable plasma drag or dynamical friction force
LOGICAL:: iDrag
! Mean location and STD of particle loading
REAL(r8) :: zp_init, zp_init_std

END MODULE collision_data


!  **************************************************************
!  MODULE rf_heating_data
!  **************************************************************
MODULE rf_heating_data
USE local

IMPLICIT NONE
REAL(r8) :: f_RF, kpar, kper, Ew
REAL(r8) :: zRes1, zRes2
INTEGER(i4) :: n_harmonic
REAL(r8),DIMENSION(:), ALLOCATABLE :: fcurr, fnew

END MODULE rf_heating_data
