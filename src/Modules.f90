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
! MODULE input_data_structure
! module containing definition of an object to contain all input data
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

END MODULE dataTYP

!  **************************************************************
! MODULE spline_fits
! module containing arrays used in spline fits
!  **************************************************************
MODULE spline_fits
USE local

IMPLICIT NONE
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
