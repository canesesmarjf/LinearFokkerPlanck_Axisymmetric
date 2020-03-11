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
REAL(r8), DIMENSION(:), ALLOCATABLE :: pcount1, pcount2, pcount3    ! Count the number of particles incident on (1) dump, (2) target, (3) EBW resonance
REAL(r8), DIMENSION(:), ALLOCATABLE :: ecount1, ecount2, ecount3    ! Record the total energy of particle incident on (1) dump, (2) target, (3) EBW resonance
REAL(r8), DIMENSION(:), ALLOCATABLE :: ecount4                      ! Record the total energy dissipated by collisional slowing down within a time step dt
REAL(r8), DIMENSION(:), ALLOCATABLE :: pcount4                      ! Record the total number of fast particles involved in the slowing down dissipated power within a time step dt

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
INTEGER :: species_a, species_b, CollOperType
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

!  **************************************************************
!  MODULE InitialParticleDistribution
!  **************************************************************
MODULE InitialParticleDistribution
  USE local
  USE ParticlePusher
  use PhysicalConstants

  IMPLICIT NONE
  ! Declare shared variables:

  ! Declare user-defined data types:
  TYPE particleLoadData
    INTEGER(i4) :: zp_InitType
    INTEGER(i4) :: kep_InitType
    INTEGER(i4) :: xip_InitType
    REAL(r8) :: zp_init
    REAL(r8) :: zp_init_std
    REAL(r8) :: zmin
    REAL(r8) :: zmax
    REAL(r8) :: kep_init
    REAL(r8) :: xip_init
  END TYPE particleLoadData

CONTAINS
  SUBROUTINE loadParticles(data0)
    IMPLICIT NONE
    ! Declare internal variables:
    REAL(r8) :: zmin, zmax, sigma_u_init
    TYPE(particleLoadData) :: data0
    REAL(r8), DIMENSION(:), ALLOCATABLE :: RmArray1, RmArray2, RmArray3
    REAL(r8), DIMENSION(:), ALLOCATABLE :: uperArray, uparArray, uArray

    ! Allocate memory:
    ALLOCATE(RmArray1(Nparts), RmArray2(Nparts), RmArray3(Nparts))
    ALLOCATE(uparArray(Nparts),uperArray(Nparts),uArray(Nparts))

    write(*,*) "zmin: ", data0%zmin
    write(*,*) "zmax: ", data0%zmax
    write(*,*) "zp_init", data0%zp_init
    write(*,*) "zp_init_std", data0%zp_init_std
    write(*,*) "kep_init", data0%kep_init
    write(*,*) "xip_init", data0%xip_init

    ! Particle position:
    zmin = data0%zmin + .01*(data0%zmax-data0%zmin)
    zmax = data0%zmax - .01*(data0%zmax-data0%zmin)
    if (data0%zp_InitType .EQ. 1) then
        ! Uniform load
        call random_number(RmArray1)
        zp = zmin + (zmax - zmin)*RmArray1
    elseif (data0%zp_InitType .EQ. 2) then
        ! Gaussian load
        call random_number(RmArray1)
        call random_number(RmArray2)
        zp = data0%zp_init_std*sqrt(-2.*log(RmArray1))*cos(2.*pi*RmArray2)  +  data0%zp_init
    end if

    ! Particle kinetic energy:
    call random_number(RmArray1)
    call random_number(RmArray2)
    call random_number(RmArray3)
    sigma_u_init = sqrt(e_c*data0%kep_init/m_t)
    uperArray = sigma_u_init*sqrt(-2.*log(RmArray1))
    uparArray = sigma_u_init*sqrt(-2.*log(RmArray2))*cos(2.*pi*RmArray3)
    uArray    = sqrt( uperArray**2 + uparArray**2 )

    if (data0%kep_InitType .EQ. 1) then
        ! Maxwellian EEDF:
        kep = (m_t*uArray**2.)/(2.*e_c)
    elseif (data0%kep_InitType .EQ. 2) then
        ! Beam EEDF:
        kep = data0%kep_init
    end if

    ! Particle pitch angle:
    if (data0%xip_InitType .EQ. 1) then
        ! Isotropic pitch angle
        xip = uparArray/uArray
    elseif (data0%xip_InitType .EQ. 2) then
        ! Beam pitch angle
        xip = data0%xip_init
    end if
  return
  END SUBROUTINE loadParticles
END MODULE InitialParticleDistribution
