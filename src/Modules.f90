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
USE LOCAL
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
  REAL(r8)    :: G
  ! Time steps to record:
  INTEGER(i4) :: jstart, jend, jincr
  ! Domain geometry:
  REAL(r8) :: zmin, zmax, dtheta, r1, r2, Area0
  INTEGER(i4) :: NZmesh
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
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: zp, kep,xip, a
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: wL, wC, wR
 INTEGER(i4), DIMENSION(:), ALLOCATABLE :: m, f1, f2, f3, f4
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: dE1, dE2, dE3, dE4, dE5
 REAL(r8)   , DIMENSION(:), ALLOCATABLE :: udE3, udErf, doppler
 REAL(r8) :: NR , NSP, alpha, tp
 REAL(r8) :: Eplus, Eminus
 REAL(r8) :: Ndot1, Ndot2, Ndot3, Ndot4
 REAL(r8) :: Edot1, Edot2, Edot3, Edot4, uEdot3
END TYPE plasmaTYP

! -----------------------------------------------------------------------------
TYPE fieldSplineTYP
 TYPE(splTYP) :: B, dB, ddB, V, dV
END TYPE fieldSplineTYP

! -----------------------------------------------------------------------------
TYPE outputTYP
 REAL(r8)   , DIMENSION(:,:), ALLOCATABLE :: zp, kep, xip, a
 INTEGER(i4), DIMENSION(:,:), ALLOCATABLE :: m
 REAL(r8)   , DIMENSION(:)  , ALLOCATABLE :: tp, jrng
 REAL(r8)   , DIMENSION(:)  , ALLOCATABLE :: NR, NSP, ER
 REAL(r8)   , DIMENSION(:)  , ALLOCATABLE :: Eplus, Eminus
 REAL(r8)   , DIMENSION(:)  , ALLOCATABLE :: Ndot1, Ndot2, Ndot3, Ndot4, Ndot5
 REAL(r8)   , DIMENSION(:)  , ALLOCATABLE :: Edot1, Edot2, Edot3, Edot4, Edot5
 REAL(r8)   , DIMENSION(:,:), ALLOCATABLE :: n, nU, unU, nUE, P11, P22, E, B, dB, ddB
 REAL(r8)   , DIMENSION(:,:), ALLOCATABLE :: U, Ppar, Pper, Tpar, Tper
 REAL(r8)   , DIMENSION(:)  , ALLOCATABLE :: zm

END TYPE outputTYP

! -----------------------------------------------------------------------------
TYPE meshTYP
 REAL(r8) , DIMENSION(:), ALLOCATABLE :: zm
 REAL(r8) , DIMENSION(:), ALLOCATABLE :: n, nU, nUE, P11, P22
 REAL(r8) , DIMENSION(:), ALLOCATABLE :: U, Ppar, Pper, Tpar, Tper
 REAL(r8) , DIMENSION(:), ALLOCATABLE :: E, B, dB, ddB
 REAL(r8) , DIMENSION(:), ALLOCATABLE :: unU
 REAL(r8) :: LZ, zmin, zmax, dzm
 INTEGER(i4) :: NZmesh
END TYPE meshTYP


CONTAINS
! ----------------------------------------------------------------------------
SUBROUTINE AllocateMesh(mesh,params)
   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(meshTYP)  , INTENT(INOUT) :: mesh
   TYPE(paramsTYP), INTENT(IN)    :: params

   ! Declare local variables:
   INTEGER(i4) :: NZmesh
 
   ! Allocate memory to mesh:
   NZmesh = params%NZmesh
   ALLOCATE(mesh%zm(NZmesh))

   ! Allocate memory for mesh defined quantities:
   ! Ghost cells are used:
   ALLOCATE(mesh%n(NZmesh + 4))
   ALLOCATE(mesh%nU(NZmesh + 4))
   ALLOCATE(mesh%unU(NZmesh + 4))
   ALLOCATE(mesh%P11(NZmesh + 4))
   ALLOCATE(mesh%P22(NZmesh + 4))
   ALLOCATE(mesh%nUE(NZmesh + 4))
   ALLOCATE(mesh%B(NZmesh + 4))
   ALLOCATE(mesh%dB(NZmesh + 4))
   ALLOCATE(mesh%ddB(NZmesh + 4))
   ALLOCATE(mesh%E(NZmesh + 4))
   ALLOCATE(mesh%U(NZmesh + 4))
   ALLOCATE(mesh%Ppar(NZmesh + 4))
   ALLOCATE(mesh%Pper(NZmesh + 4))
   ALLOCATE(mesh%Tpar(NZmesh + 4))
   ALLOCATE(mesh%Tper(NZmesh + 4))

END SUBROUTINE AllocateMesh

! ----------------------------------------------------------------------------
SUBROUTINE InitializeMesh(mesh,params)
   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(meshTYP)  , INTENT(INOUT) :: mesh
   TYPE(paramsTYP), INTENT(IN)    :: params

   ! Declare local variables:
   INTEGER(i4), DIMENSION(params%NZmesh) :: m
   INTEGER(i4) :: i

   ! Populate fields:
   mesh%NZmesh   = params%NZmesh
   mesh%zmin     = params%zmin
   mesh%zmax     = params%zmax

   ! Derived quantities:
   mesh%LZ  = params%zmax - params%zmin
   mesh%dzm = mesh%LZ/mesh%NZmesh
   m = (/ (i, i=1,mesh%NZmesh, 1) /)

   ! Create mesh:
   mesh%zm = (m-1)*mesh%dzm + 0.5*mesh%dzm + mesh%zmin

   ! Initialize all mesh-defined quantities:
   mesh%n    = 0.
   mesh%nU   = 0.
   mesh%unU  = 0.
   mesh%nUE  = 0.
   mesh%P11  = 0.
   mesh%P22  = 0.
   mesh%B    = 0.
   mesh%E    = 0.
   mesh%dB   = 0.
   mesh%ddB  = 0.
   mesh%U    = 0.
   mesh%Ppar = 0.
   mesh%Pper = 0.
   mesh%Tpar = 0.
   mesh%Tper = 0.

END SUBROUTINE InitializeMesh

! ----------------------------------------------------------------------------
SUBROUTINE AllocatePlasma(plasma,params)
   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(plasmaTYP), INTENT(INOUT) :: plasma
   TYPE(paramsTYP), INTENT(IN)    :: params

   ! Declare local variables:
   INTEGER(i4) :: NC

   ! NC: Number of computational particles
   NC = params%NC
 
   ! Allocate memory: for all computational particles
   ALLOCATE(plasma%zp(NC)  ,plasma%kep(NC) ,plasma%xip(NC) ,plasma%a(NC))
   ALLOCATE(plasma%m(NC)   ,plasma%wL(NC)  ,plasma%wC(NC)  ,plasma%wR(NC))
   ALLOCATE(plasma%f1(NC)  ,plasma%f2(NC)  ,plasma%f3(NC)  ,plasma%f4(NC))
   ALLOCATE(plasma%dE1(NC) ,plasma%dE2(NC) ,plasma%dE3(NC) ,plasma%dE4(NC), plasma%dE5(NC))
   ALLOCATE(plasma%udErf(NC))
   ALLOCATE(plasma%doppler(NC))
   ALLOCATE(plasma%udE3(NC))

END SUBROUTINE AllocatePlasma

! --------------------------------------------------------------------------
SUBROUTINE InitializePlasma(plasma,params)
   USE PhysicalConstants
   USE OMP_LIB
   IMPLICIT NONE
   ! Declare interface variables:
   TYPE(plasmaTYP), INTENT(INOUT) :: plasma
   TYPE(paramsTYP), INTENT(INOUT) :: params
   
   ! Declare local variables:
   INTEGER(i4) :: i

   ! Derived parameters: Reference cross sectional area
   params%Area0 = 0.5*params%dtheta*( params%r2**2. - params%r1**2.)

   ! Derived parameters: Select test species
   IF (params%species_a .EQ. 1) THEN
       params%qa = -e_c
       params%Ma = m_e
   ELSE
       params%qa = +params%Zion*e_c
       params%Ma = params%Aion*m_p
   END IF
   
   ! Initialize plasma: scalar quantities:
   plasma%tp      = 0.
   plasma%NR      = params%ne0*params%Area0*(params%zmax - params%zmin)
   plasma%NSP     = params%NC
   plasma%alpha   = plasma%NR/plasma%NSP
   plasma%Eplus   = 0.
   plasma%Eminus  = 0.
   plasma%Ndot1   = 0.
   plasma%Ndot2   = 0.
   plasma%Ndot3   = 0.
   plasma%Ndot4   = 0.
   plasma%Edot1   = 0.
   plasma%Edot2   = 0.
   plasma%Edot3   = 0.
   plasma%Edot4   = 0.
   plasma%uEdot3  = 0.

   !$OMP PARALLEL DO
   ! Initialize plasma: zp, kep, xip and a:
   DO i = 1,params%NC
     plasma%a(i)  = 1.
     CALL loadParticles(i,plasma,params)
   END DO
   !$OMP END PARALLEL DO

END SUBROUTINE InitializePlasma

! --------------------------------------------------------------------------
SUBROUTINE ResetFlags(i,plasma)
   USE LOCAL
   IMPLICIT NONE
   ! Declare interface variables:
   INTEGER(i4)    , INTENT(IN)    :: i
   TYPE(plasmaTYP), INTENT(INOUT) :: plasma

   plasma%f1(i)      = 0. 
   plasma%f2(i)      = 0. 
   plasma%f3(i)      = 0. 
   plasma%f4(i)      = 0.
   plasma%dE1(i)     = 0. 
   plasma%dE2(i)     = 0. 
   plasma%dE3(i)     = 0. 
   plasma%dE4(i)     = 0.
   plasma%dE5(i)     = 0.
   plasma%udE3(i)    = 0.
   plasma%udErf(i)   = 0.
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

END SUBROUTINE ComputeFieldSpline

!------------------------------------------------------------------------
SUBROUTINE AllocateOutput(output,params)
   IMPLICIT none
   ! Declare interface variables:
   TYPE(outputTYP), INTENT(INOUT) :: output
   TYPE(paramsTYP), INTENT(IN)    :: params
   
   ! Declare local variables:
   INTEGER(i4) :: j, jsize, NS, NZmesh

   ! Determine size of temporal snapshots to record:
   jsize = (params%jend-params%jstart+1)/params%jincr
   
   ! Allocate memory:
   ALLOCATE(output%jrng(jsize))
   ALLOCATE(output%zp(params%NC ,jsize))
   ALLOCATE(output%kep(params%NC,jsize))
   ALLOCATE(output%xip(params%NC,jsize))
   ALLOCATE(output%a(params%NC  ,jsize))
   ALLOCATE(output%m(params%NC  ,jsize))
   ALLOCATE(output%tp(jsize))

   ! Create array with the indices of the time steps to save:
   output%jrng = (/ (j, j=params%jstart, params%jend, params%jincr) /)

  ! Allocate memory: For all time steps
   NS = params%NS
   ALLOCATE(output%NR(NS)   ,output%NSP(NS)  ,output%Eplus(NS),output%Eminus(NS),output%ER(NS))
   ALLOCATE(output%Ndot1(NS),output%Ndot2(NS),output%Ndot3(NS),output%Ndot4(NS) ,output%Ndot5(NS))
   ALLOCATE(output%Edot1(NS),output%Edot2(NS),output%Edot3(NS),output%Edot4(NS) ,output%Edot5(NS))
  
  ! Allocate memory for mesh-defined quantities:
  NZmesh = params%NZmesh
  ALLOCATE(output%zm(NZmesh))
  ALLOCATE(output%n(NZmesh + 4  ,jsize))
  ALLOCATE(output%nU(NZmesh + 4 ,jsize))
  ALLOCATE(output%unU(NZmesh + 4,jsize))
  ALLOCATE(output%nUE(NZmesh + 4,jsize))
  ALLOCATE(output%P11(NZmesh + 4,jsize))
  ALLOCATE(output%P22(NZmesh + 4,jsize))
  ALLOCATE(output%B(NZmesh + 4  ,jsize))
  ALLOCATE(output%E(NZmesh +4   ,jsize))
  ALLOCATE(output%dB(NZmesh + 4 ,jsize))
  ALLOCATE(output%ddB(NZmesh + 4,jsize))
  ALLOCATE(output%U(NZmesh + 4,jsize))
  ALLOCATE(output%Ppar(NZmesh + 4,jsize))
  ALLOCATE(output%Pper(NZmesh + 4,jsize))
  ALLOCATE(output%Tpar(NZmesh + 4,jsize))
  ALLOCATE(output%Tper(NZmesh + 4,jsize))

END SUBROUTINE AllocateOutput

! =================================================================================
SUBROUTINE PrintParamsToTerminal(params,inputFile)
 IMPLICIT NONE
 ! Declare interface variables:
 TYPE(paramsTYP), INTENT(IN) :: params
 CHARACTER*300 :: inputFile

 ! Print to terminal:
 WRITE(*,*) '' 
 WRITE(*,*) '*********************************************************************'
 WRITE(*,*) 'Input file:         ', TRIM(inputFile)
 WRITE(*,*) 'fileDescriptor:     ', TRIM(params%fileDescriptor)
 WRITE(*,*) 'Number of particles:', params%NC
 WRITE(*,*) 'Number of steps:    ', params%NS
 WRITE(*,*) 'Particle BC:        ', params%BC_Type
 WRITE(*,*) 'dt [ns]:            ', params%dt*1E+9
 WRITE(*,*) 'iPush:              ', params%iPush
 WRITE(*,*) 'iDrag:              ', params%iDrag
 WRITE(*,*) 'iColl:              ', params%iColl
 WRITE(*,*) 'iHeat:              ', params%iHeat
 WRITE(*,*) 'iSave:              ', params%iSave
 WRITE(*,*) 'elevel:             ', params%elevel
 WRITE(*,*) 'zTarget [m]:        ', params%zmax
 WRITE(*,*) 'zDump [m]:          ', params%zmin
 WRITE(*,*) 'BC_zp_mean [m]:     ', params%BC_zp_mean
 WRITE(*,*) 'B field file:       ', TRIM(params%BFieldFile)
 WRITE(*,*) 'Prf [kW]:           ', params%Prf*1E-3
 WRITE(*,*) 'Te0:                ', params%Te0
 WRITE(*,*) 'ne0:                ', params%ne0
 IF (params%CollOperType .EQ. 1) WRITE(*,*) 'Boozer-Only collision operator'
 IF (params%CollOperType .EQ. 2) WRITE(*,*) 'Boozer-Kim collision operator'
 WRITE(*,*) '*********************************************************************'
 WRITE(*,*) ''

END SUBROUTINE PrintParamsToTerminal

! =================================================================================
SUBROUTINE SaveData(output,dir1)
  IMPLICIT NONE

 !Declare interface variables:
 TYPE(outputTYP), INTENT(IN) :: output
 CHARACTER*300  , INTENT(IN) :: dir1

 ! Declare local variables:
 CHARACTER*300 :: fileName
 
 WRITE(*,*) "Saving data ..."

    ! Save plasma quantities:
    ! --------------------------------------------------------------------------
    fileName = TRIM(TRIM(dir1)//'/'//'zp.out')
    OPEN(UNIT=8,FILE=fileName,FORM="unformatted",STATUS="unknown")
    WRITE(8) output%zp
    CLOSE(UNIT=8)
    fileName = trim(trim(dir1)//'/'//'kep.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%kep
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'xip.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%xip
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'a.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%a
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'tp.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%tp
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'m.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%m
    CLOSE(unit=8)

    ! Saving mesh-defined quantities:
    ! -------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'z_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%zm
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'n_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%n
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'nU_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%nU
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'unU_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%unU
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'nUE_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%nUE
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'Ppar_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Ppar
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'Pper_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Pper
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'Tpar_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Tpar
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'Tper_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Tper
    CLOSE(unit=8)    
    fileName = trim(trim(dir1)//'/'//'U_mesh.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%U
    CLOSE(unit=8)


    ! Saving pcount to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'pcount1.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Ndot1
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount2.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Ndot2
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount3.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Ndot3
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount4.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Ndot4
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount5.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Ndot5
    CLOSE(unit=8)


    ! Saving ecount to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'ecount1.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Edot1
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount2.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Edot2
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount3.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Edot3
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount4.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Edot4
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount5.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%Edot5
    CLOSE(unit=8)


    ! Save ER, NR and NSP:
    ! -------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'ER.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%ER
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'NR.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%NR
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'NSP.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%NSP
    CLOSE(unit=8)

 WRITE(*,*) "Data saving complete"

END SUBROUTINE SaveData

END MODULE dataTYP
