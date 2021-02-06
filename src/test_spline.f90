PROGRAM test_spline
! PURPOSE:
! ==============================================================================

! MODULES:
! ==============================================================================
USE local
USE spline_fits
USE dataTYP
USE OMP_LIB

IMPLICIT NONE

! Define local variables:
! ==============================================================================
! User-defined structures:
TYPE(splTYP) :: spline_Bz, spline_dBz
! User-defined functions:
REAL(r8) :: Interp1, diff1
REAL(r8) :: curv2, curvd
! DO loop indices:
INTEGER(i4) :: i,j
INTEGER(i4) :: nz, nq
! Pseudo random number:
REAL(r8) :: R
REAL(r8) :: dummy1, dummy2
! Data for interpolation:
REAL(R8) :: x1, x2
REAL(r8), DIMENSION(:), ALLOCATABLE :: x , Bz , dBz
REAL(r8), DIMENSION(:), ALLOCATABLE :: xq, Bzq, dBzq
! OPENMP:
INTEGER(i4) :: id, numThreads
DOUBLE PRECISION :: ostart, oend, oend_estimate
! Interface with shell:
CHARACTER*300 :: command, mpwd
INTEGER(i4) :: n_mpwd, STATUS
CHARACTER*300 :: rootDir, inputFileDir, dir0, dir1
CHARACTER*300 ::  inputFile, fileName, xpSelector
CHARACTER*300 ::  nz_str, numThreads_str,nq_str

! Read data from shell:
! ==============================================================================
CALL GET_ENVIRONMENT_VARIABLE('REPO_DIR',rootDir)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE_DIR',inputFileDir)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE',inputFile)
CALL GET_ENVIRONMENT_VARIABLE('OMP_NUM_THREADS',numThreads_str)
CALL GET_ENVIRONMENT_VARIABLE('NZ',nz_str)
CALL GET_ENVIRONMENT_VARIABLE('NQ',nq_str)

READ(numThreads_str,*) numThreads
READ(nz_str,*) nz
READ(nq_str,*) nq

! Print to terminal:
! ==============================================================================
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) 'Input file:         ', TRIM(inputFile)
WRITE(*,*) 'Input file w/ dir:  ', TRIM(inputFileDir)
WRITE(*,*) 'nz:                 ', nz
WRITE(*,*) 'Nq:                 ', nq
WRITE(*,*) 'Number of threads:  ', numThreads
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Allocate memory for interpolated data:
! =============================================================================
ALLOCATE(x(nz),Bz(nz),dBz(nz))
ALLOCATE(xq(nq),Bzq(nq),dBzq(nq))

! Allocate memory to splines:
! ==============================================================================
CALL InitSpline(spline_Bz ,nz,0._8,0._8,1,0._8)
CALL InitSpline(spline_dBz,nz,0._8,0._8,1,0._8)

! Read external file data and create spline data:
! =============================================================================
! spline_Bz:
CALL ReadSpline(spline_Bz,inputFileDir)
x1 =  MINVAL(spline_Bz%x)
x2 =  MAXVAL(spline_Bz%x)
!spline_dBz:
spline_dBz = spline_Bz
CALL diff(spline_dBz%x,spline_Bz%y,nz,spline_dBz%y)

! Based on profile data. compute spline data:
! =============================================================================
CALL ComputeSpline(spline_Bz)
CALL ComputeSpline(spline_dBz)

! Random number generator:
! ============================================================================
ostart = OMP_GET_WTIME()

! Perfom interpolation:
! ===========================================================================
!$OMP PARALLEL DO SHARED(xq,Bzq,dBzq) PRIVATE(R,dummy1,dummy2) SCHEDULE(STATIC)
do i=1,nq
do j=1,10
  ! R needs to be private to avoid race conditions:
  CALL RANDOM_NUMBER(R)

  ! Select query point:
  xq(i) = x1*(1-R) + x2*R

  ! Select between levels of abstraction:
  if (.false.) then
    ! Using fitpack functions directly:
    Bzq(i)  = curv2(xq(i),spline_Bz%n,spline_Bz%x,spline_Bz%y,spline_Bz%yp,spline_Bz%sigma)
    dBzq(i) = curvd(xq(i),spline_Bz%n,spline_Bz%x,spline_Bz%y,spline_Bz%yp,spline_Bz%sigma)
  else if (.false.) then
    ! Using user-defined functions:
    Bzq(i)  = Interp1(xq(i),spline_Bz)
    dBzq(i) = diff1(xq(i),spline_Bz)
  else if (.true.) then
    ! Using Intrinsic functions which are probably vectorized:
    Bzq(i)  = COS(xq(i))
    dBzq(i) = SIN(xq(i))
  else if (.false.) then
    ! Use another spline interpolation tool:
    CALL splint(spline_Bz%x ,spline_Bz%y ,spline_Bz%yp ,nz,xq(i),Bzq(i) )
    CALL splint(spline_dBz%x,spline_dBz%y,spline_dBz%yp,nz,xq(i),dBzq(i))
  end if
end do
end do
!$OMP END PARALLEL DO

! Calculate computational time:
! ===========================================================================
oend = OMP_GET_WTIME()
WRITE(*,*) 'Compute time: ', (oend-ostart)*1000. ,' [ms]'

! Save data:
! ==========================================================================
fileName = trim('test_spline_Bzq.out')
OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
WRITE(8) Bzq
CLOSE(unit=8)

fileName = trim('test_spline_dBzq.out')
OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
WRITE(8) dBzq
CLOSE(unit=8)

fileName = trim('test_spline_xq.out')
OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
WRITE(8) xq
CLOSE(unit=8)

fileName = trim('test_spline_Bz.out')
OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
WRITE(8) spline_Bz%y
CLOSE(unit=8)

fileName = trim('test_spline_x.out')
OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
WRITE(8) spline_Bz%x
CLOSE(unit=8)

END PROGRAM
