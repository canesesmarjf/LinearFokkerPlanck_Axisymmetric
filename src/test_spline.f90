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
TYPE(splTYP) :: spline_Bz
TYPE(spltestTYP) :: spline_Test
! User-defined functions:
REAL(r8) :: Interp1, diff1
! DO loop indices:
INTEGER(i4) :: i,j
INTEGER(i4) :: nz
! Pseudo random number seed:
INTEGER(i4) :: seed_size
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: seed
! OPENMP:
INTEGER(i4) :: id, numThreads
DOUBLE PRECISION :: ostart, oend, oend_estimate
! Interface with shell:
CHARACTER*300 :: command, mpwd
INTEGER(i4) :: n_mpwd, STATUS
CHARACTER*300 :: rootDir, inputFileDir, dir0, dir1
CHARACTER*300 ::  inputFile, fileName, xpSelector
CHARACTER*300 ::  nz_str, numThreads_str

! Read data from shell:
! ==============================================================================
CALL GET_ENVIRONMENT_VARIABLE('REPO_DIR',rootDir)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE_DIR',inputFileDir)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE',inputFile)
CALL GET_ENVIRONMENT_VARIABLE('OMP_NUM_THREADS',numThreads_str)
CALL GET_ENVIRONMENT_VARIABLE('NZ_STR',nz_str)

READ(numThreads_str,*) numThreads
READ(nz_str,*) nz

! Print to terminal:
! ==============================================================================
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) 'Input file:         ', TRIM(inputFile)
WRITE(*,*) 'nz:                 ', nz
WRITE(*,*) 'Number of threads:  ', numThreads
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''


END PROGRAM
