PROGRAM linFP
! ==============================================================================
! PURPOSE:
! Created by JF Caneses Marin on 2019-09-05
! Modified: 2019_10_23, Correct error in the RF energy kick (JFCM)
! Modified: 2020_03_05, OMP constructs added (JFCM)
! ==============================================================================

! MODULES:
! ==============================================================================
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP
USE OMP_LIB

! Define local variables:
! ==============================================================================
IMPLICIT NONE
! User-defined structures:
TYPE(paramsTYP)      :: params
TYPE(plasmaTYP)      :: plasma
TYPE(fieldSplineTYP) :: fieldspline
TYPE(outputTYP)      :: output
! DO loop indices:
INTEGER(i4) :: i,j,k
! Thread ID:
INTEGER(i4) :: id
! OPENMP computational time:
DOUBLE PRECISION :: ostart, oend, oend_estimate
! Cyclotron resonance numbers:
REAL(r8) :: dresNum, resNum0, resNum1
! Main simulation variables:
! simulation time:
REAl(r8) :: tp
! Simulation diagnotics:
! Count the number of particles incident on (1) dump, (2) target, (3) cyclotron resonance, (4) slowing down
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: pcount1, pcount2, pcount3, pcount4
! Record the total energy of particle incident on (1) dump, (2) target, (3) cyclotron resonance, (4) slowing down
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: ecount1, ecount2, ecount3, ecount4
! Local simulation diagnostic counters:
REAL(r8) :: ecnt1, ecnt2, ecnt3, ecnt4
REAL(r8) :: pcnt1, pcnt2, pcnt3, pcnt4
! To store system commands and fileNames:
CHARACTER*300 :: command, mpwd
INTEGER(i4) :: n_mpwd, STATUS
CHARACTER*300 :: inputFileDir, inputFile, fileName, xpSelector, repoDir, dir0, dir1
! Namelists:
namelist/params_nml/params

! Environment:
! ==============================================================================
CALL GET_ENVIRONMENT_VARIABLE('REPO_DIR',repoDir)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE',inputFile)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE_DIR',inputFileDir)

! Read input parameters into "params" structure:
! ==============================================================================
OPEN(unit=4,file=inputFileDir,status='old',form='formatted')
read(4,params_nml)
CLOSE(unit=4)

! Print to the terminal:
! ==============================================================================
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
WRITE(*,*) 'Ew:                 ', params%Ew
WRITE(*,*) 'Te0:                ', params%Te0
WRITE(*,*) 'ne0:                ', params%ne0
if (params%CollOperType .EQ. 1) WRITE(*,*) 'Boozer-Only collision operator'
if (params%CollOperType .EQ. 2) WRITE(*,*) 'Boozer-Kim collision operator'
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Allocate memory to splines:
! ===========================================================================
CALL AllocateFieldSpline(fieldspline,params)

! Allocate memory to main simulation variables:
! ==============================================================================
CALL AllocatePlasma(plasma,params)
ALLOCATE(pcount1(params%NS)  ,pcount2(params%NS)   ,pcount3(params%NS)  ,pcount4(params%NS))
ALLOCATE(ecount1(params%NS)  ,ecount2(params%NS)   ,ecount3(params%NS)  ,ecount4(params%NS))

! Allocate memory to output variables:
! ==============================================================================
CALL AllocateOutput(output,params)

! Create splines of electromagnetic fields:
! ===========================================================================
! Magnetic field data:
fileName = trim(adjustl(params%BFieldFile))
fileName = trim(adjustl(params%BFieldFileDir))//fileName
fileName = trim(adjustl(params%repoDir))//fileName

! Magnetic field B, dB and ddB:
CALL ReadSpline(fieldspline%B,fileName)
CALL diffSpline(fieldspline%B ,fieldspline%dB )
CALL diffSpline(fieldspline%dB,fieldspline%ddB)

! Electric potential V, dV:
fieldspline%V%x = fieldspline%B%x
fieldspline%V%y = 0
if (params%iPotential) then
  CALL PotentialProfile(fieldspline,params)
end if
CALL ComputeSpline(fieldspline%V)
CALL diffSpline(fieldspline%V,fieldspline%dV)

! Plasma flow U:
fieldspline%U%x = fieldspline%B%x
fieldspline%U%y = 0.

! Complete setting up the spline data :
CALL ComputeFieldSpline(fieldspline)

! Initialize simulation variables:
! ==============================================================================
CALL InitializePlasma(plasma,params)

!$OMP PARALLEL
if (OMP_GET_THREAD_NUM() .EQ. 0) params%threads_given = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) "Number of threads given: ", params%threads_given
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Initialize simulation diagnostics:
! ==============================================================================
! Record time index:
k = 1
! Time array:
tp = 0
! Leak current diagnostics:
pcount1 = 0; pcount2 = 0; pcount3 = 0; pcount4 = 0
! Energy leak diagnotics:
ecount1 = 0; ecount2 = 0; ecount3 = 0; ecount4 = 0

! Record start time:
! ==============================================================================
ostart = OMP_GET_WTIME()

! Loop over time:
! ==============================================================================
AllTime: do j = 1,params%NS

    !$OMP PARALLEL PRIVATE(pcnt1,pcnt2,pcnt3,pcnt4,ecnt1,ecnt2,ecnt3,ecnt4,dresNum,resNum0,resNum1)

    ! Initialize private particle counters:
    pcnt1 = 0; pcnt2 = 0; pcnt3 = 0; pcnt4 = 0
    ! Initialize private energy counters:
    ecnt1 = 0; ecnt2 = 0; ecnt3 = 0; ecnt4 = 0
    ! Initialize resonance number difference:
    dresNum = 0

    ! Loop over particles:
    ! ==============================================================================
    !$OMP DO SCHEDULE(STATIC)
    AllParticles: DO i = 1,params%NC

        ! Calculate Cyclotron resonance number:
        ! ------------------------------------------------------------------------
        IF (params%iHeat) CALL CyclotronResonanceNumber(i,plasma,resNum0,fieldspline,params)

        ! Push particles adiabatically:
        ! ------------------------------------------------------------------------
        IF (params%iPush) CALL MoveParticle(i,plasma,fieldspline,params)

        ! Re-inject particles:
        ! ------------------------------------------------------------------------
        IF (plasma%zp(i) .GE. params%zmax) THEN
           plasma%f2(i)  = 1
           plasma%dE2(i) = plasma%kep(i)
           CALL ReinjectParticles(i,plasma,params,ecnt2,pcnt2)
        ELSE IF (plasma%zp(i) .LE. params%zmin) THEN
           plasma%f1(i)  = 1
           plasma%dE1(i) = plasma%kep(i)
           CALL ReinjectParticles(i,plasma,params,ecnt1,pcnt1)
        END IF

        ! Apply Coulomb collision operator:
        ! ------------------------------------------------------------------------
        IF (params%iColl) CALL collisionOperator(i,plasma,ecnt4,pcnt4,params)

        ! Apply RF heating operator:
        ! ------------------------------------------------------------------------
        IF (params%iHeat) THEN
           CALL CyclotronResonanceNumber(i,plasma,resNum1,fieldspline,params)
           dresNum = dsign(1.d0,resNum0*resNum1)
           IF (dresNum .LT. 0 .AND. plasma%zp(i) .GT. params%zRes1 .AND. plasma%zp(i) .LT. params%zRes2)  THEN
              plasma%f3(i) = 1
              CALL RFHeatingOperator(i,plasma,ecnt3,pcnt3,fieldspline,params)
           END IF
        END IF

    END DO AllParticles
    !$OMP END DO

    !$OMP CRITICAL
     ecount1(j) = ecount1(j) + ecnt1
     ecount2(j) = ecount2(j) + ecnt2
     ecount3(j) = ecount3(j) + ecnt3
     ecount4(j) = ecount4(j) + ecnt4
     pcount1(j) = pcount1(j) + pcnt1
     pcount2(j) = pcount2(j) + pcnt2
     pcount3(j) = pcount3(j) + pcnt3
     pcount4(j) = pcount4(j) + pcnt4
    !$OMP END CRITICAL

    !$OMP END PARALLEL

    ! Select data to save:
    ! =====================================================================
    ! Check if data is to be saved
    if (params%iSave) then
       if (j .EQ. output%jrng(k)) then
          output%tp(k) = tp
          !$OMP PARALLEL DO
          do i = 1,params%NC
             ! Record "ith" particle at "kth" time
             output%zp(i,k)  = plasma%zp(i)
             output%kep(i,k) = plasma%kep(i)
             output%xip(i,k) = plasma%xip(i)
          end do
         !$OMP END PARALLEL DO
         k = k + 1
      end if
    end if

    ! Update time array:
    ! =========================================================================
    tp = tp + params%dt

    ! Estimate computational time:
    ! =====================================================================
    if (j .EQ. 150) then
       oend_estimate = OMP_GET_WTIME()
       WRITE(*,*) 'Estimated compute time: ', params%NS*(oend_estimate-ostart)/j,' [s]'
    end if

end do AllTime

! Record end time:
! =========================================================================
params%tSimTime = tp
oend = OMP_GET_WTIME()

params%tComputeTime = oend-ostart
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) 'Reached End of Program, Computational time [s] = ', params%tComputeTime
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Save data:
! ==============================================================================
if (params%iSave) then
    ! Create new directory to save output data:
    ! --------------------------------------------------------------------------
    dir1 = trim(params%repoDir)//'/OutputFiles'
    command = 'mkdir '//dir1
    CALL system(command,STATUS)
    WRITE(*,*) 'Status: ', STATUS

    ! Create subdirectory with input file name:
    ! -------------------------------------------------------------------------
    n_mpwd = lEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir1 = trim(params%repoDir)//'/OutputFiles/'//trim(dir0)
    command = 'mkdir '//dir1
    CALL system(command,STATUS)

    ! Create subdirectory based on "fileDescriptor" inside the input file:
    ! -------------------------------------------------------------------------
    dir1 = trim(dir1)//'/'//trim(params%fileDescriptor)
    command = 'mkdir '// dir1
    CALL system(command,STATUS)
    CALL getcwd(mpwd)

    ! Saving zp_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'zp.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%zp
    CLOSE(unit=8)

    ! Saving kep_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'kep.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%kep
    CLOSE(unit=8)

    ! Saving xip_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'xip.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%xip
    CLOSE(unit=8)

    ! Saving t_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'tp.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) output%tp
    CLOSE(unit=8)

    ! Saving pcount to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'pcount1.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) pcount1
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount2.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) pcount2
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount3.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) pcount3
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount4.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) pcount4
    CLOSE(unit=8)

    ! Saving ecount to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'ecount1.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) ecount1
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount2.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) ecount2
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount3.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) ecount3
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount4.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) ecount4
    CLOSE(unit=8)

    ! Write output metadata:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'metadata.txt')
    OPEN(unit=8,file=fileName,form="formatted",status="unknown")
    WRITE(8,NML = params_nml)
    CLOSE(unit=8)

    ! Copy input file to output directory:
    ! --------------------------------------------------------------------------
    dir0 = trim(params%repoDir)//'/InputFiles/'//trim(inputFile)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))
    CALL system(command)

    ! Copy magnetic field data:
    ! --------------------------------------------------------------------------
    dir0 = trim(params%repoDir)//'/BfieldData'//trim(params%BFieldFile)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))//'/Bfield.txt'
    CALL system(command)

    ! Create text file with commit Hash:
    ! --------------------------------------------------------------------------
    n_mpwd = lEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir0 = trim(params%repoDir)//'/OutputFiles/'//trim(dir0)//'/'//trim(params%fileDescriptor)
    fileName = trim(dir0)//'/commitHash.txt'
    command = 'git log --oneline -1 > '//trim(fileName)
    CALL system(command)

END if

END PROGRAM
