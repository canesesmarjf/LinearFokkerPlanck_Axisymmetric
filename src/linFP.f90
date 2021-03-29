
PROGRAM linFP
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
TYPE(meshTYP)        :: mesh
TYPE(outputTYP)      :: output
! Accumulators:
REAL(r8) :: N1,N2,N3,N4,N5
REAL(r8) :: E1,E2,E3,E4,E5
REAL(r8) :: uE3,uN1,uN2
! Conversion factor from SP to R particles:
REAL(r8) :: alpha
! Real and super particle energy, number and change:
REAL(r8) :: ER, dNSP, dNR
! New particle weight due to fueling:
REAL(r8) :: a_new
! DO loop indices:
INTEGER(i4) :: i,j,k
! Thread ID:
INTEGER(i4) :: id
! OPENMP computational time:
DOUBLE PRECISION :: ostart, oend, oend_estimate
! Cyclotron resonance numbers:
REAL(r8) :: dresNum, resNum0, resNum1
! simulation time step:
REAl(r8) :: dt
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
OPEN(UNIT=4,FILE=inputFileDir,STATUS='old',FORM='formatted')
READ(4,params_nml)
CLOSE(UNIT=4)

! Print to the terminal:
! ==============================================================================
CALL PrintParamsToTerminal(params,inputFile)

! Allocate memory:
! ===========================================================================
CALL AllocateFieldSpline(fieldspline,params)
CALL AllocatePlasma(plasma,params)
CALL AllocateMesh(mesh,params)
CALL AllocateOutput(output,params)

! Create splines of electromagnetic fields:
! ===========================================================================
! Magnetic field data:
fileName = TRIM(ADJUSTL(params%BFieldFile))
fileName = TRIM(ADJUSTL(params%BFieldFileDir))//fileName
fileName = TRIM(ADJUSTL(params%repoDir))//fileName

! Magnetic field B, dB and ddB:
CALL ReadSpline(fieldspline%B,fileName)
CALL diffSpline(fieldspline%B ,fieldspline%dB )
CALL diffSpline(fieldspline%dB,fieldspline%ddB)

! Electric field E:
fieldspline%E%x = fieldspline%B%x
fieldspline%E%y = 0
CALL ComputeSpline(fieldspline%E)

! Complete setting up the spline data :
CALL ComputeFieldSpline(fieldspline)

! Initialize simulation variables:
! ==============================================================================
CALL InitializeMesh(mesh,params,fieldspline)
CALL InitializePlasma(plasma,mesh,params)

!$OMP PARALLEL
IF (OMP_GET_THREAD_NUM() .EQ. 0) params%threads_given = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) "Number of threads given: ", params%threads_given
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Record start time:
! ==============================================================================
ostart = OMP_GET_WTIME()

! Initialize simulation diagnostics:
! ==============================================================================
! Record time index:
k = 1
! SP to R particle conversion factor:
alpha = plasma%alpha
! Time step:
dt = params%dt

! Loop over time:
! ==============================================================================
NS_loop: DO j = 1,params%NS
    ! Reset all accumulators:
    N1  = 0.; N2  = 0.; N3  = 0.; N4 = 0.; N5 = 0.;
    E1  = 0.; E2  = 0.; E3  = 0.; E4 = 0.; E5 = 0.;
    uE3 = 0.; uN1 = 0.; uN2 = 0.; ER = 0.;

    ! Reset resonance number:
    dresNum = 0.; resNum0 = 0.; resNum1 = 0.

    ! Advance particles:
    CALL AdvanceParticles(plasma,mesh,fieldspline,params)

    ! Assign charge to cells:
    CALL AssignCell(plasma,mesh,params)

    ! Calculate moments and extrapolate to mesh:
    CALL ExtrapolateMomentsToMesh(plasma,mesh,params)

    ! Apply collision operator:
    CALL ApplyCollisionOperator(plasma,mesh,params)

    !$OMP PARALLEL DO REDUCTION(+:uN1,uN2,N1,N2,N3,N4,E1,E2,uE3,E4)
    ! Calculate particle and energy rates:
    NC_loop2: DO i = 1,params%NC
       uN1 = uN1  +             (alpha/dt)* plasma%f1(i)
       uN2 = uN2  +             (alpha/dt)* plasma%f2(i)
       N1  = N1   +  (alpha/dt)*plasma%a(i)*plasma%f1(i)
       N2  = N2   +  (alpha/dt)*plasma%a(i)*plasma%f2(i)
       N3  = N3   +  (alpha/dt)*plasma%a(i)*plasma%f3(i)
       N4  = N4   +  (alpha/dt)*plasma%a(i)*plasma%f4(i)
       E1  = E1   +  (alpha/dt)*plasma%a(i)*plasma%f1(i)*plasma%dE1(i)*e_c
       E2  = E2   +  (alpha/dt)*plasma%a(i)*plasma%f2(i)*plasma%dE2(i)*e_c
       uE3 = uE3  +  (alpha/dt)*plasma%a(i)*plasma%f3(i)*plasma%udE3(i)*e_c
       E4  = E4   +  (alpha/dt)*plasma%a(i)*plasma%f4(i)*plasma%dE4(i)*e_c
    END DO NC_loop2
    !$OMP END PARALLEL DO

    a_new = params%G/(uN1 + uN2)
    IF (a_new .GT. 100) THEN
       a_new = 1.
    END IF

    !$OMP PARALLEL DO
    ! Apply particle re-injection and RF operator:
    NC_loop3: DO i = 1,params%NC
       IF (plasma%f1(i) .EQ. 1 .OR. plasma%f2(i) .EQ. 1) THEN
           plasma%a(i) = a_new
           CALL ReinjectParticles(i,plasma,params)
           plasma%dE5(i) = plasma%kep(i)
       END IF
       IF (plasma%f3(i) .EQ. 1) THEN
           CALL RFOperator(i,plasma,fieldspline,params,uE3)
       END IF
    END DO NC_loop3
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO REDUCTION(+:E3,N5,E5,ER)
    ! Calculate RF power absorption and real particle source:
    NC_loop4: DO i = 1,params%NC
       E3 = E3 + (alpha/dt)*plasma%a(i)*plasma%f3(i)*plasma%dE3(i)*e_c
       N5 = N5 + (alpha/dt)*plasma%a(i)*(plasma%f1(i) + plasma%f2(i))
       E5 = E5 + (alpha/dt)*plasma%a(i)*(plasma%f1(i) + plasma%f2(i))*plasma%dE5(i)*e_c
       ER = ER + alpha*plasma%a(i)*plasma%kep(i)*e_c
    END DO NC_loop4
    !$OMP END PARALLEL DO

    ! Field solve:
    CALL AdvanceEfield(mesh,params)

    ! Interpolate electromagnetic fields to particle positions:
    CALL InterpolateElectromagneticFields(plasma,mesh,params)

    ! Calculate new number of real particles NR and super-particles NSP:
    ! ==============================================================================
    dNR  = a_new*(uN1 + uN2)*dt - (N1 + N2)*dt
    dNSP = dNR/alpha
    plasma%NR  = plasma%NR  + dNR
    plasma%NSP = plasma%NSP + dNSP

    ! Assign new NR and NSP:
    ! ============================================================================
    output%NR(j)  = plasma%NR
    output%NSP(j) = plasma%NSP
    output%ER(j)  = ER
    
    ! Assign particle and energy rates in physical units [P/s] and [J/s]
    ! ==============================================================================
    output%Ndot1(j) = N1
    output%Ndot2(j) = N2
    output%Ndot3(j) = N3
    output%Ndot4(j) = N4
    output%Ndot5(j) = N5
    output%Edot1(j) = E1
    output%Edot2(j) = E2
    output%Edot3(j) = E3
    output%Edot4(j) = E4
    output%Edot5(j) = E5

    ! Select data to save:
    ! =====================================================================
    ! Check if data is to be saved
    IF (params%iSave) THEN
       IF (j .EQ. output%jrng(k)) THEN
	  ! Record particle quantities:
          output%tp(k) = plasma%tp
          !$OMP PARALLEL DO
          DO i = 1,params%NC
             output%zp(i,k)    = plasma%zp(i)
             output%kep(i,k)   = plasma%kep(i)
             output%xip(i,k)   = plasma%xip(i)
             output%a(i,k)     = plasma%a(i)
             output%m(i,k)     = plasma%m(i)
             output%np(i,k)    = plasma%np(i)
             output%Up(i,k)    = plasma%Up(i)
             output%Tparp(i,k) = plasma%Tparp(i)
             output%Tperp(i,k) = plasma%Tperp(i)
             output%Bp(i,k)    = plasma%Bp(i)
             output%dBp(i,k)   = plasma%dBp(i)
             output%ddBp(i,k)  = plasma%ddBp(i)
          END DO
         !$OMP END PARALLEL DO

	 ! Record mesh quantities:
	 output%n(:,k)    = mesh%n
	 output%nU(:,k)   = mesh%nU
	 output%unU(:,k)  = mesh%unU
	 output%nUE(:,k)  = mesh%nUE
	 output%Ppar(:,k) = mesh%Ppar
	 output%Pper(:,k) = mesh%Pper
	 output%Tpar(:,k) = mesh%Tpar
	 output%Tper(:,k) = mesh%Tper
	 output%U(:,k)    = mesh%U
	 output%ddB(:,k)  = mesh%ddB
	 output%E(:,k)    = mesh%E
         ! Increment counter:
         k = k + 1
      END IF
    END IF

    ! Update time array:
    ! =========================================================================
    plasma%tp = plasma%tp + dt

    ! Estimate computational time:
    ! =====================================================================
    IF (j .EQ. 150) THEN
       oend_estimate = OMP_GET_WTIME()
       WRITE(*,*) 'Estimated compute time: ', params%NS*(oend_estimate-ostart)/j,' [s]'
    END IF

END DO NS_loop

! Record end time:
! =========================================================================
params%tSimTime = plasma%tp
output%zm = mesh%zm
oend = OMP_GET_WTIME()

params%tComputeTime = oend-ostart
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) 'Reached End of Program, Computational time [s] = ', params%tComputeTime
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Save data:
! ==============================================================================
IF (params%iSave) THEN
    ! Create new directory to save output data:
    ! --------------------------------------------------------------------------
    dir1 = TRIM(params%repoDir)//'/OutputFiles'
    command = 'mkdir '//dir1
    CALL SYSTEM(command,STATUS)
    WRITE(*,*) 'Status: ', STATUS

    ! Create subdirectory with input file name:
    ! -------------------------------------------------------------------------
    n_mpwd = LEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir1 = TRIM(params%repoDir)//'/OutputFiles/'//TRIM(dir0)
    command = 'mkdir '//dir1
    CALL SYSTEM(command,STATUS)

    ! Create subdirectory based on "fileDescriptor" inside the input file:
    ! -------------------------------------------------------------------------
    dir1 = TRIM(dir1)//'/'//TRIM(params%fileDescriptor)
    command = 'mkdir '// dir1
    CALL SYSTEM(command,STATUS)
    CALL GETCWD(mpwd)

    ! Save output data:
    ! -------------------------------------------------------------------------
    CALL SaveData(output,dir1)

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
    n_mpwd = LEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir0 = TRIM(params%repoDir)//'/OutputFiles/'//trim(dir0)//'/'//trim(params%fileDescriptor)
    fileName = TRIM(dir0)//'/commitHash.txt'
    command = 'git log --oneline -1 > '//TRIM(fileName)
    CALL SYSTEM(command)

END IF

END PROGRAM
