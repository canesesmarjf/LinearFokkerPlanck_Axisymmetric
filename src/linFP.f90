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
TYPE(outputTYP)      :: output
! Accumulators:
REAL(r8) :: n1,n2,n3,n4
REAL(r8) :: e1,e2,e3,e4
REAL(r8) :: e3_hat
! Super particle to real particle conversion factor:
REAL(r8) :: SP2RP
! DO loop indices:
INTEGER(i4) :: i,j,k
! Thread ID:
INTEGER(i4) :: id
! OPENMP computational time:
DOUBLE PRECISION :: ostart, oend, oend_estimate
! Cyclotron resonance numbers:
REAL(r8) :: dresNum, resNum0, resNum1
! simulation time:
REAl(r8) :: tp
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

! Allocate memory:
! ===========================================================================
CALL AllocateFieldSpline(fieldspline,params)
CALL AllocatePlasma(plasma,params)
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

! Electric potential V, dV:
fieldspline%V%x = fieldspline%B%x
fieldspline%V%y = 0
IF (params%iPotential) THEN
  CALL PotentialProfile(fieldspline,params)
END IF
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
! Time array:
tp = 0

! Loop over time:
! ==============================================================================
NS_loop: DO j = 1,params%NS
    ! Reset all accumulators:
    n1 = 0. ; n2 = 0. ; n3 = 0. ; n4 = 0.;
    e1 = 0. ; e2 = 0. ; e3 = 0. ; e4 = 0.;

    ! Reset resonance number:
    dresNum = 0.; resNum0 = 0.; resNum1 = 0.

    ! Reset reference RF power:
    e3_hat = 0.

    ! Calculate scale factor: super-particle to real-particle
    SP2RP = plasma%NR(j)/plasma%NSP(j)

    !$OMP PARALLEL FIRSTPRIVATE(dresNum,resNum0,resNum1)
    !$OMP DO SCHEDULE(STATIC)
    NC_loop1: DO i = 1,params%NC

        ! Initialize plasma: flags and Energy increments
        ! ------------------------------------------------------------------------
        CALL ResetFlags(i,plasma)

        ! Calculate Cyclotron resonance number:
        ! ------------------------------------------------------------------------
        IF (params%iHeat) CALL CyclotronResonanceNumber(i,plasma,resNum0,fieldspline,params)

        ! Push particles adiabatically:
        ! ------------------------------------------------------------------------
        IF (params%iPush) CALL MoveParticle(i,plasma,fieldspline,params)

        ! Check boundaries:
        ! ------------------------------------------------------------------------
        IF (plasma%zp(i) .GE. params%zmax) THEN
           plasma%f2(i) = 1; plasma%E2(i) = plasma%kep(i)
        ELSE IF (plasma%zp(i) .LE. params%zmin) THEN
           plasma%f1(i) = 1; plasma%E1(i) = plasma%kep(i)
        END IF

        ! Apply Coulomb collision operator:
        ! ------------------------------------------------------------------------
        IF (params%iColl) CALL collisionOperator(i,plasma,params)

        ! Check resonance condition:
        ! ------------------------------------------------------------------------
        IF (params%iHeat) THEN
           CALL CyclotronResonanceNumber(i,plasma,resNum1,fieldspline,params)
           dresNum = dsign(1.d0,resNum0*resNum1)
           IF (dresNum .LT. 0 .AND. plasma%zp(i) .GT. params%zRes1 .AND. plasma%zp(i) .LT. params%zRes2)  THEN
              plasma%f3(i) = 1
              CALL RFoperatorTerms(i,plasma,fieldspline,params)
           END IF
        END IF

    END DO NC_loop1
    !$OMP END DO

    !$OMP DO REDUCTION(+:n1,n2,n3,n4,e1,e2,e3_hat,e4)
    ! Calculate particle and energy rates:
    NC_loop2: DO i = 1,params%NC
       n1     = n1     +     (SP2RP/params%dt)*plasma%a(i)*plasma%f1(i)
       n2     = n2     +     (SP2RP/params%dt)*plasma%a(i)*plasma%f2(i)
       n3     = n3     +     (SP2RP/params%dt)*plasma%a(i)*plasma%f3(i)
       n4     = n4     +     (SP2RP/params%dt)*plasma%a(i)*plasma%f4(i)
       e1     = e1     + e_c*(SP2RP/params%dt)*plasma%a(i)*plasma%f1(i)*plasma%E1(i)
       e2     = e2     + e_c*(SP2RP/params%dt)*plasma%a(i)*plasma%f2(i)*plasma%E2(i)
       e3_hat = e3_hat + e_c*(SP2RP/params%dt)*plasma%a(i)*plasma%f3(i)*plasma%E3_hat(i)
       e4     = e4     + e_c*(SP2RP/params%dt)*plasma%a(i)*plasma%f4(i)*plasma%E4(i)
    END DO NC_loop2
    !$OMP END DO
  
    !$OMP DO
    ! Apply particle re-injection and RF operator:
    NC_loop3: DO i = 1,params%NC
       IF (plasma%f1(i) .EQ. 1 .OR. plasma%f2(i) .EQ. 1) THEN
           CALL ReinjectParticles(i,plasma,params)
       END IF
       IF (plasma%f3(i) .EQ. 1) THEN
           CALL RFOperator(i,plasma,fieldspline,params,e3_hat)
       END IF
    END DO NC_loop3
    !$OMP END DO

    !$OMP DO REDUCTION(+:e3)
    ! Calculate RF power absorption:
    NC_loop4: DO i = 1,params%NC
       e3 = e3 + e_c*(SP2RP/params%dt)*plasma%a(i)*plasma%f3(i)*plasma%E3(i)
    END DO NC_loop4
    !$OMP END DO

    !$OMP END PARALLEL

    ! Assign particle and energy rates in physical units [P/s] and [J/s]
    ! ==============================================================================
    plasma%Ndot1(j) = n1
    plasma%Ndot2(j) = n2
    plasma%Ndot3(j) = n3
    plasma%Ndot4(j) = n4
    plasma%Edot1(j) = e1
    plasma%Edot2(j) = e2
    plasma%Edot3(j) = e3
    plasma%Edot4(j) = e4

    ! Select data to save:
    ! =====================================================================
    ! Check if data is to be saved
    IF (params%iSave) THEN
       IF (j .EQ. output%jrng(k)) THEN
          output%tp(k) = tp
          !$OMP PARALLEL DO
          DO i = 1,params%NC
             ! Record "ith" particle at "kth" time
             output%zp(i,k)  = plasma%zp(i)
             output%kep(i,k) = plasma%kep(i)
             output%xip(i,k) = plasma%xip(i)
          END DO
         !$OMP END PARALLEL DO
         k = k + 1
      END IF
    END IF

    ! Update time array:
    ! =========================================================================
    tp = tp + params%dt

    ! Estimate computational time:
    ! =====================================================================
    IF (j .EQ. 150) THEN
       oend_estimate = OMP_GET_WTIME()
       WRITE(*,*) 'Estimated compute time: ', params%NS*(oend_estimate-ostart)/j,' [s]'
    END IF

END DO NS_loop

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

    ! Saving zp_hist to file:
    ! --------------------------------------------------------------------------
    fileName = TRIM(TRIM(dir1)//'/'//'zp.out')
    OPEN(UNIT=8,FILE=fileName,FORM="unformatted",STATUS="unknown")
    WRITE(8) output%zp
    CLOSE(UNIT=8)

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
    WRITE(8) plasma%Ndot1
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount2.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) plasma%Ndot2
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount3.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) plasma%Ndot3
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount4.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) plasma%Ndot4
    CLOSE(unit=8)

    ! Saving ecount to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'ecount1.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) plasma%Edot1
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount2.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) plasma%Edot2
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount3.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) plasma%Edot3
    CLOSE(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount4.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) plasma%Edot4
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
    n_mpwd = LEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir0 = TRIM(params%repoDir)//'/OutputFiles/'//trim(dir0)//'/'//trim(params%fileDescriptor)
    fileName = TRIM(dir0)//'/commitHash.txt'
    command = 'git log --oneline -1 > '//TRIM(fileName)
    CALL SYSTEM(command)

END IF

END PROGRAM
