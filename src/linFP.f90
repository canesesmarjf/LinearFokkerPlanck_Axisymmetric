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
TYPE(inTYP)  :: in
TYPE(splTYP) :: spline_B
TYPE(splTYP) :: spline_dB
TYPE(splTYP) :: spline_ddB
TYPE(splTYP) :: spline_V
TYPE(splTYP) :: spline_dV
TYPE(plasmaTYP) :: plasma
! DO loop indices:
INTEGER(i4) :: i,j,k
! Size of time interval:
INTEGER(i4) :: jsize
! Indices of time steps to save:
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: jrng
! Thread ID:
INTEGER(i4) :: id
! OPENMP computational time:
DOUBLE PRECISION :: ostart, oend, oend_estimate
! Cyclotron resonance numbers:
REAL(r8) :: dresNum, resNum0, resNum1
! Main simulation variables:
! simulation time:
REAl(r8) :: tp
! subset of zp, kep and xip to save:
REAL(r8), DIMENSION(:,:), ALLOCATABLE :: zp_hist, kep_hist, xip_hist
! Subset of tp:
REAl(r8), DIMENSION(:)  , ALLOCATABLE :: t_hist
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

! Create input namelist from the user-defined structures:
! ==============================================================================
namelist/in_nml/in

! Get root directory:
! ==============================================================================
CALL GET_ENVIRONMENT_VARIABLE('REPO_DIR',repoDir)

! Get input file name and directory:
! =============================================================================
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE',inputFile)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE_DIR',inputFileDir)

! Read input data into "in" structure:
! ==============================================================================
OPEN(unit=4,file=inputFileDir,status='old',form='formatted')
read(4,in_nml)
CLOSE(unit=4)

! Select the test species:
! ==============================================================================
if (in%species_a .eq. 1) then
    in%qa = -e_c
    in%Ma = m_e
else
    in%qa = +in%Zion*e_c
    in%Ma = in%Aion*m_p
end if

! Print to the terminal:
! ==============================================================================
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) 'Input file:         ', TRIM(inputFile)
WRITE(*,*) 'fileDescriptor:     ', TRIM(in%fileDescriptor)
WRITE(*,*) 'Number of particles:', in%Nparts
WRITE(*,*) 'Number of steps:    ', in%Nsteps
WRITE(*,*) 'Particle BC:        ', in%BC_Type
WRITE(*,*) 'dt [ns]:            ', in%dt*1E+9
WRITE(*,*) 'iPush:              ', in%iPush
WRITE(*,*) 'iDrag:              ', in%iDrag
WRITE(*,*) 'iColl:              ', in%iColl
WRITE(*,*) 'iHeat:              ', in%iHeat
WRITE(*,*) 'iSave:              ', in%iSave
WRITE(*,*) 'elevel:             ', in%elevel
WRITE(*,*) 'zTarget [m]:        ', in%zmax
WRITE(*,*) 'zDump [m]:          ', in%zmin
WRITE(*,*) 'BC_zp_mean [m]:     ', in%BC_zp_mean
WRITE(*,*) 'B field file:       ', TRIM(in%BFieldFile)
WRITE(*,*) 'Ew:                 ', in%Ew
WRITE(*,*) 'Te0:                ', in%Te0
WRITE(*,*) 'ne0:                ', in%ne0
if (in%CollOperType .EQ. 1) WRITE(*,*) 'Boozer-Only collision operator'
if (in%CollOperType .EQ. 2) WRITE(*,*) 'Boozer-Kim collision operator'
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Allocate memory to splines:
! ===========================================================================
! InitSpline takes in a variable of splType and performs the following:
! - Populates it scalar fields
! - Allocates memory for the profile data "y", indepedent coordinate "x", etc
CALL InitSpline(spline_B  ,in%nz,0._8,0._8)
CALL InitSpline(spline_dB ,in%nz,0._8,0._8)
CALL InitSpline(spline_ddB,in%nz,0._8,0._8)
CALL InitSpline(spline_V  ,in%nz,0._8,0._8)
CALL InitSpline(spline_dV ,in%nz,0._8,0._8)

! Allocate memory to main simulation variables:
! ==============================================================================
ALLOCATE(plasma%zp(in%Nparts),plasma%kep(in%Nparts),plasma%xip(in%Nparts))
ALLOCATE(pcount1(in%Nsteps)  ,pcount2(in%Nsteps)   ,pcount3(in%Nsteps)  ,pcount4(in%Nsteps))
ALLOCATE(ecount1(in%Nsteps)  ,ecount2(in%Nsteps)   ,ecount3(in%Nsteps)  ,ecount4(in%Nsteps))

! Allocate memory to output variables:
! ==============================================================================
! Determine size of temporal snapshots to record:
jsize = (in%jend-in%jstart+1)/in%jincr
ALLOCATE(jrng(jsize))
ALLOCATE(zp_hist(in%Nparts,jsize),kep_hist(in%Nparts,jsize),xip_hist(in%Nparts,jsize),t_hist(jsize))

! Create array with the indices of the time steps to save:
jrng = (/ (j, j=in%jstart, in%jend, in%jincr) /)

! Create splines of electromagnetic fields:
! ===========================================================================
! Magnetic field data:
fileName = trim(adjustl(in%BFieldFile))
fileName = trim(adjustl(in%BFieldFileDir))//fileName
fileName = trim(adjustl(in%repoDir))//fileName

! Populate spline profile data using external file:
! B:
CALL ReadSpline(spline_B,fileName)
! dB and ddB:
CALL diffSpline(spline_B ,spline_dB )
CALL diffSpline(spline_dB,spline_ddB)

! Electric potential V:
spline_V%x = spline_B%x
spline_V%y = 0
if (in%iPotential) then
  CALL PotentialProfile(spline_V,in)
end if
CALL ComputeSpline(spline_V)

! dV:
CALL diffSpline(spline_V,spline_dV)

! Complete setting up the spline data :
CALL ComputeSpline(spline_B  )
CALL ComputeSpline(spline_dB )
CALL ComputeSpline(spline_ddB)
CALL ComputeSpline(spline_V  )
CALL ComputeSpline(spline_dV )

! Inititalize zp, kep, xip
! ==============================================================================
plasma%kep = 0.; plasma%xip = 0.; plasma%zp = 0.;
WRITE(*,*) "Initializing PDF..."
!$OMP PARALLEL
if (OMP_GET_THREAD_NUM() .EQ. 0) then
    in%threads_given = OMP_GET_NUM_THREADS()
    WRITE(*,*) ''
    WRITE(*,*) '*********************************************************************'
    WRITE(*,*) "Number of threads given: ", in%threads_given
    WRITE(*,*) '*********************************************************************'
    WRITE(*,*) ''
end if
!$OMP DO
DO i = 1,in%Nparts
  CALL loadParticles(plasma%zp(i),plasma%kep(i),plasma%xip(i),in)
END DO
!$OMP END DO
!$OMP END PARALLEL
WRITE(*,*) "Initialization complete"

! Test initial distribution:
! ==========================================================================
fileName = "LoadParticles.dat"
OPEN(unit=8,file=fileName,form="formatted",status="unknown")
do i = 1,in%Nparts
    WRITE(8,*) plasma%zp(i), plasma%kep(i), plasma%xip(i)
end do
CLOSE(unit=8)

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
AllTime: do j = 1,in%Nsteps

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
    AllParticles: do i = 1,in%Nparts

        ! Calculate Cyclotron resonance number:
        ! ------------------------------------------------------------------------
        if (in%iHeat) CALL CyclotronResonanceNumber(i,plasma,resNum0,in,spline_B)

        ! Push particles adiabatically:
        ! ------------------------------------------------------------------------
        if (in%iPush) CALL MoveParticle(i,plasma,in,spline_B,spline_dB,spline_dV)

        ! Re-inject particles:
        ! ------------------------------------------------------------------------
        if (plasma%zp(i) .GE. in%zmax) CALL ReinjectParticles(i,plasma,in,ecnt2,pcnt2)
        if (plasma%zp(i) .LE. in%zmin) CALL ReinjectParticles(i,plasma,in,ecnt1,pcnt1)

        ! Apply Coulomb collision operator:
        ! ------------------------------------------------------------------------
        if (in%iColl) CALL collisionOperator(i,plasma,ecnt4,pcnt4,in)

        ! Apply RF heating operator:
        ! ------------------------------------------------------------------------
        if (in%iHeat) then
           CALL CyclotronResonanceNumber(i,plasma,resNum1,in,spline_B)
           dresNum = dsign(1.d0,resNum0*resNum1)
           if (dresNum .LT. 0 .AND. plasma%zp(i) .GT. in%zRes1 .AND. plasma%zp(i) .LT. in%zRes2)  then
              !WRITE(*,*) 'Resonance at zp: ', zp(i)
        CALL RFHeatingOperator(plasma%zp(i),plasma%kep(i),plasma%xip(i),ecnt3,pcnt3,in,spline_B,spline_dB,spline_ddB,spline_dV)
           end if
        end if

    end do AllParticles
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
    if (in%iSave) then
       if (j .EQ. jrng(k)) then
          t_hist(k) = tp
          !$OMP PARALLEL DO
          do i = 1,in%Nparts
             ! Record "ith" particle at "kth" time
             zp_hist(i,k)  = plasma%zp(i)
             kep_hist(i,k) = plasma%kep(i)
             xip_hist(i,k) = plasma%xip(i)
          end do
         !$OMP END PARALLEL DO
         k = k + 1
      end if
    end if

    ! Update time array:
    ! =========================================================================
    tp = tp + in%dt

    ! Estimate computational time:
    ! =====================================================================
    if (j .EQ. 150) then
       oend_estimate = OMP_GET_WTIME()
       WRITE(*,*) 'Estimated compute time: ', in%Nsteps*(oend_estimate-ostart)/j,' [s]'
    end if

end do AllTime

! Record end time:
! =========================================================================
in%tSimTime = tp
oend = OMP_GET_WTIME()

in%tComputeTime = oend-ostart
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) 'Reached End of Program, Computational time [s] = ', in%tComputeTime
WRITE(*,*) '*********************************************************************'
WRITE(*,*) ''

! Save data:
! ==============================================================================
if (in%iSave) then
    ! Create new directory to save output data:
    ! --------------------------------------------------------------------------
    dir1 = trim(in%repoDir)//'/OutputFiles'
    command = 'mkdir '//dir1
    CALL system(command,STATUS)
    WRITE(*,*) 'Status: ', STATUS

    ! Create subdirectory with input file name:
    ! -------------------------------------------------------------------------
    n_mpwd = lEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir1 = trim(in%repoDir)//'/OutputFiles/'//trim(dir0)
    command = 'mkdir '//dir1
    CALL system(command,STATUS)

    ! Create subdirectory based on "fileDescriptor" inside the input file:
    ! -------------------------------------------------------------------------
    dir1 = trim(dir1)//'/'//trim(in%fileDescriptor)
    command = 'mkdir '// dir1
    CALL system(command,STATUS)
    CALL getcwd(mpwd)

    ! Saving zp_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'zp.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) zp_hist
    CLOSE(unit=8)

    ! Saving kep_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'kep.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) kep_hist
    CLOSE(unit=8)

    ! Saving xip_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'xip.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) xip_hist
    CLOSE(unit=8)

    ! Saving t_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'tp.out')
    OPEN(unit=8,file=fileName,form="unformatted",status="unknown")
    WRITE(8) t_hist
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
    WRITE(8,NML = in_nml)
    CLOSE(unit=8)

    ! Copy input file to output directory:
    ! --------------------------------------------------------------------------
    dir0 = trim(in%repoDir)//'/InputFiles/'//trim(inputFile)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))
    CALL system(command)

    ! Copy magnetic field data:
    ! --------------------------------------------------------------------------
    dir0 = trim(in%repoDir)//'/BfieldData'//trim(in%BFieldFile)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))//'/Bfield.txt'
    CALL system(command)

    ! Create text file with commit Hash:
    ! --------------------------------------------------------------------------
    n_mpwd = lEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir0 = trim(in%repoDir)//'/OutputFiles/'//trim(dir0)//'/'//trim(in%fileDescriptor)
    fileName = trim(dir0)//'/commitHash.txt'
    command = 'git log --oneline -1 > '//trim(fileName)
    CALL system(command)

END if

END PROGRAM
