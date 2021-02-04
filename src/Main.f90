PROGRAM LinearFokkerPlanckSolver_1D2V
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
TYPE(splTYP) :: spline_Bz
TYPE(splTYP) :: spline_ddBz
TYPE(splTYP) :: spline_Phi
TYPE(spltestTYP) :: spline_Test
! DO loop indices:
INTEGER(i4) :: i,j,k,tt
! Pseudo random number seed:
INTEGER(i4) :: seed_size
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: seed
! Size of time interval:
INTEGER(i4) :: jsize
! Indices of time steps to save:
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: jrng
! Thread ID:
INTEGER(i4) :: id
! CPU time at start and end of computation:
REAL(r8) :: tstart, tend
! OPENMP computational time:
DOUBLE PRECISION :: ostart, oend, oend_estimate
! Cyclotron resonance number change:
REAL(r8) :: df
! Main simulation variables:
! simulation time:
REAl(r8) :: tp
! Particle position (zp), kinetic energy (kep), pitch angle (xip):
REAL(r8), DIMENSION(:,:)  , ALLOCATABLE :: xip, zp, kep
! subset of zp, kep and xip to save:
REAL(r8), DIMENSION(:,:), ALLOCATABLE :: zp_hist, kep_hist, xip_hist
! Subset of tp:
REAl(r8), DIMENSION(:)  , ALLOCATABLE :: t_hist
! Simulation diagnotics:
! Count the number of particles incident on (1) dump, (2) target, (3) cyclotron resonance, (4) slowing down
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: pcount1, pcount2, pcount3, pcount4
! Record the total energy of particle incident on (1) dump, (2) target, (3) cyclotron resonance, (4) slowing down
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: ecount1, ecount2, ecount3, ecount4
! Cyclotron resonance number:
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: fcurr, fnew
! Local simulation diagnostic counters:
REAL(r8) :: ecnt1, ecnt2, ecnt3, ecnt4
REAL(r8) :: pcnt1, pcnt2, pcnt3, pcnt4
! To store system commands and fileNames:
CHARACTER*300 :: command, mpwd
INTEGER(i4) :: n_mpwd, STATUS
CHARACTER*300 :: inputFileDir, inputFile, fileName, xpSelector, rootDir, dir0, dir1
! Declare user-defined functions:
REAL(r8) :: Interp1, diff1

! Create input namelist from the user-defined structures:
! ==============================================================================
namelist/in_nml/in

! Get root directory:
! ==============================================================================
CALL GET_ENVIRONMENT_VARIABLE('REPO_DIR',rootDir)

! Get input file name and directory:
! =============================================================================
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE',inputFile)
CALL GET_ENVIRONMENT_VARIABLE('INPUT_FILE_DIR',inputFileDir)

! Read input data into "in" structure:
! ==============================================================================
OPEN(unit=4,file=inputFileDir,status='old',form='formatted')
read(4,in_nml)
CLOSE(unit=4)

! Populate the in%rootDir field:
! =============================================================================
in%rootDir = trim(adjustl(rootDir))

if (in%species_a .eq. 1) then
    in%q_t = -e_c
    in%m_t = m_e
else
    in%q_t = +in%Zion*e_c
    in%m_t = in%Aion*m_p
end if

! Print to the terminal:
! ==============================================================================
WRITE(*,*) ''
WRITE(*,*) '*********************************************************************'
WRITE(*,*) 'Input file:         ', TRIM(inputFile)
WRITE(*,*) 'fileDescriptor:     ', TRIM(in%fileDescriptor)
WRITE(*,*) 'Number of particles:', in%Nparts
WRITE(*,*) 'Number of steps:    ', in%Nsteps
WRITE(*,*) 'Particle BC:        ', in%particleBC
WRITE(*,*) 'dt [ns]:            ', in%dt*1E+9
WRITE(*,*) 'iPush:              ', in%iPush
WRITE(*,*) 'iDrag:              ', in%iDrag
WRITE(*,*) 'iColl:              ', in%iColl
WRITE(*,*) 'iHeat:              ', in%iHeat
WRITE(*,*) 'iSave:              ', in%iSave
WRITE(*,*) 'elevel:             ', in%elevel
WRITE(*,*) 'zTarget [m]:        ', in%zmax
WRITE(*,*) 'zDump [m]:          ', in%zmin
WRITE(*,*) 'zp_init [m]:        ', in%zp_init
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
CALL InitSpline(spline_Bz  ,in%nz,0._8,0._8,1,0._8)
CALL InitSpline(spline_ddBz,in%nz,0._8,0._8,1,0._8)
CALL InitSpline(spline_Phi ,in%nz,0._8,0._8,1,0._8)

! Allocate memory for spline_test variable:
! ==============================================================================
! Spline_test is a variable that holds up to 6 different test profiles "y" and
! one "x" coordinate.
CALL InitSplineTest(spline_Test,in%nz)

! Allocate memory to main simulation variables:
! ==============================================================================
ALLOCATE(zp(in%Nparts), kep(in%Nparts), xip(in%Nparts))
ALLOCATE(pcount1(in%Nsteps),pcount2(in%Nsteps),pcount3(in%Nsteps),pcount4(in%Nsteps))
ALLOCATE(ecount1(in%Nsteps),ecount2(in%Nsteps),ecount3(in%Nsteps),ecount4(in%Nsteps))

! Allocate memory to output variables:
! ==============================================================================
! Determine size of temporal snapshots to record:
jsize = (in%jend-in%jstart+1)/in%jincr
ALLOCATE(jrng(jsize))
ALLOCATE(zp_hist(in%Nparts,jsize),kep_hist(in%Nparts,jsize),xip_hist(in%Nparts,jsize),t_hist(jsize))

! Create array with the indices of the time steps to save:
jrng = (/ (j, j=in%jstart, in%jend, in%jincr) /)

! Allocate memory to local variables:
! ==============================================================================
! Cyclotron resonance number:
ALLOCATE(fcurr(in%Nparts),fnew(in%Nparts))

! Create splines of electromagnetic fields:
! ===========================================================================
! Magnetic field data:
fileName = trim(adjustl(in%BFieldFile))
fileName = trim(adjustl(in%BFieldFileDir))//fileName
fileName = trim(adjustl(in%rootDir))//fileName

! Populate profile data based on external file:
CALL ReadSpline(spline_Bz,fileName)

! Based on profile data. compute spline data:
CALL ComputeSpline(spline_Bz)

! Second derivative of the magnetic field:
spline_ddBz%x = spline_Bz%x
spline_ddBz%y = spline_Bz%yp
CALL ComputeSpline(spline_ddBz)

! Electric potential:
spline_Phi%x = spline_Bz%x
spline_Phi%y = 0
if (in%iPotential) then
  CALL PotentialProfile(spline_Phi,in)
end if
CALL ComputeSpline(spline_Phi)

! Test the splines:
! ==========================================================================
if (.true.) then
    do i=1,in%nz
        !spline_Test%x(i)  = in%zmin + i*0.03
        spline_Test%x(i)  = 0.0 + i*0.02
        spline_Test%y1(i) = Interp1(spline_Test%x(i),spline_Bz)
        spline_Test%y2(i) = diff1(spline_Test%x(i),spline_Bz)
        spline_Test%y3(i) = Interp1(spline_Test%x(i),spline_ddBz)
        spline_Test%y4(i) = Interp1(spline_Test%x(i),spline_Phi)
        spline_Test%y5(i) = BESSEL_JN(0,spline_Test%x(i))
        spline_Test%y6(i) = BESSEL_JN(1,spline_Test%x(i))
    end do
    fileName = "B_spline.dat"
    OPEN(unit=8,file=fileName,form="formatted",status="unknown")
    do j = 1,in%nz
        WRITE(8,*) spline_Test%x(j), spline_Test%y1(j), spline_Test%y2(j),&
         spline_Test%y3(j), spline_Test%y4(j), spline_Test%y5(j), spline_Test%y6(j)
    end do
    CLOSE(unit=8)
end if

! Initialize pseudo random number generator:
! ===========================================================================
CALL random_seed(size=seed_size)
ALLOCATE(seed(seed_size))
CALL random_seed(get=seed)
seed = 314159565
CALL random_seed(put=seed)

! Inititalize zp, kep, xip
! ==============================================================================
kep = 0.; xip = 0.; zp = 0.;
WRITE(*,*) "Initializing PDF..."
!$OMP PARALLEL DO
DO i = 1,in%Nparts
  CALL loadParticles(zp(i),kep(i),xip(i),in)
END DO
!$OMP END PARALLEL DO
WRITE(*,*) "Initialization complete"

! Test initial distribution:
! ==========================================================================
fileName = "LoadParticles.dat"
OPEN(unit=8,file=fileName,form="formatted",status="unknown")
do i = 1,in%Nparts
    WRITE(8,*) zp(i), kep(i), xip(i)
end do
CLOSE(unit=8)

! Initialize simulation diagnostics:
! ==============================================================================
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
TimeStepping: do j = 1,in%Nsteps

!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP SHARED(j,zp,kep,xip,fcurr,fnew,in,ecount1,ecount2,ecount3,ecount4,pcount1,pcount2,pcount3,pcount4) , &
!$OMP& SHARED(spline_ddBz,spline_Bz,spline_Phi)

          ! ! PRIVATE(pcnt1,pcnt2,pcnt3,pcnt4,ecnt1,ecnt2,ecnt3,ecnt4,i,df)

          ! Initialize private particle counters:
          pcnt1 = 0; pcnt2 = 0; pcnt3 = 0; pcnt4 = 0
          ! Initialize private energy counters:
          ecnt1 = 0; ecnt2 = 0; ecnt3 = 0; ecnt4 = 0
          ! Initialize resonance number difference:
          df = 0

          ! Loop over particles:
          ! ==============================================================================
          !$OMP DO
              AllParticles: do i = 1,in%Nparts

            		if (j .EQ. 1 .AND. i .EQ. 1) then
	                WRITE(*,*) "spline_ddBz", spline_ddBz%n
	                WRITE(*,*) "spline_Phi", spline_Phi%n
	                WRITE(*,*) "spline_Bz", spline_Bz%n
			            in%threads_given = OMP_GET_NUM_THREADS()
            		  WRITE(*,*) ''
            		  WRITE(*,*) '*********************************************************************'
            		  WRITE(*,*) "Number of threads given: ", in%threads_given
            		  WRITE(*,*) "Padding: ", in%padding, " Real(r8)"
            		  WRITE(*,*) '*********************************************************************'
            		  WRITE(*,*) ''
            		end if

                ! Calculate Cyclotron resonance number:
                ! ------------------------------------------------------------------------
                if (in%iHeat) CALL CyclotronResonanceNumber(zp(i),kep(i),xip(i),fcurr(i),in,spline_Bz)
                ! fcurr and fnew could be declared private

                ! Push particles adiabatically:
                ! ------------------------------------------------------------------------
		            if (in%iPush) CALL MoveParticle(zp(i),kep(i),xip(i),in,spline_Bz,spline_Phi)

                ! Re-inject particles:
                ! ------------------------------------------------------------------------
                if (zp(i) .GE. in%zmax) CALL ReinjectParticles(zp(i),kep(i),xip(i),in,ecnt2,pcnt2)
                if (zp(i) .LE. in%zmin) CALL ReinjectParticles(zp(i),kep(i),xip(i),in,ecnt1,pcnt1)

                ! Apply Coulomb collision operator:
                ! ------------------------------------------------------------------------
                if (in%iColl) then
                    ! "in" needs to be private to avoid race condition. This can be
                    ! fixed by looping over species inside the subroutine "collisionOperator"
                    in%species_b = 1
                    CALL collisionOperator(zp(i),kep(i),xip(i),ecnt4,pcnt4,in)
                    in%species_b = 2
                    CALL collisionOperator(zp(i),kep(i),xip(i),ecnt4,pcnt4,in)
                end if

                ! Apply RF heating operator:
                ! ------------------------------------------------------------------------
                if (in%iHeat) then
                  CALL CyclotronResonanceNumber(zp(i),kep(i),xip(i),fnew(i),in,spline_Bz)
                  df = dsign(1.d0,fcurr(i)*fnew(i))
                  if (df .LT. 0 .AND. zp(i) .GT. in%zRes1 .AND. zp(i) .LT. in%zRes2)  then
                    CALL RFHeatingOperator(zp(i),kep(i),xip(i),ecnt3,pcnt3,in,spline_Bz,spline_ddBz,spline_Phi)
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

    ! Update time array:
    ! =========================================================================
    tp = tp + in%dt

    ! Select data to save:
    ! =====================================================================
    ! Check if data is to be saved
    if (in%iSave) then
        do k = 1,jsize
            if (j .EQ. jrng(k)) then
                t_hist(k) = tp
                    !$OMP PARALLEL DO PRIVATE(i)
                    do i = 1,in%Nparts
                            ! Record "ith" particle position at "kth" time
                            zp_hist(i,k) = zp(i)
                            ! Record "ith" particle KE at "kth" time
                            kep_hist(i,k) = kep(i)
                            ! Record "ith" particle pitch angle at "kth" time
                            xip_hist(i,k) = xip(i)
                    end do
                   !$OMP END PARALLEL DO
            endif
        end do
    end if

    ! Estimate computational time:
    ! =====================================================================
    id = OMP_GET_THREAD_NUM()
    if (id .EQ. 0) then
      if (j .EQ. 50) then
	       oend_estimate = OMP_GET_WTIME()
         WRITE(*,*) 'Estimated compute time: ', in%Nsteps*(oend_estimate-ostart)/j,' [s]'
      end if
   end if

end do TimeStepping

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
    dir1 = trim(in%rootDir)//'/OutputFiles'
    command = 'mkdir '//dir1
    CALL system(command,STATUS)
    WRITE(*,*) 'Status: ', STATUS

    ! Create subdirectory with input file name:
    ! -------------------------------------------------------------------------
    n_mpwd = lEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir1 = trim(in%rootDir)//'/OutputFiles/'//trim(dir0)
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
    dir0 = trim(in%rootDir)//'/InputFiles/'//trim(inputFile)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))
    CALL system(command)

    ! Copy magnetic field data:
    ! --------------------------------------------------------------------------
    dir0 = trim(in%rootDir)//'/BfieldData'//trim(in%BFieldFile)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))//'/Bfield.txt'
    CALL system(command)

    ! Create text file with commit Hash:
    ! --------------------------------------------------------------------------
    n_mpwd = lEN_TRIM(inputFile)-3
    dir0 = inputFile
    dir0 = dir0(1:n_mpwd)
    dir0 = trim(in%rootDir)//'/OutputFiles/'//trim(dir0)//'/'//trim(in%fileDescriptor)
    fileName = trim(dir0)//'/commitHash.txt'
    command = 'git log --oneline -1 > '//trim(fileName)
    CALL system(command)

end if

End PROGRAM
