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
INTEGER(i4) :: i,j,k
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
! openMP computational time:
DOUBLE PRECISION :: ostart, oend
! Cyclotron resonance number change:
REAL(r8) :: df
! Main simulation variables:
! simulation time:
REAl(r8) :: tp
! Particle position (zp), kinetic energy (kep), pitch angle (xip):
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: xip, zp, kep
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
! Local simulation diagnostics:
REAL(r8) :: ecnt, ecnt1, ecnt2
REAL(r8) :: pcnt, pcnt1, pcnt2
! To store system commands and fileNames:
CHARACTER*300 :: command, mpwd
INTEGER(i4) :: n_mpwd, STATUS
CHARACTER*300 :: fileName, xpSelector, rootDir, dir0, dir1

! Create input namelist from the user-defined structures:
! ==============================================================================
namelist/in_nml/in
namelist/xp_nml/xpSelector

! Get root directory:
! ==============================================================================
call getcwd(mpwd)
! Get length of directory string:
n_mpwd = LEN(trim(adjustl(mpwd)))
! Remove parts of string:
n_mpwd = n_mpwd - 4
rootDir = mpwd(1:n_mpwd)

! Read input data into "in" structure:
! ==============================================================================
! Read contents of xpSelector:
fileName = "/InputFiles/xpSelector.in"
fileName = trim(adjustl(rootDir))//fileName
open(unit=4,file=fileName,status='old',form='formatted')
read(4,xp_nml)
close(unit=4)

! Read the file with name given by contents of xpSelector:
fileName = trim(adjustl(rootDir))//"/InputFiles/"//trim(adjustl(xpSelector))
open(unit=4,file=fileName,status='old',form='formatted')
read(4,in_nml)
close(unit=4)

! populate the in%rootDir field:
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
print *, 'Input file:         ', xpSelector
print *, 'fileDescriptor:     ', in%fileDescriptor
print *, 'Number of particles:', in%Nparts
print *, 'Number of steps:    ', in%Nsteps
print *, 'dt [ns]:            ', in%dt*1E+9
print *, 'iPush:              ', in%iPush
print *, 'iDrag:              ', in%iDrag
print *, 'iColl:              ', in%iColl
print *, 'iHeat:              ', in%iHeat
print *, 'iSave:              ', in%iSave
print *, 'elevel:             ', in%elevel
print *, 'zTarget [m]:        ', in%zmax
print *, 'zDump [m]:          ', in%zmin
print *, 'zp_init [m]:        ', in%zp_init
print *, 'B field file:       ', in%BFieldFile
print *, 'Ew:                 ', in%Ew
print *, 'Te0:                 ', in%Te0
print *, 'ne0:                ', in%ne0
if (in%CollOperType .EQ. 1) print *, 'Boozer-Only collision operator'
if (in%CollOperType .EQ. 2) print *, 'Boozer-Kim collision operator'

! Allocate memory to splines:
! ===========================================================================
call InitSpline(spline_Bz  ,in%nz,0._8,0._8,1,0._8)
call InitSpline(spline_ddBz,in%nz,0._8,0._8,1,0._8)
call InitSpline(spline_Phi ,in%nz,0._8,0._8,1,0._8)
call InitSplineTest(spline_Test,in%nz)

! Allocate memory to main simulation variables:
! ==============================================================================
ALLOCATE(zp(in%Nparts), kep(in%Nparts), xip(in%Nparts))
ALLOCATE(pcount1(in%Nsteps),pcount2(in%Nsteps),pcount3(in%Nsteps),pcount4(in%Nsteps))
ALLOCATE(ecount1(in%Nsteps),ecount2(in%Nsteps),ecount3(in%Nsteps),ecount4(in%Nsteps))

! Allocate memory to output variables:
! ==============================================================================
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
CALL ReadSpline(spline_Bz,fileName)
CALL ComputeSpline(spline_Bz)

! Second derivative of the magnetic field:
spline_ddBz%x = spline_Bz%x
spline_ddBz%y = spline_Bz%yp
CALL ComputeSpline(spline_ddBz)

! Predefined electric potential:
spline_Phi%x = spline_Bz%x
spline_Phi%y = 0
if (in%iPotential) then
  call PotentialProfile(spline_Phi,in)
end if
CALL ComputeSpline(spline_Phi)

! Test the splines:
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
    open(unit=8,file=fileName,form="formatted",status="unknown")
    do j = 1,in%nz
        write(8,*) spline_Test%x(j), spline_Test%y1(j), spline_Test%y2(j),&
         spline_Test%y3(j), spline_Test%y4(j), spline_Test%y5(j), spline_Test%y6(j)
    end do
    close(unit=8)
end if

! Record start time:
! ==============================================================================
ostart = OMP_GET_WTIME()

! Initialize pseudo random number generator:
! ===========================================================================
call random_seed(size=seed_size)
ALLOCATE(seed(seed_size))
call random_seed(get=seed)
seed = 314159565
call random_seed(put=seed)

! Inititalize zp, kep, xip
! ==============================================================================
kep = 0.; xip = 0.; zp = 0.;
call loadParticles(zp,kep,xip,in)

! Test initial distribution:
fileName = "LoadParticles.dat"
open(unit=8,file=fileName,form="formatted",status="unknown")
do i = 1,in%Nparts
    write(8,*) zp(i), kep(i), xip(i)
end do
close(unit=8)

! Initialize simulation diagnostics:
! ==============================================================================
! Time array:
tp = 0
! Leak current diagnostics:
pcount1 = 0; pcount2 = 0; pcount3 = 0; pcount4 = 0
! Energy leak diagnotics:
ecount1 = 0; ecount2 = 0; ecount3 = 0; ecount4 = 0

! OMP setup:
! ==============================================================================
! Set the number of threads:
call OMP_SET_NUM_THREADS(in%threads_request)
!$OMP PARALLEL PRIVATE(id)
    id = OMP_GET_THREAD_NUM()
    in%threads_given = OMP_GET_NUM_THREADS()
    if (id .EQ. 0) write(*,*) "number of threads given: ", in%threads_given
    if (id .EQ. 0) write(*,*) "chunk size: ", in%chunk
!$OMP END PARALLEL

! Start of simulation:
! ==============================================================================
! Begin time stepping:
TimeStepping: do j = 1,in%Nsteps

    ! Calculate Cyclotron resonance number:
    ! ==========================================================================
    if (in%iHeat) then
        do i = 1,in%Nparts
            call CyclotronResonanceNumber(zp(i),kep(i),xip(i),fcurr(i),in,spline_Bz)
        end do
    end if

    ! Push particles adiabatically:
    ! ==========================================================================
    if (in%iPush) then
      !$OMP PARALLEL DO PRIVATE(i)
        do i = 1,in%Nparts
            call MoveParticle(zp(i),kep(i),xip(i),in,spline_Bz,spline_Phi)
        end do
      !$OMP END PARALLEL DO
    end if

    ! Re-inject particles:
    ! =========================================================================
    if (.true.) then
      !$OMP PARALLEL PRIVATE(ecnt1, ecnt2, pcnt1, pcnt2)
          ecnt1 = 0; ecnt2 = 0; pcnt1 = 0; pcnt2 = 0
          !$OMP DO
              do i = 1,in%Nparts
                  if (zp(i) .GE. in%zmax) then
                      call ReinjectParticles(zp(i),kep(i),xip(i),in,ecnt2,pcnt2)
                  else if (zp(i) .LE. in%zmin) then
                      call ReinjectParticles(zp(i),kep(i),xip(i),in,ecnt1,pcnt1)
                  end if
              end do
          !$OMP END DO
          !$OMP CRITICAL
              ecount1(j) = ecount1(j) + ecnt1
              pcount1(j) = pcount1(j) + pcnt1
              ecount2(j) = ecount2(j) + ecnt2
              pcount2(j) = pcount2(j) + pcnt2
          !$OMP END CRITICAL
      !$OMP END PARALLEL
    end if

    ! Apply Coulomb collision operator:
    ! =========================================================================
    if (in%iColl) then
    		!$OMP PARALLEL PRIVATE(i, id, ecnt, pcnt)
              ecnt = 0; pcnt = 0
              id = OMP_GET_THREAD_NUM()

              in%species_b = 1
          		!$OMP DO
              		do i = 1,in%Nparts
                    !if (id .EQ. 0) write(*,*) "Thread", id, " at i = ", i
              			call collisionOperator(zp(i),kep(i),xip(i),ecnt,pcnt,in)
              		end do
          		!$OMP END DO

              in%species_b = 2
          		!$OMP DO
              		do i = 1,in%Nparts
              		    call collisionOperator(zp(i),kep(i),xip(i),ecnt,pcnt,in)
              		end do
          		!$OMP END DO

          		!$OMP CRITICAL
              		ecount4(j) = ecount4(j) + ecnt
              		pcount4(j) = pcount4(j) + pcnt
          		!$OMP END CRITICAL

    		!$OMP END PARALLEL
    end if

    ! Apply RF heating operator:
    ! =========================================================================
    if (in%iHeat) then
      !$OMP PARALLEL PRIVATE(i, ecnt, pcnt, df, fnew)
          ecnt = 0; pcnt = 0; df = 0;
          !$OMP DO
              do i = 1,in%Nparts
                      call CyclotronResonanceNumber(zp(i),kep(i),xip(i),fnew(i),in,spline_Bz)
                      df = dsign(1.d0,fcurr(i)*fnew(i))
                      if (df .LT. 0 .AND. zp(i) .GT. in%zRes1 .AND. zp(i) .LT. in%zRes2)  then
                        call RFHeatingOperator(zp(i),kep(i),xip(i),ecnt,pcnt,in,spline_Bz,spline_ddBz,spline_Phi)
                      end if
              end do
          !$OMP END DO
          !$OMP CRITICAL
              ecount3(j) = ecount3(j) + ecnt
              pcount3(j) = pcount3(j) + pcnt
          !$OMP END CRITICAL
      !$OMP END PARALLEL
    end if

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
end do TimeStepping

! Record end time:
! =========================================================================
in%tSimTime = tp
oend = OMP_GET_WTIME()

in%tComputeTime = oend-ostart
print *, 'Reached End of Program, Computational time [s] = ', in%tComputeTime

! Save data:
! ==============================================================================
if (in%iSave) then
    ! Create new directory to save output data:
    ! --------------------------------------------------------------------------
    dir1 = trim(in%rootDir)//'/OutputFiles'
    command = 'mkdir '//dir1
    call system(command)

    n_mpwd = lEN_TRIM(xpSelector)-3
    dir0 = xpSelector
    dir0 = dir0(1:n_mpwd)
    dir1 = trim(in%rootDir)//'/OutputFiles/'//trim(dir0)
    command = 'mkdir '//dir1
    call system(command,STATUS)
    dir1 = trim(dir1)//'/'//trim(in%fileDescriptor)
    command = 'mkdir '// dir1
    call system(command,STATUS)
    call getcwd(mpwd)

    ! Saving zp_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'zp.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) zp_hist
    close(unit=8)

    ! Saving kep_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'kep.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) kep_hist
    close(unit=8)

    ! Saving xip_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'xip.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) xip_hist
    close(unit=8)

    ! Saving t_hist to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'tp.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) t_hist
    close(unit=8)

    ! Saving pcount to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'pcount1.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount1
    close(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount2.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount2
    close(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount3.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount3
    close(unit=8)
    fileName = trim(trim(dir1)//'/'//'pcount4.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount4
    close(unit=8)

    ! Saving ecount to file:
    ! --------------------------------------------------------------------------
    fileName = trim(trim(dir1)//'/'//'ecount1.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount1
    close(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount2.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount2
    close(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount3.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount3
    close(unit=8)
    fileName = trim(trim(dir1)//'/'//'ecount4.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount4
    close(unit=8)

    ! Write output data:
    fileName = trim(trim(dir1)//'/'//'metadata.txt')
    open(unit=8,file=fileName,form="formatted",status="unknown")
    write(8,NML = in_nml)
    close(unit=8)

    ! Copy inputdata to output file:
    dir0 = trim(in%rootDir)//'/InputFiles/'//trim(xpSelector)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))
    call system(command)

    ! copy magnetic field data:
    dir0 = trim(in%rootDir)//'/BfieldData'//trim(in%BFieldFile)
    command = 'cp '//trim(trim(dir0)//' '//trim(dir1))//'/Bfield.txt'
    WRITE(*,*) 'command: ', command 
    call system(command)

    ! create text file with commit Hash:
    n_mpwd = lEN_TRIM(xpSelector)-3
    dir0 = xpSelector
    dir0 = dir0(1:n_mpwd)
    dir0 = trim(in%rootDir)//'/OutputFiles/'//trim(dir0)//'/'//trim(in%fileDescriptor)
    fileName = trim(dir0)//'/commitHash.txt'
    command = 'git log --oneline -1 > '//trim(fileName)
    call system(command)

end if

End PROGRAM
