PROGRAM LinearFokkerPlanckSolver_1D2V
! ===========================================================================
! PURPOSE:
! Created by JF Caneses Marin on 2019-09-05
! Modified: 2019_10_23, Correct error in the RF energy kick (JFCM)
! Modified: 2020_03_05, OMP constructs added (JFCM)
! ===========================================================================

! ===========================================================================
! Define modules to use
USE local
USE spline_fits
USE PhysicalConstants
USE dataTYP
USE OMP_LIB

! ===========================================================================
! Define local variables
! ===========================================================================
IMPLICIT NONE
TYPE(inTYP)  :: in
TYPE(derTYP) :: der
TYPE(splTYP) :: spline_Bz
TYPE(splTYP) :: spline_ddBz
TYPE(splTYP) :: spline_Phi
TYPE(splTYP) :: spline_j0
TYPE(splTYP) :: spline_j1
TYPE(spltestTYP) :: spline_Test
INTEGER(i4) :: i,j,k                                                          ! Indices for do loops
INTEGER(i4) :: seed_size                                                      ! Random number generator variable
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: seed                                ! Store the random num gen seed
INTEGER(i4) :: jsize
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: jrng                                ! Indices of time steps to save
INTEGER(i4) :: id
REAL(r8) :: tstart, tend, tComputeTime, tSimTime                              ! Variables to hold cpu time at start and end of computation
REAL(r8) :: df
REAl(r8) :: tp                                                      ! Hold simulation time
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: xip, zp, kep                 ! Particle position (zp), kinetic energy (KEp), pitch angle (Xip)
REAL(r8), DIMENSION(:,:), ALLOCATABLE :: zp_hist, kep_hist, xip_hist
REAl(r8), DIMENSION(:)  , ALLOCATABLE :: t_hist
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: pcount1, pcount2, pcount3, pcount4   ! Count the number of particles incident on (1) dump, (2) target, (3) EBW resonance
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: ecount1, ecount2, ecount3, ecount4   ! Record the total energy of particle incident on (1) dump, (2) target, (3) EBW resonance
REAL(r8), DIMENSION(:)  , ALLOCATABLE :: fcurr, fnew
REAL(r8) :: ecnt, ecnt1, ecnt2
REAL(r8) :: pcnt, pcnt1, pcnt2
CHARACTER*150 :: command, mpwd
CHARACTER*250 :: fileName                                            ! Are used to name the output files

!INTEGER, EXTERNAL :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
!INTEGER, EXTERNAL :: OMP_SET_NUM_THREADS, OMP_GET_WTIME

! Create input and output namelists from the user-defined structures:
namelist/in_nml/in

! Record start time:
call cpu_time(tstart)

! Read input data into in structure:
fileName = "data.in"
open(unit=4,file=fileName,status='old',form='formatted')
read(4,in_nml)
close(unit=4)

if (in%species_a .eq. 1) then
    der%q   = -e_c
    der%m_t = m_e
else
    der%q   = +in%Zion*e_c
    der%m_t = in%Aion*m_p
end if

print *, 'fileDescriptor     ', in%fileDescriptor
print *, 'Number of particles', in%Nparts
print *, 'Number of steps    ', in%Nsteps
print *, 'dt [ns]            ', in%dt*1E+9
print *, 'iPush              ', in%iPush
print *, 'iDrag              ', in%iDrag
print *, 'iColl              ', in%iColl
print *, 'iHeat              ', in%iHeat
print *, 'iSave              ', in%iSave
print *, 'elevel             ', in%elevel
print *, 'zTarget [m]        ', in%zmax
print *, 'zDump [m]          ', in%zmin
print *, 'zp_init [m]        ', in%zp_init
print *, 'B field file       ', in%BFieldFile
print *, 'Ew                 ', in%Ew
print *, 'Te0                ', in%Te0
print *, 'ne0                ', in%ne0

if (in%CollOperType .EQ. 1) print *, 'Boozer-Only collision operator'
if (in%CollOperType .EQ. 2) print *, 'Boozer-Kim collision operator'

! ===========================================================================
call InitSpline(spline_Bz  ,in%nz,0._8,0._8,1,0._8)
call InitSpline(spline_ddBz,in%nz,0._8,0._8,1,0._8)
call InitSpline(spline_Phi ,in%nz,0._8,0._8,1,0._8)
call InitSpline(spline_j0  ,in%nz,0._8,0._8,1,3._8)
call InitSpline(spline_j1  ,in%nz,0._8,0._8,1,3._8)
call InitSplineTest(spline_Test,in%nz)
call InitDer(der,in)

! Allocate memory to "allocatable" variables
ALLOCATE(zp(in%Nparts), kep(in%Nparts), xip(in%Nparts))
ALLOCATE(pcount1(in%Nsteps),pcount2(in%Nsteps),pcount3(in%Nsteps),pcount4(in%Nsteps))
ALLOCATE(ecount1(in%Nsteps),ecount2(in%Nsteps),ecount3(in%Nsteps),ecount4(in%Nsteps))
ALLOCATE(fcurr(in%Nparts),fnew(in%Nparts))

! Define variables for the saving process:
jsize = (in%jend-in%jstart+1)/in%jincr
ALLOCATE(jrng(jsize))
ALLOCATE(zp_hist(in%Nparts,jsize),kep_hist(in%Nparts,jsize),xip_hist(in%Nparts,jsize),t_hist(jsize))

! Create array with the indices of the time steps to save:
jrng = (/ (j, j=in%jstart, in%jend, in%jincr) /)

! ===========================================================================
! Magnetic field data:
fileName = trim(adjustl(in%BFieldFile))
fileName = trim(adjustl(in%BFieldFileDir))//fileName
fileName = trim(adjustl(in%rootDir))//fileName
CALL ReadSpline(spline_Bz,fileName)
CALL ComputeSpline(spline_Bz)

! ===========================================================================
! Second derivative of the magnetic field:
spline_ddBz%x = spline_Bz%x
spline_ddBz%y = spline_Bz%yp
CALL ComputeSpline(spline_ddBz)

! ===========================================================================
! Predefined electric potential:
spline_Phi%x = spline_Bz%x
spline_Phi%y = 0
if (in%iPotential) then
  call PotentialProfile(spline_Phi,in)
end if
CALL ComputeSpline(spline_Phi)

! ===========================================================================
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

! ===========================================================================
! Initialize the random number generator
call random_seed(size=seed_size)
ALLOCATE(seed(seed_size))
call random_seed(get=seed)
seed = 314159565
call random_seed(put=seed)

! ==============================================================================
! Initialize particle position zp, kinetic energy kep and pitch angle xip
kep = 0.; xip = 0.; zp = 0.;
call loadParticles(zp,kep,xip,in,der)

fileName = "LoadParticles.dat"
open(unit=8,file=fileName,form="formatted",status="unknown")
do i = 1,in%nz
    write(8,*) zp(i), kep(i), xip(i)
end do
close(unit=8)

! Initialize the time array
tp = 0
! Initialize leak current diagnostics
pcount1 = 0; pcount2 = 0; pcount3 = 0; pcount4 = 0
! Initialize the energy leak diagnotics
ecount1 = 0; ecount2 = 0; ecount3 = 0; ecount4 = 0

! Set the number of threads:
call OMP_SET_NUM_THREADS(in%threads_request)
!$OMP PARALLEL PRIVATE(id)
    id = OMP_GET_THREAD_NUM()
    in%threads_given = OMP_GET_NUM_THREADS()
    if (id .EQ. 0) write(*,*) "number of threads given: ", in%threads_given
!$OMP END PARALLEL

! ==============================================================================
! BEGIN TIME STEPPING
TimeStepping: do j = 1,in%Nsteps
    ! Calculate Cyclotron resonance number
    ! =========================================================================

    if (in%iHeat) then
        do i = 1,in%Nparts
            call CyclotronResonanceNumber(zp(i),kep(i),xip(i),fcurr(i),in,der,spline_Bz)
        end do
    end if

    ! =========================================================================
    ! PUSH PARTICLES ADIABATICALLY
        !WRITE(*,*) 'kep(3), before', kep(3)
    if (in%iPush) then
      !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
        do i = 1,in%Nparts
            call MoveParticle(zp(i),kep(i),xip(i),in,der,spline_Bz,spline_Phi)
        end do
      !$OMP END PARALLEL DO
    end if
    !WRITE(*,*) 'kep(3), after', kep(3)

    ! =========================================================================
    ! RE-INJECT PARTICLES
    if (.true.) then
      !$OMP PARALLEL PRIVATE(ecnt1, ecnt2, pcnt1, pcnt2)
          ecnt1 = 0; ecnt2 = 0; pcnt1 = 0; pcnt2 = 0
          !$OMP DO
              do i = 1,in%Nparts
                  if (zp(i) .GE. in%zmax) then
                      call ReinjectParticles(zp(i),kep(i),xip(i),in,der,ecnt2,pcnt2)
                  else if (zp(i) .LE. in%zmin) then
                      call ReinjectParticles(zp(i),kep(i),xip(i),in,der,ecnt1,pcnt1)
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

    ! =========================================================================
    ! APPLY COULOMB COLLISION OPERATOR
    !write(*,*) "time", j
    if (in%iColl) then
    		!$OMP PARALLEL PRIVATE(i, id, ecnt, pcnt)
              ecnt = 0; pcnt = 0
              id = OMP_GET_THREAD_NUM()

              der%species_b = 1 ! test particle on field electrons
          		!$OMP DO SCHEDULE(STATIC)
              		do i = 1,in%Nparts
                    !if (id .EQ. 0) write(*,*) "Thread", id, " at i = ", i
              			call collisionOperator(zp(i),kep(i),xip(i),ecnt,pcnt,in,der)
              		end do
          		!$OMP END DO

              der%species_b = 2 ! test particle on field ions
          		!$OMP DO
              		do i = 1,in%Nparts
              		    call collisionOperator(zp(i),kep(i),xip(i),ecnt,pcnt,in,der)
              		end do
          		!$OMP END DO

          		!$OMP CRITICAL
              		ecount4(j) = ecount4(j) + ecnt
              		pcount4(j) = pcount4(j) + pcnt
          		!$OMP END CRITICAL

    		!$OMP END PARALLEL
    end if

    ! =========================================================================
    ! APPLY RF HEATING OPERATOR
    ! 	Modify: kep(i), xip(i)
    ! 	Conserved: z(i)
    if (in%iHeat) then
      !$OMP PARALLEL PRIVATE(i, ecnt, pcnt, df, fnew)
          ecnt = 0; pcnt = 0; df = 0;
          !$OMP DO SCHEDULE(STATIC)
              do i = 1,in%Nparts
                      call CyclotronResonanceNumber(zp(i),kep(i),xip(i),fnew(i),in,der,spline_Bz)
                      df = dsign(1.d0,fcurr(i)*fnew(i))
                      if (df .LT. 0 .AND. zp(i) .GT. in%zRes1 .AND. zp(i) .LT. in%zRes2)  then
                        call RFHeatingOperator(zp(i),kep(i),xip(i),ecnt,pcnt,in,der,spline_Bz,spline_ddBz,spline_Phi)
                      end if
              end do
          !$OMP END DO
          !$OMP CRITICAL
              ecount3(j) = ecount3(j) + ecnt
              pcount3(j) = pcount3(j) + pcnt
          !$OMP END CRITICAL
      !$OMP END PARALLEL
    end if

    ! =========================================================================
    ! Update time array
    tp = tp + in%dt

    ! =====================================================================
    ! SELECT DATA TO BE SAVED TO OUTPUT FILE
    ! Check if data is to be saved
    if (in%iSave) then
        do k = 1,jsize
            if (j .EQ. jrng(k)) then
                t_hist(k) = tp
				        !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
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
tSimTime = tp
call cpu_time(tend)
tComputeTime = (tend-tstart)/in%threads_given
print *, 'Reached End of Program, Computational time = ', tComputeTime
! Save data:
! =========================================================================
if (in%iSave) then

    ! Create new directory to save output data:
    ! ---------------------------------------------------------------------
	command = 'mkdir '// trim(in%fileDescriptor)
	write(*,*) command
    call system(command)
    call getcwd(mpwd)

    ! Saving zp_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'zp.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) zp_hist
    close(unit=8)


    ! Saving kep_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'kep.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) kep_hist
    close(unit=8)

    ! Saving xip_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'xip.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) xip_hist
    close(unit=8)

    ! Saving t_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'tp.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) t_hist
    close(unit=8)

    ! Saving pcount to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'pcount1.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount1
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'pcount2.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount2
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'pcount3.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount3
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'pcount4.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount4
    close(unit=8)

    ! Saving ecount to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'ecount1.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount1
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'ecount2.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount2
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'ecount3.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount3
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'ecount4.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount4
    close(unit=8)

    ! Write output data:
    fileName = trim(trim(mpwd)//'/'//trim(in%fileDescriptor)//'/'//'data.out')
    open(unit=8,file=fileName,form="formatted",status="unknown")
    write(8,NML = in_nml)
    close(unit=8)

end if

End PROGRAM
