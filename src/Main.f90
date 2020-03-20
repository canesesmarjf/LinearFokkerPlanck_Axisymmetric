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
USE ParticlePusher
USE plasma_params
USE collision_data
USE rf_heating_data
USE dataTYP
USE OMP_LIB

! ===========================================================================

! Define local variables
! ===========================================================================
IMPLICIT NONE
TYPE(inTYP)  :: in

REAL(r8) :: tstart, tend, tComputeTime, tSimTime                              ! Variables to hold cpu time at start and end of computation
INTEGER(i4) :: i,j,k                                                          ! Indices for do loops
REAL(r8) :: curv2, curvd                                                      ! Declare functions from fitpack
INTEGER(i4) :: seed_size                                                      ! Random number generator variable
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: seed                                ! Store the random num gen seed

REAL(r8), DIMENSION(:), ALLOCATABLE :: pcount1, pcount2, pcount3    ! Count the number of particles incident on (1) dump, (2) target, (3) EBW resonance
REAL(r8), DIMENSION(:), ALLOCATABLE :: ecount1, ecount2, ecount3    ! Record the total energy of particle incident on (1) dump, (2) target, (3) EBW resonance
REAL(r8), DIMENSION(:), ALLOCATABLE :: ecount4                      ! Record the total energy dissipated by collisional slowing down within a time step dt
REAL(r8), DIMENSION(:), ALLOCATABLE :: pcount4                      ! Record the total number of fast particles involved in the slowing down dissipated power within a time step dt

REAL(r8) :: ecnt, ecnt1, ecnt2
REAL(r8) :: pcnt, pcnt1, pcnt2
INTEGER(i4) :: jsize                                                          ! Total umber of time steps to save
INTEGER(i4), DIMENSION(:), ALLOCATABLE :: jrng                                ! Indices of time steps to save
REAL(r8) :: df
CHARACTER*100 :: besselj_data                                                 ! Name of file containing the bessel function data

! Create local variables to hold some of the data from InputFile
LOGICAL :: iColl, iHeat, iSave, iPush, iPotential
INTEGER(i4) :: zp_InitType, kep_InitType, xip_InitType                                      !
CHARACTER*150 :: BFieldFile, BFieldFileDir, B_data, rootDir
CHARACTER*250 :: fileExt, fileName                                            ! Are used to name the output files
INTEGER(i4) :: jstart, jend, jincr                                            ! To define time steps to save
REAL(r8) :: zTarget, zDump
CHARACTER*150 :: fileDescriptor
CHARACTER*150 :: command, mpwd

!INTEGER(i4) ::  threads_request

! Create input and output namelists from the user-defined structures:
namelist/in_nml/in

! Create namelist:
! ==================================================================================================
! Namelist for input data:
! ------------------------
namelist /indata/ BFieldFile, BFieldFileDir, nz, Nparts, Nsteps, iColl, &
                dt, ne0, Te0, ti0, Zeff, Zion, species_a, zp_InitType, &
                zp_init, zp_init_std, iHeat, Aion, s1, s2, s3, phi1, phi2, phi3, jstart, jend, jincr, &
                f_RF, zRes1, zRes2, kper, kpar, Ew, n_harmonic, kep_init, iDrag, iSave, iPush, elevel, CollOperType, &
                zTarget, zDump, fileDescriptor, iPotential, kep_InitType, xip_InitType, xip_init, rootDir, threads_request

! Namelist for output metadata file:
! ----------------------------------
namelist /metadata/ fileDescriptor, &
                    ne0, Te0, Ti0, Aion, Zeff, Zion, species_a, zDump, zTarget, &
                    Nparts, Nsteps, dt, tComputeTime, tSimTime, threads_request, &
                    jstart, jend, jincr, &
                    f_RF, n_harmonic, Ew, kper, kpar, zRes1, zRes2, &
                    CollOperType, elevel, &
                    iSave, iPush, iColl, iHeat, iPotential, iDrag, &
                    zp_InitType, zp_init, zp_init_std, &
                    kep_InitType, kep_init, &
                    xip_InitType, xip_init, &
                    BFieldFileDir, BFieldFile, nz, &
                    s1, s2, s3, phi1, phi2, phi3

! Record start time:
call cpu_time(tstart)

! Read input data:
fileName = "data.in"
open(unit=4,file=fileName,status='old',form='formatted')
read(4,in_nml)
close(unit=4)

! ===========================================================================
! Read input files
open(unit=4,file='inputfile.in',status='old',form='formatted')
read(4,indata)
close(unit=4)

if (in%species_a .eq. 1) then
    q = -e_c
    m_t =m_e
    print *, 'Test particles: Electrons'
else
    q = +in%Zion*e_c
    m_t = in%Aion*m_p
    print *, 'Test particles: Ions'
end if

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
print *, 'zp_init [m]       ', in%zp_init
print *, 'B field file       ', in%BFieldFile
print *, 'Ew                 ', in%Ew
print *, 'Te0                ', in%Te0
print *, 'ne0                 ', in%ne0

if (in%CollOperType .EQ. 1) print *, 'Boozer-Only collision operator'
if (in%CollOperType .EQ. 2) print *, 'Boozer-Kim collision operator'

! ===========================================================================
! Allocate memory to "allocatable" variables
ALLOCATE(zp(in%Nparts), kep(in%Nparts), xip(in%Nparts))
ALLOCATE(pcount1(in%Nsteps),pcount2(in%Nsteps),pcount3(in%Nsteps),pcount4(in%Nsteps))
ALLOCATE(ecount1(in%Nsteps),ecount2(in%Nsteps),ecount3(in%Nsteps),ecount4(in%Nsteps))
! Variables for spline fits
ALLOCATE(z_Ref(in%nz), B_Ref(in%nz), Phi_Ref(in%nz), ddB_Ref(in%nz))
ALLOCATE(b_spl(in%nz), b_temp(in%nz), phi_spl(in%nz), phi_temp(in%nz), ddb_spl(in%nz), ddb_temp(in%nz))
ALLOCATE(zz(in%nz), b1(in%nz), ddb1(in%nz))
! Variables for RF operator
ALLOCATE(fcurr(in%Nparts),fnew(in%Nparts))
! Define variables for the saving process
jsize = (in%jend-in%jstart+1)/in%jincr
ALLOCATE(jrng(jsize))
ALLOCATE(zp_hist(in%Nparts,jsize),kep_hist(in%Nparts,jsize),xip_hist(in%Nparts,jsize),t_hist(jsize))

! Create array with the indices of the time steps to save
jrng = (/ (j, j=jstart, jend, jincr) /)

! ===========================================================================
! Read Bessel function data
besselj_data = trim(adjustl("besselj01_0_to_40.txt"))
open(unit=8,file=besselj_data,status="old")
do i=1,501
    read(8,*) x_j_ref(i),j0_ref(i), j1_ref(i)
end do
close(unit=8)
slp1 = 0.; slpn = 0.; sigma = 1.;islpsw = 3
call curv1(501,x_j_ref,j0_ref,slp1,slpn,islpsw,j0_spl,j0_temp,sigma,ierr1)
call curv1(501,x_j_ref,j1_ref,slp1,slpn,islpsw,j1_spl,j1_temp,sigma,ierr1)

! ===========================================================================
! Read magnetic field data
B_data = trim(adjustl(in%BFieldFile))
B_data = trim(adjustl(in%BFieldFileDir))//B_data
B_data = trim(adjustl(in%rootDir))//B_data
open(unit=8,file=B_data,status="old")
do i=1,501
    read(8,*) z_Ref(i),B_Ref(i)
end do
close(unit=8)
slp1 = 0.; slpn = 0.; sigma = 1.; islpsw = 3
call curv1(in%nz,z_Ref,B_Ref  ,slp1,slpn,islpsw,b_spl  ,b_temp  ,sigma,ierr1)

! nz,z_Ref,B_Ref,slp1,slpn,islpsw and sigma are unaltered.
! Outputs are b_spl, b_temp, ierr1
! b_spl represents the 2nd spatial derivative of B_ref
! b_spl is to be used to calculate the second derivative of Omega for the interaction time

! Second spatial derivative of the magnetic field
ddB_Ref = b_spl
call curv1(in%nz,z_Ref,ddB_Ref,slp1,slpn,islpsw,ddb_spl,ddb_temp,sigma,ierr2)
! need to confirm correctness of second derivative

! ===========================================================================
! Setup predefined electric potential
call PotentialProfile(in%iPotential)

! Electric potential
slp1 = 0.; slpn = 0.; sigma = 1. ;islpsw = 3
call curv1(in%nz,z_Ref,Phi_Ref,slp1,slpn,islpsw,phi_spl,phi_temp,sigma,ierr1)

if (.false.) then
! Test the spline
    do i=1,500
        zz(i) = z_Ref(i) + 1.e-3*(-z_Ref(i) + z_Ref(i+1))/2.
        b1(i)   = curv2(zz(i),in%nz,z_Ref,B_Ref  ,b_spl  ,sigma)
        ddb1(i) = curv2(zz(i),in%nz,z_Ref,ddB_Ref,ddb_spl,sigma)
    end do
    ! Save data to file to test spline
    if (.true.) then
    fileName = "B_spline.dat"
        open(unit=8,file=fileName,form="formatted",status="unknown")
        do j = 1,500
            write(8,*) zz(j), b1(j), ddb1(j)
        end do
        close(unit=8)
    end if
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
call loadParticles(in)

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
            call CyclotronResonanceNumber(zp(i),kep(i),xip(i),fcurr(i))
        end do
    end if

    ! =========================================================================
    ! PUSH PARTICLES ADIABATICALLY
    if (in%iPush) then
      !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
        do i = 1,in%Nparts
            call MoveParticle(zp(i),kep(i),xip(i))
        end do
      !$OMP END PARALLEL DO
    end if

    ! =========================================================================
    ! RE-INJECT PARTICLES
    if (.true.) then
      !$OMP PARALLEL PRIVATE(ecnt1, ecnt2, pcnt1, pcnt2)
          ecnt1 = 0; ecnt2 = 0; pcnt1 = 0; pcnt2 = 0
          !$OMP DO
              do i = 1,in%Nparts
                  if (zp(i) .GE. in%zmax) then
                      call ReinjectParticles(zp(i),kep(i),xip(i),ecnt2,pcnt2)
                  else if (zp(i) .LE. in%zmin) then
                      call ReinjectParticles(zp(i),kep(i),xip(i),ecnt1,pcnt1)
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

              species_b = 1 ! test particle on field electrons
          		!$OMP DO SCHEDULE(STATIC)
              		do i = 1,in%Nparts
                    !if (id .EQ. 0) write(*,*) "Thread", id, " at i = ", i
              			call collisionOperator(zp(i),kep(i),xip(i),ecnt,pcnt)
              		end do
          		!$OMP END DO

              species_b = 2 ! test particle on field ions
          		!$OMP DO
              		do i = 1,in%Nparts
              		    call collisionOperator(zp(i),kep(i),xip(i),ecnt,pcnt)
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
                      call CyclotronResonanceNumber(zp(i),kep(i),xip(i),fnew(i))
                      df = dsign(1.d0,fcurr(i)*fnew(i))
                      if (df .LT. 0 .AND. zp(i) .GT. in%zRes1 .AND. zp(i) .LT. in%zRes2)  then
                        call RFHeatingOperator(zp(i),kep(i),xip(i),ecnt,pcnt)
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
	command = 'mkdir '// trim(fileDescriptor)
	write(*,*) command
    call system(command)
    call getcwd(mpwd)

    ! Saving zp_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'zp.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) zp_hist
    close(unit=8)


    ! Saving kep_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'kep.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) kep_hist
    close(unit=8)

    ! Saving xip_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'xip.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) xip_hist
    close(unit=8)

    ! Saving t_hist to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'tp.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) t_hist
    close(unit=8)

    ! Saving pcount to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'pcount1.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount1
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'pcount2.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount2
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'pcount3.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount3
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'pcount4.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) pcount4
    close(unit=8)

    ! Saving ecount to file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'ecount1.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount1
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'ecount2.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount2
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'ecount3.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount3
    close(unit=8)
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'ecount4.out')
    open(unit=8,file=fileName,form="unformatted",status="unknown")
    write(8) ecount4
    close(unit=8)

    ! Create Metadata file
    ! ---------------------------------------------------------------------
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'Metadata.out')
    open(unit=8,file=fileName,form="formatted",status="unknown")
    write(8,NML = metadata)
    close(unit=8)

    ! Write output data:
    fileName = trim(trim(mpwd)//'/'//trim(fileDescriptor)//'/'//'data.out')
    open(unit=8,file=fileName,form="formatted",status="unknown")
    write(8,NML = in_nml)
    close(unit=8)

end if

End PROGRAM
