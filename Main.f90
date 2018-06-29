    !  Console1.f90
    !
    !  FUNCTIONS:
    !  Console1 - Entry point of console application.
    !

    !****************************************************************************
    !
    !  PROGRAM: Console1
    !
    !  PURPOSE:  Entry point for the console application.
    !
    !****************************************************************************
    !#include "fintrf.h"
    program Console1

    use Header
    use routines

    implicit none

#ifdef mpi  !mpi
    include 'mpif.h'
#endif        

    ! Variables
    type (structparamstype) :: params
    type (gridsType) :: grids
    !real (kind=rk) :: minInc(Tperiods)
    real (kind=rk) :: V(Tperiods+1, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: policyA1(Tperiods, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: policyC(Tperiods, numPointsA, numAIME, numPointsY)
    integer :: policyL(Tperiods, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: EV(Tperiods+1, numPointsA, numAIME, numPointsY)

    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    real (kind=rk) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    real (kind=rk) :: start, finish, moments(2,24),weights(2,24), y(dimEstimation+1), p(dimEstimation+1,dimEstimation)
    integer :: action, ios, requiredl
    INTEGER(kind=4) :: iter
    
#ifdef mpi  !mpi    
    integer :: provided
#endif       
    call cpu_time(start)

    params%r= 0.02
    if (fullLifeCycle) then
        params%startA =0
    else
        params%startA = maxval(grids%initialAssets)
    end if
    params%mu = 0
    params%sigma = 0.0922
    params%rho = 0.96
    params%hrsWrk = 0.3159
    params%delta(3) = 9.083298
    params%delta(2)=0.029333
    params%delta(1)=-0.00033
    params%spouseInc = 6235.8998884204320
    params%tol = 1e-10
    params%minCons = 1e-5

    params%nu =   0.372042538693115
    params%beta = 0.976948640333028
    params%gamma = 2.46102724400491
    params%db(1) = 0.764846423121866
    params%db(2) =  -5.472985265008665E-005
    params%thetab = 0.0375
    params%k=650000

    
    call setupMisc(params,grids)
    
#ifdef mpi    
    call MPI_Init_thread(MPI_THREAD_MULTIPLE,provided,ierror)!mpi_init
    if (ierror.ne.0) stop 'mpi problem171'
    call mpi_comm_rank(mpi_comm_world, rank, ierror)
    call mpi_comm_size(mpi_comm_world, procSize, ierror)
    if (ierror.ne.0) stop 'mpi problem172'
    if (rank.eq.0) write (*, *) 'Using MPI in solution. Using', procSize, 'cores'
    if (rank.eq.0) write(*,*)
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem173'
#else
    rank = 0
    procSize = 1
#endif

    action =1
    if (action .EQ. 1) then
        params%nu = 0.400970471984035 !0.379444754253488
        params%beta = 0.974287833949847  !0.999958017173274
        params%gamma = 2.65238277512422  !2.50999213311961
        params%db(1) = 0.824316602045162 !0.780063950700460
        params%db(2) =  -5.898533683743610E-005 !-5.581876523249634E-005    
        
        call getassetgrid( params, grids%maxInc, grids%Agrid)
        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, .TRUE. )
        !simulate
        if (rank == 0) then
            !call cpu_time(starsim)
            call simWithUncer(params, grids, policyA1,policyL,EV, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME )
            !call cpu_time(finish)
            !print '("Time = ",f11.3," seconds.")',finish-starsim
            call writetofile(ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        end if
    else
        if (Tretire+startAge .NE. 60) then
            write (*,*) "Only estimate with SPA = 60"
        end if
#ifdef win
        open (unit = 1001,file='..\\..\\moments\\moments.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file='..\\..\\moments\\assetmom.txt', action='read', IOSTAT = ios)
        open (unit = 1003,file='..\\..\\moments\\Weight.txt', action='read', IOSTAT = ios)
        open (unit = 1004,file='..\\..\\moments\\weightAsset.txt', action='read', IOSTAT = ios)
#else
        open (unit = 1001,file='../moments/moments.txt', action='read', IOSTAT = ios)
        open (unit = 1002,file='../moments/assetmom.txt', action='read', IOSTAT = ios)
#endif
        read (1001, *) moments(1,:)
        read (1002, *) moments(2,:)
        read (1003, *) weights(1,:)
        read (1004, *) weights(2,:)
        close (unit=1001)
        close (unit=1002)
        close (unit=1003)
        close (unit=1004)
        
        if (rank==0) then
            print '("Setting up initial guess for hilling climbing algorithm")'
        end if

        call initialGuess(rank,params,grids,moments,weights,p,y)

        call amoeba(p,y,0.0001_rk,gmm_criteria,iter) !0.001_rk

        if (rank==0) then
            print '("P = ",f6.3)',P(1,:)
            print '("Y = ",f16.3)',Y
#ifdef win
            inquire (iolength=requiredl)  P(1,:)
            open (unit=201, file='..\\out\params.txt', status='unknown',recl=requiredl, action='write')
            write (201, * ) P(1,:)
#else
            inquire (iolength=requiredl)  P
            open (unit=201, file='./out/params.txt', status='unknown',recl=requiredl, action='write')
            write (201, * ) P(1,1)
            write (201, * ) P(1,2)
            write (201, * ) P(1,3)
            write (201, * ) P(1,4)
            write (201, * ) P(1,5)
#endif
            close (unit=201)

            print '("Generating files")'
        end if

        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, .FALSE. )

        if (rank==0) then
            call simWithUncer(params, grids, policyA1,policyL,EV, grids%Simy, cpath, apath, vpath, lpath, ypath ,AIME)
            call writetofile(ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
        end if
#ifdef mpi 
        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'
#endif        
    end if

    if (rank==0) then
        call cpu_time(finish)
        print '("Time = ",f11.3," seconds.")',finish-start
    end if
#ifdef mpi
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
    if (rank.eq.0) then
        call mpi_finalize(ierror)
    end if
    if (ierror.ne.0) stop 'mpi problem190'
#endif

    contains
    function gmm_criteria(control)
    implicit none
    !inputs
    real (kind=rk), intent(in) ::control(:)
    !output
    real (kind=rk) :: gmm_criteria
    if (maxval(control(1:2)) > 1)  then
        gmm_criteria = huge(gmm_criteria)
        return
    end if
    if (minval(control(1:4)) < 0 .or.  control(5)> 0 )  then !.OR. control(6)< 0
        gmm_criteria = huge(gmm_criteria)
        return
    end if
    params%nu = control(1)
    params%beta = control(2)
    params%gamma = control(3)
    params%db(1)= control(4)
    params%db(2)= control(5)
    gmm_criteria = gmm(params,grids,moments,weights) !*-1.0

    end function

    end program Console1

