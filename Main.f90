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

#ifdef mpi
    include 'mpif.h'
#endif        

    ! Variables
    type (structparamstype) :: params
    type (gridsType) :: grids
    !real (kind=rk) :: minInc(Tperiods)
    real (kind=rk) :: V(Tperiods+1, numPointsType,numPointsA, numAIME, numPointsY)
    real (kind=rk) :: policyA1(Tperiods, numPointsType, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: policyC(Tperiods, numPointsType, numPointsA, numAIME, numPointsY)
    integer :: policyL(Tperiods, numPointsType, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: EV(Tperiods+1, numPointsType, numPointsA, numAIME, numPointsY)

    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    real (kind=rk) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    real (kind=rk) :: start, finish, moments(2,24),weights(2,24), y(dimEstimation+1), p(dimEstimation+1,dimEstimation)
    integer :: action, ios, requiredl, typeSim
    INTEGER(kind=4) :: iter

#ifdef mpi
    integer :: provided
#endif       
    call cpu_time(start)
    
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
    
    params%r= 0.02

    params%mu = 0
    params%sigma = 0.0922
    params%rho = 0.96
    params%hrsWrk = 0.3159
    !params%delta(3) = 9.083298
    !params%delta(2)=0.029333
    !params%delta(1)=-0.00033

    !Earning type 1 - uneducated no db
    params%delta(1,3) = 9.083298
    params%delta(1,2)=0.029333
    params%delta(1,1)=-0.00023

    !Earning type 2 - uneducated db
    if (numPointsType .GE. 2) then
        params%delta(2,3) = 8.5256042
        params%delta(2,2)= 0.061582208
        params%delta(2,1)= -0.0005901
    end if
    
    !Earning type 3 - educated no db
    if (numPointsType .GE. 3) then
        params%delta(3,3) = 7.906901
        params%delta(3,2)= 0.094386
        params%delta(3,1)= -0.00091
    end if
    !Earning type 4 - educated db
    if (numPointsType .GE. 4) then
        params%delta(4,3) = 7.392151
        params%delta(4,2)= 0.121426
        params%delta(4,1)= -0.00117
    end if

    params%spouseInc(1) = 5121!6235.8998884204320
    params%spouseInc(2) = 5282
    params%spouseInc(3) = 6840
    params%spouseInc(4) = 7684
    
    params%tol = 1e-10
    params%minCons = 1e-5

    params%nu =   0.372042538693115
    params%beta = 0.999!0.976948640333028
    params%gamma = 3 !2.46102724400491
    params%db(1) = 0.764846423121866
    params%db(2) =  -5.472985265008665E-006 !-5.472985265008665E-005
    params%thetab = 0.0375
    params%k=650000


    call setupMisc(params,grids)
    if (fullLifeCycle) then
        params%startA =0
    else
        params%startA = maxval(grids%initialAssets)
    end if


    action =2
    if (action .EQ. 1) then
        params%nu = 0.268145043046353 !0.279542353537685 !0.281392371818158 ! 0.249813174520667 ! 0.371590333431661 !0.400970471984035 !0.379444754253488
        params%beta = 0.886870656649138! 0.917822713705556!0.915901849357468 !0.907953667733027  !  0.981318536670500 !0.974287833949847  !0.999958017173274
        params%gamma =  2.16221278404756! 2.25411606955360! 2.26903386602418 !1.65249116180330 !1.73187978845077 ! 2.65238277512422  !2.50999213311961
        params%db(1) = 0.551253586443696 !0.574684220151412 ! 0.578487494679671  ! 0.513566862288025  !0.883861989502770 ! 0.824316602045162 !0.780063950700460
        params%db(2) = -3.944586276588851E-006! -4.112248054123012E-006!  -4.139463186468392E-007  !-3.674912562169277E-005 !-4.939458831206524E-005! -5.898533683743610E-005 !-5.581876523249634E-005
       
       
       
       
 
        
        do typeSim = 1, numPointsType
            call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
        end do
        
        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, .TRUE. )
        !simulate
        if (rank == 0) then
            !call cpu_time(starsim)
            call simWithUncer(params, grids, policyA1,policyL,EV, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME )
            !call cpu_time(finish)
            !print '("Time = ",f11.3," seconds.")',finish-starsim
            call writetofile(ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
            call writetofileByType(ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
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

        call amoeba(p,y, 0.07_rk,gmm_criteria,iter) !0.0001_rk !0.001_rk!

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
            call writetofileByType(ypath, cpath, apath, vpath, lpath, grids%Simy,AIME)
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
    if (minval(control(1:4)) < 0 .or.  control(5)> 0 .OR. control(6)< 0 )  then !.OR. control(6)< 0
        gmm_criteria = huge(gmm_criteria)
        return
    end if
    params%nu = control(1)
    params%beta = control(2)
    params%gamma = control(3)
    params%db(1)= control(4)
    params%db(2)= control(5)
    params%thetab = control(6)
    write (*,*) "Calcualting GMM for",  control(1),  control(2),  control(3),  control(4),  control(5),  control(6)
    gmm_criteria = gmm(params,grids,moments,weights) !*-1.0

    end function

    end program Console1

