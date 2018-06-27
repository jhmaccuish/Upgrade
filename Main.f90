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
    real (kind=rk) :: minInc(Tperiods), maxA(Tperiods+1)
    real (kind=rk) :: span, loggrid(numPointsA)
    integer :: ixT, i
    real (kind=rk) :: V(Tperiods+1, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: policyA1(Tperiods, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: policyC(Tperiods, numPointsA, numAIME, numPointsY)
    integer :: policyL(Tperiods, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: EV(Tperiods+1, numPointsA, numAIME, numPointsY)
    real (kind=rk) :: EdU(Tperiods,   numPointsA, numAIME, numPointsY)

    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: yemp(Tperiods, numSims)
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    real (kind=rk) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    real (kind=rk) :: start, finish, moments(2,24), error, test, starsim
    integer :: action, ios, inter, requiredl, provided
    INTEGER(kind=4) :: iter
    REAL(kind=rk) :: y(dimEstimation+1)
    REAL(kind=rk) :: p(dimEstimation+1,dimEstimation), temp(85)

    params%system=1

    call cpu_time(start)
    if (params%system == 1 ) then
        open (unit = 1001,file='..\\..\\moments\\InitialAssets.txt', action='read', IOSTAT = ios)
    else
        open (unit = 1001,file='../moments/InitialAssets.txt', action='read', IOSTAT = ios)
    end if
    read (1001, *,  IOSTAT = ios) grids%initialAssets
    close (unit=1001)

    temp = (/-0.01374, -0.01069, -0.00767, -0.00467, -0.00169, 0.00126, 0.00418,  0.00708, 0.00996, 0.01281, 0.01563, 0.01843, &
        0.02121, 0.02396, 0.02668,  0.02938,  0.03206, 0.03471, 0.03733,  0.03994, 0.04251, 0.04506, 0.04759, 0.05009, &
        0.05256, 0.05501,  0.05744,  0.05984,  0.06221, 0.06457, 0.06689,  0.06919, 0.07147, 0.07372, 0.07594, 0.07814, &
        0.08032, 0.08247,  0.08460,  0.08670,  0.08877, 0.09082, 0.09285,  0.09485, 0.09683, 0.09878, 0.10070, 0.10261, &
        0.10448, 0.10633, 0.10816,  0.10996,  0.11174, 0.11349, 0.11521,  0.11691, 0.11859, 0.12024, 0.12187, 0.12347, &
        0.12505, (I*0.0,I=1,24)/)
    grids%fc = temp(startAge-20+1:85)

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

    params%nu =   0.372042538693115 !0.5!
    params%beta = 0.976948640333028 !0.99!
    params%gamma = 2.46102724400491 !0.4!
    params%db(1) = 0.764846423121866
    params%db(2) =  -5.472985265008665E-005
    params%thetab = 0.0375
    params%k=650000

    call getIncomeGrid(params, grids%Ygrid, grids%incTransitionMrx, minInc, grids%maxInc, grids%AIMEgrid, grids%benefit,grids%fc)

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

    action =2
    if (action .EQ. 1) then
        params%nu = 0.379444754253488     
        params%beta = 0.999958017173274     
        params%gamma = 2.50999213311961     
        params%db(1) = 0.780063950700460     
        params%db(2) =  -5.581876523249634E-005
        call getassetgrid( params, grids%maxInc, grids%Agrid)
        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, .TRUE. )
        !simulate
        if (rank == 0) then
            !call cpu_time(starsim)
            call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp, AIME )
            !call cpu_time(finish)
            !print '("Time = ",f11.3," seconds.")',finish-starsim
            call writetofile(params, ypath, cpath, apath, vpath, lpath, yemp,AIME)
        end if
    else

        if (params%system == 1 ) then
            open (unit = 1001,file='..\\..\\moments\\moments.txt', action='read', IOSTAT = ios)
        else
            open (unit = 1001,file='../moments/moments.txt', action='read', IOSTAT = ios)
        end if
        read (1001, *) moments(1,:)
        if (params%system == 1 ) then
            open (unit = 1002,file='..\\..\\moments\\assetmom.txt', action='read', IOSTAT = ios)
        else
            open (unit = 1002,file='../moments/assetmom.txt', action='read', IOSTAT = ios)
        end if
        read (1002, *) moments(2,:)

        if (rank==0) then
            print '("Setting up initial guess for hilling climbing algorithm")'
        end if

        call initialGuess(rank,params,grids,moments,p,y)

        call amoeba(p,y,0.001_rk,gmm_criteria,iter)

        if (rank==0) then
            print '("P = ",f6.3)',P(1,:)
            print '("Y = ",f16.3)',Y
            if (params%system == 1 ) then !ifort
                inquire (iolength=requiredl)  P(1,:)
                open (unit=201, file='..\\out\params.txt', status='unknown',recl=requiredl, action='write')
                write (201, * ) P(1,:)
            else !Gfort
                inquire (iolength=requiredl)  P
                open (unit=201, file='./out/params.txt', status='unknown',recl=requiredl, action='write')
                write (201, * ) P(1,1)
                write (201, * ) P(1,2)
                write (201, * ) P(1,3)
                write (201, * ) P(1,4)
                write (201, * ) P(1,5)
            end if
            close (unit=201)

            print '("Generating files")'
        end if

        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, .FALSE. )

        if (rank==0) then
            call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp ,AIME)
            call writetofile(params, ypath, cpath, apath, vpath, lpath, yemp,AIME)
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
    gmm_criteria = gmm(params,grids,moments) !*-1.0

    end function

    end program Console1

