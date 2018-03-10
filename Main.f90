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

    ! Variables
    type (structparamstype) :: params
    type (gridsType) :: grids
    real (kind=rk) :: minInc(Tperiods), maxInc(Tperiods), maxA(Tperiods+1)
    real (kind=rk) :: startA, span, loggrid(numPointsA)
    integer :: ixT, i
    real (kind=rk) :: V(Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: policyA1(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: policyC(Tperiods, numPointsA, numPointsY, numAIME)
    integer :: policyL(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: EV(Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: EdU(Tperiods,   numPointsA, numPointsY, numAIME)

    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: yemp(Tperiods, numSims)
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    real (kind=rk) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    !real (kind=rk) :: AIME(Tperiods + 1,numSims)

    !mwPointer matOpen, mxCreateDoubleMatrix, mxCreateString
    !    character (10) :: tempformat_numpeepcols
    !    character (20) :: format_numpeepcols
    !    character (20) :: format_numpeepcols_int
    real (kind=rk) :: start, finish, moments(2,24), error, nu
    integer :: action, ios, inter, requiredl
    INTEGER(kind=4) :: iter
    REAL(kind=8) :: y(dimEstimation+1)
    REAL(kind=8) :: p(dimEstimation+1,dimEstimation)
    
    CHARACTER(len=255) :: cwd
    CALL getcwd(cwd)
    WRITE(*,*) TRIM(cwd)

    params%system=0


    call cpu_time(start)

    params%r=0.016
    startA =0
    params%beta = 0.992!0.988164484822505 just labour!             ! Discount factor
    params%gamma = 3!             ! Coefficient of relative risk aversion
    params%nu =  0.256!0.294310003024158 just labour
    !0.287 elsa momement leisure!0.322 cormac's moemts!cormac published 0.4637 !lesuire
    
    params%mu = 0 !;                      ! mean of initial log income
    params%sigma = 0.0922 !mean([0.007 0.005 0.010 0.012])^0.5;! variance of innovations to log income
    params%rho = 0.96 !                ! persistency of log income
    params%hrsWrk = 0.3159 !1840/5824;
    params%delta(3) = 8.2270 !mean([9.083298    8.5256042   7.906901    7.392151]);
    params%delta(2)=0.0767 !mean([0.029333  0.061582208 0.094386    0.121426]);
    params%delta(1)=-7.2503e-04 !mean([-0.00023 -0.0005901  -0.00091    -0.00117]);
    params%spouseInc = 6.2318e+03!mean([5121, 5282, 6840, 7684]);
    params%tol = 1e-10 !                ! max allowed error
    params%minCons = 1e-5 !              ! min allowed consumption

    call getIncomeGrid(params, grids%Ygrid, grids%incTransitionMrx, minInc, maxInc, grids%AIMEgrid, grids%benefit)

    !Set maximum assets
    maxA(1) = startA;
    do ixt = 2, Tperiods+1
        maxA(ixt) = (maxA(ixt - 1) + maxInc(ixt-1) ) * (1+params%r)
    end do

    !Create asset grid
    do ixt = 1, Tperiods+1
        span =  (log(1.0+log(1.0+log(1+maxA(ixt)))) - log(1.0+log(1.0+log(1.0))) )/ (numPointsA-1)
        loggrid = log(1.0+log(1.0+log(1.0))) + span*(/(i,i=0,numPointsA-1)/)
        grids%Agrid(ixt, :) = (/(exp(exp(exp(loggrid(i))-1.0)-1.0)-1.0,i=1,numPointsA)/) !exp(exp(exp(loggrid)-1)-1)-1
    end do

    action = 0
    if (action .EQ. 1) then

        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, EdU, .TRUE. )
        !simulate
        call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp )
        !
        call writetofile(params, ypath, cpath, apath, vpath, lpath, yemp)

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

        !error = golden_generic(0.0_rk, 1.0_rk, nu, gmm_criteria,0.001_rk,.TRUE.)
        !print '("Nu = ",f6.3)',nu
        
        print '("Setting up initial guess for hilling climbing algorithm")'
        !Setting up initial guess for hilling climbing algorithm
        p(1,1) = params%nu
        p(1,2) = params%beta
        p(1,3) = params%gamma
        y(1) = gmm_criteria(p(1,:))
        p(2,1) = 0.4637
        p(2,2) = 0.970 
        P(2,3) = 1
        y(2) = gmm_criteria(p(2,:))
        p(3,1) = 0.322
        p(3,2) = 0.9843
        P(3,3) = 2
        y(3) = gmm_criteria(p(3,:))
        p(4,1) = 0.55
        p(4,2) = 0.96
        P(4,3) = 0.5
        y(4) = gmm_criteria(p(3,:))
        call amoeba(p,y,0.002_8,gmm_criteria,iter)
        print '("P = ",f6.3)',P
        print '("Y = ",f16.3)',Y
        !params%nu =
        if (params%system == 1 ) then !ifort
            inquire (iolength=requiredl)  P(1,:)
            open (unit=201, file='..\\out\params.txt', status='unknown',recl=requiredl, action='write')
            write (201, * ) P(1,:)
        else !Gfort
            inquire (iolength=requiredl)  P(1,:)
            open (unit=201, file='./out/params.txt', status='unknown',recl=requiredl, action='write')
            write (201, * ) P(1,:)
        end if
        print '("Generating files")'
        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, EdU, .FALSE. )
        call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp )
        call writetofile(params, ypath, cpath, apath, vpath, lpath, yemp)

    end if
    call cpu_time(finish)
    print '("Time = ",f7.3," seconds.")',finish-start


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
    if (minval(control) < 0)  then
        gmm_criteria = huge(gmm_criteria)
        return
    end if    
    params%nu = control(1)
    params%beta = control(2)
    params%gamma = control(3)
    gmm_criteria = gmm(params,grids,moments) !*-1.0

    end function

    end program Console1

