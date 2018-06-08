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
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    !mwPointer matOpen, mxCreateDoubleMatrix, mxCreateString
    !    character (10) :: tempformat_numpeepcols
    !    character (20) :: format_numpeepcols
    !    character (20) :: format_numpeepcols_int
    real (kind=rk) :: start, finish, moments(2,24), error, test, starsim
    integer :: action, ios, inter, requiredl
    INTEGER(kind=4) :: iter
    REAL(kind=rk) :: y(dimEstimation+1)
    REAL(kind=rk) :: p(dimEstimation+1,dimEstimation)

    !CHARACTER(len=255) :: cwd
    !CALL getcwd(cwd)
    !WRITE(*,*) TRIM(cwd)

    params%system=1


    call cpu_time(start)

    !grids%fc = (/(I*0.0,I=1,85)/)
    grids%fc = (/-0.01374, -0.01069, -0.00767, -0.00467, -0.00169, 0.00126, 0.00418,  0.00708, 0.00996, 0.01281, 0.01563, 0.01843, &
        0.02121, 0.02396, 0.02668,  0.02938,  0.03206, 0.03471, 0.03733,  0.03994, 0.04251, 0.04506, 0.04759, 0.05009, &
        0.05256, 0.05501,  0.05744,  0.05984,  0.06221, 0.06457, 0.06689,  0.06919, 0.07147, 0.07372, 0.07594, 0.07814, &
        0.08032, 0.08247,  0.08460,  0.08670,  0.08877, 0.09082, 0.09285,  0.09485, 0.09683, 0.09878, 0.10070, 0.10261, &
        0.10448, 0.10633, 0.10816,  0.10996,  0.11174, 0.11349, 0.11521,  0.11691, 0.11859, 0.12024, 0.12187, 0.12347, &
        0.12505, (I*0.0,I=1,24)/)
    params%r= 0.02  !0.016
    params%startA =0 !415.87*20!
    params%beta =  0.98 !0.999999849491823 !0.98716208644227588 !0.99485273766567861 !0.992!0.988164484822505 just labour!             ! Discount factor
    params%gamma = 3 !1.610884853321348E-003 !2.16705246913580    ! 6.3251268401011834E-003 !3!             ! Coefficient of relative risk aversion
    params%nu =  0.466352772226276  !0.308994215349724!0.28685133136511892 !0.256!0.294310003024158 just labour
    !0.287 elsa momement leisure!0.322 cormac's moemts!cormac published 0.4637 !lesuire
    params%db(1) = 0.893279806195867     !0.8334
    params%db(2) = -(6.9447e-06)!-8.223557117746168E-002 !


    params%mu = 0 !;                      ! mean of initial log income
    params%sigma = 0.0922 !0.0922 !mean([0.007 0.005 0.010 0.012])^0.5;! variance of innovations to log income
    params%rho = 0.96 !0.96 !                ! persistency of log income
    params%hrsWrk = 0.3159 !1840/5824;
    params%delta(3) = 9.083298 !8.2270 !mean([9.083298    8.5256042   7.906901    7.392151]);
    params%delta(2)=0.029333 !0.0767 !mean([0.029333  0.061582208 0.094386    0.121426]);
    params%delta(1)=-0.00033!-0.00023 !-7.2503e-04 !mean([-0.00023 -0.0005901  -0.00091    -0.00117]);
    params%spouseInc = 6235.8998884204320 !6.2318e+03!mean([5121, 5282, 6840, 7684]);
    params%tol = 1e-10 !                ! max allowed error
    params%minCons = 1e-5 !              ! min allowed consumption

    call getIncomeGrid(params, grids%Ygrid, grids%incTransitionMrx, minInc, grids%maxInc, grids%AIMEgrid, grids%benefit,grids%fc)

    !        params%nu =    0.38022456150504280 !0.339785209413524 !0.466352772226276 !
    !        params%beta = 0.97235322545400193 !0.958731419050955 ! 0.999999849491823 !
    !        params%gamma =  2.0920418173800615  !1.89794757374393 !1.610884853321348E-003 !
    !        !0.287 elsa momement leisure!0.322 cormac's moemts!cormac published 0.4637 !lesuire
    !        params%db(1) = 0.80643840150485868  !0.893686680412417 !0.893279806195867 !
    !        params%db(2) = -4.7393420983952571E-005 !-4.050600090249978E-005 !-8.223557117746168E-002!
    params%nu =      0.38022456150504280 !0.339785209413524 !0.466352772226276 !
    params%beta =  0.97235322545400193 !0.958731419050955 ! 0.999999849491823 !
    params%gamma =  2.0920418173800615  !1.89794757374393 !1.610884853321348E-003 !
    !0.287 elsa momement leisure!0.322 cormac's moemts!cormac published 0.4637 !lesuire
    params%db(1) = 0.91387622628345655 !0.80643840150485868  !0.893686680412417 !0.893279806195867 !
    params%db(2) = -4.7393420983952571E-005 !-4.050600090249978E-005 !-8.223557117746168E-002!


    action =2
    if (action .EQ. 1) then
        !        params%nu =     0.37490145517514262 !0.33102935844275810! 0.38022456150504280 !0.339785209413524 !0.466352772226276 !
        !        params%beta = 0.97401163475971553 !0.99917272385048772! 0.97235322545400193 !0.958731419050955 ! 0.999999849491823 !
        !        params%gamma =   1.9925202627061784 !4.2876825152294789 !2.0920418173800615  !1.89794757374393 !1.610884853321348E-003 !
        !        !0.287 elsa momement leisure!0.322 cormac's moemts!cormac published 0.4637 !lesuire
        !        params%db(1) = 0.92384410362795921! 0.91387622628345655 !0.80643840150485868  !0.893686680412417 !0.893279806195867 !
        !        params%db(2) = -4.6263715941694290E-005!-1.5160228027272263E-005!-4.7393420983952571E-005 !-4.050600090249978E-005 !-8.223557117746168E-002!
        params%nu =     0.37490145517514262!0.38022456150504280 !
        params%beta = 0.97401163475971553!0.97235322545400193 !
        params%gamma = 1.9925202627061784!2.0920418173800615 !
        params%db(1) = 0.92384410362795921!0.80643840150485868  !
        params%db(2) = -4.6263715941694290E-005!-4.7393420983952571E-005 !

        call getassetgrid( params, grids%maxInc, grids%Agrid)
        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, EdU, .TRUE. )
        !simulate
        call cpu_time(starsim)
        call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp, AIME )
        call cpu_time(finish)
        print '("Time = ",f11.3," seconds.")',finish-starsim

        call writetofile(params, ypath, cpath, apath, vpath, lpath, yemp,AIME)

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

        !!$omp parallel default(shared)   
        !!$omp sections
        !!$omp section
        p(1,1) = params%nu
        p(1,2) = params%beta
        p(1,3) = params%gamma
        !p(1,4) = params%spouseinc
        !p(1,5) = params%StartA
        p(1,4) = params%db(1)
        p(1,5) = params%db(2)
        !p(1,6) = params%StartA
        y(1) = gmm_criteria(p(1,:))

        !!$omp section
        p(2,1) = 0.4637
        p(2,2) = 0.970
        P(2,3) = 1
        !p(2,4) = params%spouseinc*1.1
        !p(2,5) = params%StartA*1.2
        p(2,4) = params%db(1)*0.9
        p(2,5) = params%db(2)*1.2
        !p(2,6) = params%StartA*1.1
        y(2) = gmm_criteria(p(2,:))

        !!$omp section
        p(3,1) = 0.322
        p(3,2) = 0.9843
        P(3,3) = 2
        !p(3,4) = params%spouseinc*0.9
        !p(3,5) = params%StartA*0.7
        p(3,4) = params%db(1)*1.1
        p(3,5) = params%db(2)*0.7
        !p(3,6) = params%StartA*1.67
        y(3) = gmm_criteria(p(3,:))

        !!$omp section
        p(4,1) = 0.55
        p(4,2) = 0.96
        P(4,3) = 0.5
        !p(4,4) = params%spouseinc*0.95
        !p(4,5) = params%StartA*0.95
        p(4,4) = params%db(1)*1.3
        p(4,5) = params%db(2)*0.95
        !p(4,6) = params%StartA*0.89
        y(4) = gmm_criteria(p(4,:))

        !!$omp section
        p(5,1) = 0.15
        p(5,2) = 0.9999
        P(5,3) = 4
        !p(5,4) = params%spouseinc*1.05
        !p(5,5) = params%StartA*1.15
        p(5,4) = params%db(1)*0.85
        p(5,5) = params%db(2)*1.15
        !p(5,6) = params%StartA*1.3
        y(5) = gmm_criteria(p(5,:))

        !!$omp section
        p(6,1) = 0.27
        p(6,2) = 0.986
        P(6,3) = 0.9
        !p(6,4) = params%spouseinc*1.01
        !p(6,5) = params%StartA*0.87
        p(6,4) = params%db(1)*0.99
        p(6,5) = params%db(2)*0.87
        !p(6,6) = params%StartA*0.4
        y(6) = gmm_criteria(p(6,:))

        !!$omp section
        !        p(7,1) = 0.27
        !        p(7,2) = 0.986
        !        P(7,3) = 0.9
        !        !p(7,4) = params%spouseinc*1.01
        !        !p(7,5) = params%StartA*0.87
        !        p(7,4) = params%db(1)*2
        !        p(7,5) = params%db(2)*0.4
        !        p(7,6) = params%StartA*1.8
        !        y(7) = gmm_criteria(p(7,:))
        !!$omp end sections
        !!$omp end parallel
        call amoeba(p,y,0.02_rk,gmm_criteria,iter) !0.002_8

        print '("P = ",f6.3)',P(1,:)
        print '("Y = ",f16.3)',Y
        !params%nu =
        if (params%system == 1 ) then !ifort
            inquire (iolength=requiredl)  P(1,:)
            open (unit=201, file='..\\out\params.txt', status='unknown',recl=requiredl, action='write')
            write (201, * ) P(1,:)
        else !Gfort
            inquire (iolength=requiredl)  P
            open (unit=201, file='./out/params.txt', status='unknown',recl=requiredl, action='write')
            !write (201, * ) P(1,:)
            write (201, * ) P(1,1)
            write (201, * ) P(1,2)
            write (201, * ) P(1,3)
            write (201, * ) P(1,4)
            write (201, * ) P(1,5)
            !write (201, * ) P(1,6)
        end if
        close (unit=201)

        print '("Generating files")'
        call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, EdU, .FALSE. )
        call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp ,AIME)
        call writetofile(params, ypath, cpath, apath, vpath, lpath, yemp,AIME)

    end if
    call cpu_time(finish)
    print '("Time = ",f11.3," seconds.")',finish-start


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
    params%StartA= control(6)
    !params%spouseinc = control(3)
    gmm_criteria = gmm(params,grids,moments) !*-1.0

    end function

    end program Console1

