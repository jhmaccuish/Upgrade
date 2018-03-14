    module routines
    use Header
    use routines_generic

    implicit none

    contains

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Get Income grid
    subroutine getIncomeGrid(params, Ygrid, incTransitionMrx, minInc, maxInc, AIMEgrid, benefit, fc)
    implicit none

    !inputs
    type (structparamstype), intent(inout) :: params
    real (kind=rk):: fc(:)

    !outputs
    real (kind=rk) :: Ygrid(:,:), incTransitionMrx(:,:), minInc(:), maxInc(:), AIMEgrid(:,:), benefit(:)

    !local
    real (kind=rk) :: sig_inc, ly(numPointsY), upper(Tperiods+1), a !Q(numPointsY,numPointsY)
    integer :: t, i

    sig_inc = params%sigma/((1-params%rho**2)**0.5)

    !Why 3
    call tauchen(numPointsY,params%mu,params%rho,params%sigma,3,ly,incTransitionMrx)

    !Ygrid = exp(repmat(ly', T, 1)+repmat(polyval(delta,1:T)'-fc(1:T),1,numPointsY));
    !minInc = exp(repmat((-normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((-normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));
    !maxInc = exp(repmat((normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));

    params%pension = 107.45*52
    upper(1) = 0
    do t=1 , Tperiods
        Ygrid(t,:)= exp(ly+params%delta(1)*t**2+params%delta(2)*t+params%delta(3)-fc(t))
        minInc(t) = exp((-normBnd * sig_inc)+params%delta(1)*t**2+params%delta(2)*t+params%delta(3)-fc(t))
        maxInc(t) = exp((normBnd * sig_inc)+params%delta(1)*t**2+params%delta(2)*t+params%delta(3)-fc(t))
        upper(t+1) = upper(t) + maxInc(t)
        if (t <=spouseretire) then
            a = (upper(t+1)/t)/(numAIME-1)
            AIMEgrid(t,:) = a*(/(i,i=0,numAIME-1)/) !linspace(0,upper(t+1)/t,numAIME);
        else
            AIMEgrid(t,:) = AIMEgrid(t-1,:)
        end if
        if (t <=5) then
            benefit(t) = 57.90*52
        else
            benefit(t) = 73.10*52
        end if
    end do
    AIMEgrid(Tperiods+1,:) = AIMEgrid(Tperiods,:)

    !benefit = [57.90*52*ones(5,1);73.10*52*ones(length(minInc)-5,1)]

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Inpliments tauchen to discretise AR1
    subroutine tauchen(N,mu,rho,sigma,m,Z,Zprob)

    implicit none
    !Function TAUCHEN
    !
    !Purpose:    Finds a Markov chain whose sample paths
    !            approximate those of the AR(1) process
    !                z(t+1) = (1-rho)*mu + rho * z(t) + eps(t+1)
    !            where eps are normal with stddev sigma
    !
    !Format:     {Z, Zprob} = Tauchen(N,mu,rho,sigma,m)
    !
    !Input:      N       scalar, number of nodes for Z
    !            mu      scalar, unconditional mean of process
    !            rho     scalar
    !            sigma   scalar, std. dev. of epsilons
    !            m       max +- std. devs.
    !
    !Output:     Z       N*1 vector, nodes for Z
    !            Zprob   N*N matrix, transition probabilities
    !
    !    Martin Floden
    !    Fall 1996
    !
    !    This procedure is an implementation of George Tauchen's algorithm
    !    described in Ec. Letters 20 (1986) 177-181.
    !

    !inputs
    integer, intent(in) :: N, m
    real (kind=rk), intent(in):: mu, rho, sigma

    !outputs
    real (kind=rk), intent(out) :: Z(N), Zprob(N,N)

    !local
    real (kind=rk) :: a, zstep
    integer :: i, j, k

    a     = (1-rho)*mu;

    Z(N)  = m * sqrt(sigma**2 / (1 - rho**2))
    Z(1)  = -Z(N)
    zstep = (Z(N) - Z(1)) / (N - 1)

    do i = 2, (N-1)
        Z(i) = Z(1) + zstep * (i - 1)
    end do

    Z = Z + a / (1-rho);

    do j = 1, N
        do k = 1, N
            if (k == 1) then
                Zprob(j,k) = cdf_normal((Z(1) - a - rho * Z(j) + zstep / 2) / sigma)
            elseif (k == N) then
                Zprob(j,k) = 1 - cdf_normal((Z(N) - a - rho * Z(j) - zstep / 2) / sigma)
            else
                Zprob(j,k) = cdf_normal((Z(k) - a - rho * Z(j) + zstep / 2) / sigma) - &
                    cdf_normal((Z(k) - a - rho * Z(j) - zstep / 2) / sigma);
            end if
        end do
    end do

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Normal CDF
    function cdf_normal(x)
    REAL(KIND=rk), INTENT(in)  :: X
    REAL(KIND=rk) :: cdf_normal

    cdf_normal =  0.5d0*( 1+erf(x/SQRT(2.d0)) )

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Solve value function by backward induction
    subroutine solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, EdU, show )
    implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    logical :: show

    !outputs
    real (kind=rk), intent(out) :: V(Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(out) :: policyA1(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(out) :: policyC(Tperiods, numPointsA, numPointsY, numAIME)
    integer, intent(out) :: policyL(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(out) :: EV(Tperiods+1, numPointsA, numPointsY, numAIME);
    real (kind=rk), intent(out) :: EdU(Tperiods,   numPointsA, numPointsY, numAIME);

    !local
    integer :: ixt, ixAIME, ixA, ixY, ixL
    real (kind=rk) :: negV, A, Y, lbA1, ubA1, AIME, EV1(numPointsA ,numAIME) ! Agrid1(numPointsA)
    real (kind=rk) :: AIME1grid(numAIME), policyA1temp, negVtemp, realisedV(numPointsY)

    !Test
    real (kind=rk) :: testC

    !Set the terminal value function and expected value function to 0
    EV(Tperiods + 1, :,:,:)  = 0;          ! continuation value at T-1
    V(Tperiods + 1,:,:,:) = 0;

    do ixt=Tperiods,1, -1                               ! Loop from time T-1 to 1
        !Agrid1 = grids%Agrid(ixt + 1, :);               ! The grid on assets tomorrow
        AIME1grid = grids%AIMEgrid(ixt + 1, :);
        do ixAIME = 1, numAIME
            do ixA = 1, numPointsA                   ! points on asset grid
                if (ixt < stopwrok) then
                    !Although doesn't recieve income still need to loop around
                    !hypothetical income because participation doesn't effect
                    !earning potential
                    ! STEP 1. solve problem at grid points in assets, income + labour choices
                    ! ---------------------------------------------------------
                    do ixY = 1, numPointsY               ! points on income grid
                        !negV = -huge(negv) !-log(0.0) !inf
                        lbA1 = grids%Agrid(ixt + 1, 1);          ! lower bound: assets tomorrow
                        EV1  = EV(ixt + 1,:, ixY,:);  ! relevant section of EV matrix (in assets tomorrow)
                        call solvePeriod(params, grids, grids%Ygrid(ixt, ixY),grids%Agrid(ixt, ixA), grids%AIMEgrid(ixt,ixAIME), &
                                ixt, lbA1, EV1, grids%benefit(ixt), policyA1(ixt,ixA,ixY,ixAIME), &
                                policyC(ixt, ixA, ixY,ixAIME), policyL(ixt,ixA,ixY,ixAIME), V(ixt, ixA, ixY,ixAIME))
!                        do ixL = 0,(numPointsL-1),1           ! points of income choice
!                            ! Value of income and information for optimisation
!                            A    = grids%Agrid(ixt, ixA);            ! assets today
!                            Y    = ixL*grids%Ygrid(ixt, ixY)+ (1-ixL)*grids%benefit(ixt);
!                            AIME = grids%AIMEgrid(ixt,ixAIME);
!                            !end if
!                            !AIME only depends on earned income so add spousal
!                            !and pensions after calculating it
!                            call gross2net(params,Y,ixt,ixl,AIME)
!                            ! Next peridos AIME
!                            if (ixL==1 .AND.  ixt < spouseretire) then
!                                    AIME =  grids%Ygrid(ixt, ixY)/ixt + AIME * (ixt-1)/ixt
!                            else if ((ixL==0 .AND.  ixt < spouseretire)) then
!                                    AIME = AIME * (ixt-1)/ixt
!                            end if
!
!
!                            ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
!
!
!                            ! Compute solution
!                            if (ubA1 - lbA1 < params%minCons) then         ! if liquidity constrained
!                                negVtemp = objectivefunc(params, grids,lbA1, A, Y,ixL,ixt, AIME,EV1);
!                                policyA1temp = lbA1;
!                            else                               ! if interior solution
!                                ![policyA1temp, negVtemp] = ...
!                                !    fminbnd(@(A1) objectivefunc(A1, A, Y,ixL,ixt, AIME), lbA1, ubA1, optimset('TolX',tol));
!                                negVtemp = golden_generic(lbA1, ubA1, policyA1temp, func,params%tol,.FALSE.)
!                            end if! if (ubA1 - lbA1 < minCons)
!                            if (negVtemp > negV) then
!                                negV = negVtemp;
!                                policyA1(ixt,ixA,ixY,ixAIME)=policyA1temp;
!                                policyL(ixt,ixA,ixY,ixAIME)=ixL;
!                                ! Store solution and its value
!                                policyC(ixt, ixA, ixY,ixAIME) = A + Y - policyA1(ixt, ixA, ixY,ixAIME)/(1+params%r);
!                            end if
!                        end do
!                        !end
!                        testC = policyC(ixt, ixA, ixY,ixAIME)
!                        V(ixt, ixA, ixY,ixAIME)       = negV;
                        !dU(ixt, ixA, ixY)      = getmargutility(policyC(ixt, ixA, ixY),policyL(ixt,ixA,ixY));
                    end do
                else
                    negV = -huge(negv) !-log(0.0) !inf
                    ixL = 0
                    ! Value of income and information for optimisation
                    A    = grids%Agrid(ixt, ixA)            ! assets today
                    Y    = grids%benefit(ixt)
                    AIME = grids%AIMEgrid(ixt,ixAIME);
                    call gross2net(params,Y,ixt,ixl,AIME)
                    if (y < 0) then
                        print *, 'Error Y< 0'
                    end if

                    lbA1 = grids%Agrid(ixt + 1, 1);          ! lower bound: assets tomorrow
                    ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
                    EV1  = EV(ixt + 1,:, 1,:);  ! relevant section of EV matrix (in assets tomorrow)
                    ! Compute solution
                    if (ubA1 - lbA1 < params%minCons) then         ! if liquidity constrained
                        negVtemp = objectivefunc(params, grids,lbA1, A, Y,ixL,ixt, AIME,EV1);
                        policyA1temp = lbA1;
                    else                               ! if interior solution
                        ![policyA1temp, negVtemp] = ...
                        !    fminbnd(@(A1) objectivefunc(A1, A, Y,ixL,ixt, AIME), lbA1, ubA1, optimset('TolX',tol));
                        negVtemp = golden_generic(lbA1, ubA1, policyA1temp, func,params%tol,.FALSE.)
                    end if! if (ubA1 - lbA1 < minCons)
                    if (negVtemp > negV) then
                        negV = negVtemp;
                        policyA1(ixt,ixA,(/(ixY, ixY=1, numPointsY)/),ixAIME)=policyA1temp;
                        policyL(ixt,ixA,(/(ixY, ixY=1, numPointsY)/),ixAIME)=ixL;
                        ! Store solution and its value
                        policyC(ixt, ixA, (/(ixY, ixY=1, numPointsY)/),ixAIME) = A + Y - policyA1(ixt, ixA, 1,ixAIME)/(1+params%r);
                    end if
                    !end
                    testC = policyC(ixt, ixA, 1,ixAIME)
                    V(ixt, ixA, (/(ixY, ixY=1, numPointsY)/),ixAIME)       = negV;
                    !dU(ixt, ixA, ixY)      = getmargutility(policyC(ixt, ixA, ixY),policyL(ixt,ixA,ixY));
                end if

                ! STEP 2. integrate out income today conditional on income
                ! yesterday to get EV and EdU
                ! --------------------------------------------------------
                realisedV(:) = V(ixt, ixA, :, ixAIME);
                !realiseddU(:,:) = dU(ixt, ixA, :);
                do ixY = 1,numPointsY,1
                    EV(ixt, ixA, ixY, ixAIME)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV);
                    !EdU(ixt, ixA, ixY) = incTransitionMrx(ixY,:)*realiseddU;
                end do !ixY

            end do!ixA
        end do!ixAIME

        if (show) WRITE(*,*)  'Passed period ', ixt, ' of ',  Tperiods
    end do!ixt
    contains
    function func(x)
    real (kind = rk), intent(in) :: x
    real (kind = rk) :: func
    func = objectivefunc(params, grids,x, A, Y,ixL,ixt, AIME,EV1)
    end function

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Objective function
    function objectivefunc(params, grids, A1, A0, Y,L,ixP, AIME, EV1)
    implicit none
    ! This function returns the following quantity:
    ! - (u(c) +  b V( A1))
    ! where c is calculated from today's assets and tomorrow's assets

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    real (kind=rk), intent(in) :: A1, A0, Y, AIME, EV1(:,:)
    integer, intent(in) :: ixP, L
    !ouptut
    real (kind=rk) :: objectivefunc
    !local
    real (kind=rk) :: cons, VA1, mortal(Tperiods+1)
    ! Declare global we need this file have access to
    !mortal = (/1,2,&
    !        3,4/)

    mortal = (/0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011071, 0.011907, 0.012807, &
        0.013676, 0.01475, 0.015818, 0.016846, 0.018079, 0.019343, 0.020659, 0.0218, 0.023505, 0.025202, 0.02696, &
        0.028831, 0.031017, 0.033496, 0.036024, 0.038912, 0.042054, 0.045689, 0.049653, 0.054036, 0.05886, 0.064093, &
        0.069636, 0.07533, 0.081069, 0.086912, 0.093067, 0.099807, 0.107483, 0.116125, 0.125196, 0.134361, 0.143881, &
        0.1542, 0.165675, 0.17842, 0.192363, 0.2117, 0.1672,0.1565, 0.1485,0.1459, 1.0/);

    !Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    !total value (u(c) + b * VA1
    cons = A0  + Y - (A1)/(1+params%r);
    !VA1 = interp1(Agrid1,EV1,A1, interpMethod, 'extrap');
    ![X,Y] = meshgrid(Agrid1,AIME1grid);
    !VA1 = interp2(Agrid1, AIME1grid , EV1, A1, AIME);

    call linearinterp2_withextrap(grids%Agrid(ixP + 1, :), grids%AIMEgrid(ixP + 1, :), &
        numPointsA, numAIME, A1, AIME, VA1, 1, 1, EV1)
    !interp2(EV1, A1/Agrid1(20), AIME/AIME1grid(10));
    objectivefunc = utility(params,cons,L) + params%beta * (1- mortal(ixP))* VA1;

    !! ------------------------------------------------------------------------
    !The optimisation routine that we will use searches for the minimum of the
    !function. We want the maximum. So we multiply out function here by -1 so
    !that the optimiser will fill the minimum of the negative of our function,
    !i.e. the maximum of our functino

    !objectivefunc = - objectivefunc;

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Utility Function
    function utility(params,cons,L)
    implicit none
    !This function takes consumption as an argument and returns utility. The
    !utility functin is CRRA except that we add a very small number (eps) to
    !consumption so that the computer can deal wi
    !inputs
    type (structparamstype), intent(in) :: params
    integer, intent(in) :: L
    real (kind=rk), intent(in) :: cons
    !real (kind=rk), intent(in) :: A1
    !outpus
    real (kind=rk) :: utility, les

    if (cons<=0) then
        print *, 'Error in utility! Consumption is LE 0'
        stop
    end if
    !10/112 comuniting time -10/112
    les=(L)*(1-params%hrsWrk -10/112)+(1-L);
    if (params%gamma == 1) then
        utility = log(cons**params%nu*les**(1-params%nu));
    else
        utility= ((cons**params%nu*les**(1-params%nu))**(1-params%gamma)  )/(1-params%gamma);
    end if

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Defined Benefite Pension function
    function dbPension(params,AIME)
    implicit none
    !calibrated to give 25,000 max from 60,000 AIME upwards
    !db2 = level/(point^2-2point), db1 = 6.9447e-06*points
    !point

    !inputs
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: AIME
    !outputs
    real (kind=rk) :: dbPension
    !local
    real (kind=rk):: bound
    bound = -params%db(1) /(2.0*params%db(2))
    if (AIME < bound) then
        dbPension = params%db(2)*AIME**2 + params%db(1)*AIME;
    else
        dbPension = params%db(2)*bound**2 + params%db(1)*bound !25000
    endif
    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Simulation Subroutine
    subroutine simWithUncer(params, grids, policyA1,policyL,EV,  y, c, a, v, l, yemp, AIME )
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    real (kind=rk), intent(in) :: policyA1(Tperiods, numPointsA, numPointsY, numAIME)
    !real (kind=rk), intent(in) :: policyC(Tperiods, numPointsA, numPointsY, numAIME)
    integer, intent(in) :: policyL(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk), intent(in) :: EV(Tperiods+1, numPointsA, numPointsY, numAIME);

    !outputs
    real (kind=rk), intent(out) :: y(Tperiods, numSims) !income
    real (kind=rk), intent(out) :: yemp(Tperiods, numSims)
    real (kind=rk), intent(out) :: c(Tperiods, numSims)  !consumption
    integer, intent(out) :: l(Tperiods, numSims) !labour supply
    real (kind=rk), intent(out) :: v(Tperiods, numSims)  !value
    real (kind=rk), intent(out) :: a(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(out) :: AIME(Tperiods + 1,numSims)

    !real (kind=rk) :: AIME(Tperiods + 1,numSims)

    !local
    real (kind=rk) :: startingA(numSims), startAIME, sig_inc, sig_initial
    real (kind=rk) :: e(Tperiods, numSims), temp(numSims, Tperiods)    ! the innovations to log income
    real (kind=rk) :: logy1(1,numSims)        ! draws for initial income
    real (kind=rk) :: ly(Tperiods, numSims)           ! log income
    !real (kind=rk) :: ypathIndex(Tperiods, numSims)   ! holds the index (location) in the vector
    integer :: s, t, seed1, seed2, ios, idxA(1), idxAIME(1), idxY(1)
    !real (kind=rk) :: temp2(numPointsY)
    !CHARACTER(len=255) :: cwd
    !CALL getcwd(cwd)

    integer :: seedIn

    !local
    real (kind=rk) ::  uniformRand(numSims), ltemp, lbA1, EV1(numPointsA ,numAIME)
    INTEGER :: n, i, uniformInt(numSims)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    seedIn = 16101988
    !Set seed
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    seed = seedIn * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)

    !!get uniform random number
    CALL RANDOM_NUMBER(uniformRand)

    do i=1, numsims
        uniformInt(i) = nint(uniformRand(i)*numPointsA)
        startingA(i) =grids%Agrid(1,uniformInt(i))
    end do
    startAIME = 0

    ! Obtain time series of incomes for our simulated individuals
    ! Draw random draws for starting income and for innovations
    seed1 = 1223424; ! For the innovations
    seed2 = 234636;  ! For initial income
    sig_inc = params%sigma/ ((1-params%rho**2)**0.5)
    sig_initial = 0.2950 !mean([0.073 0.053 0.110 0.112])^0.5;
    !e = getNormalDraws( 0.0_rk, params%sigma,  T, numSims, seed1);  ! normally distributed random draws for the innovation
    if (params%system == 1 ) then
        open (unit = 1001,file='..\\Errors.txt', action='read', IOSTAT = ios)
    else
        open (unit = 1001,file='./data/Errors.txt', action='read', IOSTAT = ios)
        !open (unit = 1001,file='..\\..\\..\\Errors.txt', action='read', IOSTAT = ios)
        !open (unit = 1001,file='./data/Errors.txt', action='read', IOSTAT = ios)
    end if
    read (1001, *) temp
    e = transpose(temp)
    !logy1 =  getNormalDraws( params%mu, sig_initial,  1, numSims, seed2); ! a random draw for the initial income
    if (params%system == 1 ) then
        open (unit = 1001,file='..\\IntProd.txt', action='read', IOSTAT = ios)
    else
        open (unit = 1001,file='./data/IntProd.txt', action='read', IOSTAT = ios)
        !open (unit = 1001,file='..\\..\\..\\IntProd.txt', action='read', IOSTAT = ios)
        !open (unit = 1001,file='./data/IntProd.txt', action='read', IOSTAT = ios)
    end if

    read (1001, *) logy1

    ! Get all the incomes, recursively
    do s = 1, numSims, 1                           ! loop through individuals
        ly(1, s) = truncate(logy1(1, s), -normBnd*sig_inc,normBnd*sig_inc )
        y(1, s) = exp(ly(1, s)+params%delta(1)+params%delta(2)+params%delta(3)-grids%fc(1))
        do t = 1,Tperiods-1,1                              ! loop through time periods for a particular individual
            ly(t+1, s) = (1 -params%rho) * params%mu + params%rho * ly(t, s) + e(t + 1, s)
            ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
            y(t+1, s) = exp( ly(t+1, s) + params%delta(1)*(t+1)**2+params%delta(2)*(t+1)+params%delta(3)-grids%fc(t+1) )
        end do ! t
    end do! s

    do s = 1,numSims,1
        a(1, s) = startingA(s)
        AIME(1,s)=startAIME
        do t = 1,Tperiods,1                              ! loop through time periods for a particular individual
            !Should improve this nearest neigbhour
            idxY=minloc(abs(y(t, s)-grids%Ygrid(t,:)))
            idxA=minloc(abs(a(t, s)-grids%Agrid(t,:)))
            idxAIME=minloc(abs(AIME(t, s)-grids%AIMEgrid(t,:)))
            l(t,s)=policyL(t,idxA(1),idxY(1),idxAIME(1))
            yemp(t,s) = y(t, s)

            call linearinterp3_withextrap(grids%Agrid(t,:), grids%Ygrid(t, :), grids%AIMEgrid(t,:), &
                numPointsA, numPointsY, numAIME,  a(t, s), y(t, s),AIME(t, s), ltemp,  real(policyL(t, :, :,:),rk))
            if (abs(l(t,s) - ltemp) > 0.01) then
                lbA1 = grids%Agrid(t + 1, 1);          ! lower bound: assets tomorrow
                ev1 = EV(t + 1,:, idxY(1),:)
                call solvePeriod(params, grids, y(t, s), a(t, s), AIME(t, s) ,t, lbA1, ev1, &
                grids%benefit(t),a(t+1, s),c(t, s),l(t,s),v(t  , s))
            else
                !a(t+1, s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tA1, a(t, s), (y(t, s)));
                !interp3(Agrid(t,:)', Ygrid(t, :)', AIMEgrid(t,:)',tA1, a(t, s), y(t, s),AIME(t, s));
                call linearinterp3_withextrap(grids%Agrid(t,:), grids%Ygrid(t, :), grids%AIMEgrid(t,:), &
                    numPointsA, numPointsY, numAIME,  a(t, s), y(t, s),AIME(t, s), a(t+1, s),  policyA1(t, :, :,:))
                !v(t  , s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tV , a(t, s), (y(t, s)));
                !interp3(Agrid(t,:)', Ygrid(t, :)', AIMEgrid(t,:)',tV, a(t, s), y(t, s),AIME(t, s));
                call linearinterp3_withextrap(grids%Agrid(t,:), grids%Ygrid(t, :), grids%AIMEgrid(t,:),&
                    numPointsA, numPointsY, numAIME,  a(t, s), y(t, s),AIME(t, s), v(t  , s),  EV(t, :, :,:))
                                ! Get consumption from today's assets, today's income and
                ! Check whether next period's asset is below the lowest
                ! permissable
                !if ( a(t+1, s) < Agrid(t+1, 1))
                !   [ a(t+1, s) ] = checkSimExtrap( Agrid(t+1, 1),y(t, s), t );
                !end
            end if
            if (l(t,s) .EQ. 0) then
                y(t, s)=grids%benefit(t)
                AIME(t+1, s) =   AIME(t, s) * (t-1)/t
            else
                if  (t < 65) then
                    AIME(t+1, s) =  y(t,s)/t + AIME(t, s) * (t-1)/t
                else
                    AIME(t+1, s) =   AIME(t, s)
                end if
            end if

            call gross2net(params,y(t, s),t,l(t,s), AIME(t, s))
            ! tomorrow's optimal assets
            c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+params%r))

        end do   !t
    end do! s
    !yemp = yemp*l


    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Get Draws from a normal distribution
    function getNormalDraws( mu, sigma,  dim1, dim2, seedIn)
    implicit none
    !Inputs
    real (kind=rk), intent(in):: mu, sigma
    integer, intent(in) :: dim1, dim2, seedIn

    !output
    real (kind=rk) :: getNormalDraws(dim1,dim2)

    !local
    real (kind=rk) ::  uniformRand(dim1,dim2,2) !,StdRandNormal(dim1,dim2)
    INTEGER :: i, n, indx1, indx2
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    real (kind=rk) ::  r, theta
    real (kind=rk), PARAMETER :: PI=3.141592653589793238462

    !Set seed
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    seed = seedIn * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)

    !!get uniform random number
    CALL RANDOM_NUMBER(uniformRand)

    do indx1 = 1,dim1
        do indx2 = 1, dim2
            r = (-2.0d0*log(uniformRand(indx1,indx2,1)))**0.5
            theta = 2.0d0*PI*uniformRand(indx1,indx2,2)
            getNormalDraws(indx1, indx2) = mu+sigma*r*sin(theta)
        end do
    end do

    !Convert to normal
    !getNormalDraws = mu  + sigma * StdRandNormal

    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!truncate obsrvtion
    function truncate(y, negtrunc, postrunc)
    implicit none
    !input
    real (kind=rk), intent(in) :: y, negtrunc, postrunc

    !output
    real (kind=rk) :: truncate

    ! Truncate input if its too big or too small
    if (y < negtrunc) then ! If y is less than the value negtrunc
        truncate = negtrunc;
    elseif (y > postrunc) then
        truncate = postrunc; ! If y is greater than the value postrunc
    else
        truncate = y;
    end if
    end function

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Convert Groos to net income
    subroutine gross2net(params,Y,ixt,ixl, AIME)
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    integer, intent(in) :: ixt,ixl
    real (kind=rk), intent(in) :: AIME

    !changing
    real (kind =rk), intent(inout) :: Y

    !Y = Y +params%spouseInc
    !Add own and husbands pension
    if (ixt >= Tretire) then
        Y = Y + params%pension
        if  (ixt >= spouseretire) then
            Y = Y + params%pension
            if (ixL==0) then
                Y =   Y + dbPension(params,AIME);
            end if
        else
            Y = Y +params%spouseInc
        end if
    else
        Y = Y +params%spouseInc
    end if

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Returns GMM Cretieria
    function gmm(params,grids,target)
    implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(inout) :: grids
    real (kind=rk), intent(in) :: target(:,:)
    !output
    real (kind=rk) :: gmm
    !local
    real (kind=rk) :: V(Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: policyA1(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: policyC(Tperiods, numPointsA, numPointsY, numAIME)
    integer :: policyL(Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: EV(Tperiods+1, numPointsA, numPointsY, numAIME);
    real (kind=rk) :: EdU(Tperiods,   numPointsA, numPointsY, numAIME);
    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    real (kind=rk) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) ::  meanL(Tperiods), meanA(Tperiods)
    real (kind=rk) :: yemp(Tperiods, numSims)
    real (kind=rk) :: AIME(Tperiods + 1,numSims)
    
    integer :: n

    !Set asset grid
    call getassetgrid( params, grids%maxInc, grids%Agrid)
    !solve
    call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, EdU, .FALSE. )
    !simulate
    call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp, AIME )
    do n=1,Tperiods
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
        meanA(n)=sum(real(Apath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
    end do

    gmm = dot_product(abs(meanL(32:32+23)-target(1,:)),abs(meanL(32:32+23)-target(1,:))) + &
          dot_product(abs(meanA(32:32+23)-target(2,:)),abs(meanA(32:32+23)-target(2,:)))
    
    end function
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Write to file

    subroutine writetofile(params, ypath, cpath, apath, vpath, lpath, yemp, AIME)
    implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: ypath(Tperiods, numSims) !income
    real (kind=rk), intent(in) :: cpath(Tperiods, numSims)  !consumption
    integer, intent(in) :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk), intent(in) :: vpath(Tperiods, numSims)  !value
    real (kind=rk), intent(in) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(in)  :: yemp(Tperiods, numSims)
    real (kind=rk), intent(in)  :: AIME(Tperiods + 1,numSims)
    !local
    integer :: n, requiredl 
    real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanA(Tperiods), meanC(Tperiods), meanY(Tperiods)
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods)
    
    if (params%system == 1 ) then !ifort
        do n=1,Tperiods
            meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
            !write (201, * ) meanL(n)
            meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
            !write (202, * ) meanV(n)
            meanA(n)=sum(real(apath(n,:),rk))/real(numSims,rk)
            !write (203, * ) meanA(n)
            meanC(n)=sum(real(cpath(n,:),rk))/real(numSims,rk)
            !write (204, * ) meanC(n)
            meanY(n)=sum(real(ypath(n,:),rk))/real(numSims,rk)
            !write (205, * ) meanY(n)
            meanYemp(n)=sum(real(yemp(n,:),rk))/real(numSims,rk)            
            !write (205, * ) meanY(n)
            meanAIME(n)=sum(real(AIME(n,:),rk))/real(numSims,rk)
        end do
        !L
        inquire (iolength=requiredl)  meanL
        open (unit=201, file='..\\out\lpath.txt', status='unknown',recl=requiredl, action='write')
        write (201, '(6E15.3)' ) meanL
        close( unit=201)
        !V
        inquire (iolength=requiredl)  meanV
        open (unit=202, file='..\\out\Vpath.txt', status='unknown',recl=requiredl, action='write')
        write (202, '(6E15.3)' ) meanV
        close( unit=202)
        !A
        inquire (iolength=requiredl)  meanA
        open (unit=203, file='..\\out\Apath.txt', status='unknown',recl=requiredl, action='write')
        write (203, '(6E15.3)' ) meanA
        close( unit=203)
        !C
        inquire (iolength=requiredl)  meanC
        open (unit=204, file='..\\out\Cpath.txt', status='unknown',recl=requiredl, action='write')
        write (204, '(6E15.3)' ) meanC
        close( unit=204)
        !Y
        inquire (iolength=requiredl)  meanY
        open (unit=205, file='..\\out\Ypath.txt', status='unknown',recl=requiredl, action='write')
        write (205, '(6E15.3)' ) meanY
        close( unit=205)

        inquire (iolength=requiredl)  meanYemp
        open (unit=206, file='..\\out\YempPath.txt', status='unknown',recl=requiredl, action='write')
        write (206, '(6E15.3)' ) meanYemp
        close( unit=206)
        
        inquire (iolength=requiredl)  meanAIME
        open (unit=207, file='..\\out\AIMEPath.txt', status='unknown',recl=requiredl, action='write')
        write (207, '(6E15.3)' ) meanAIME
        close( unit=207)        
        
    else !Gfort
        inquire (iolength=requiredl)  meanL
        open (unit=201, file='./out/lpath', status='unknown',recl=requiredl, action='write')
        write (201, * ) 'Header'

        inquire (iolength=requiredl)  meanV
        open (unit=202, file='./out/Vpath', status='unknown',recl=requiredl, action='write')
        write (202, * ) 'Header'

        inquire (iolength=requiredl)  meanA
        open (unit=203, file='./out/Apath', status='unknown',recl=requiredl, action='write')
        write (203, * ) 'Header'
    
        inquire (iolength=requiredl)  meanC
        open (unit=204, file='./out/Cpath', status='unknown',recl=requiredl, action='write')
        write (204, * ) 'Header'
    
        inquire (iolength=requiredl)  meany
        open (unit=205, file='./out/Ypath', status='unknown',recl=requiredl, action='write')
        write (205, * ) 'Header'

        inquire (iolength=requiredl)  meanYemp
        open (unit=206, file='./out/YempPath', status='unknown',recl=requiredl, action='write')
        write (206, * ) 'Header'

        inquire (iolength=requiredl)  meanAIME
        open (unit=207, file='./out/AIMEPath', status='unknown',recl=requiredl, action='write')
        write (207, * ) 'Header'
    

        do n=1,Tperiods
            meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
            write (201, * ) meanL(n)
            meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
            write (202, * ) meanV(n)
            meanA(n)=sum(real(apath(n,:),rk))/real(numSims,rk)
            write (203, * ) meanA(n)
            meanC(n)=sum(real(cpath(n,:),rk))/real(numSims,rk)
            write (204, * ) meanC(n)
            meanY(n)=sum(real(ypath(n,:),rk))/real(numSims,rk)
            write (205, * ) meanY(n)
            !meanYemp(n)=sum(real(yemp(n,:),rk))/real(meanL(n)*numSims,rk)
            meanYemp(n)=sum(real(yemp(n,:),rk))/real(numSims,rk)
            write (206, * ) meanYemp(n)
            meanAIME(n)=sum(real(AIME(n,:),rk))/real(numSims,rk)
            write (207, * ) meanAIME(n)
         end do
        close( unit=201)
        close( unit=202)
        close( unit=203)
        close( unit=204)
        close( unit=205)
        close( unit=206)
        close( unit=207)
    end if

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!get asset grid
    subroutine getassetgrid( params, maxInc, Agrid)
        implicit none
    !inputs
    type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) ::  maxInc(Tperiods)
    !outputs
    real (kind=rk), intent(out) :: Agrid(Tperiods+1, numPointsA)
    !local
    real (kind=rk) :: maxA(Tperiods+1), loggrid(numPointsA), span, test
    integer :: ixt, i

    !Set maximum assets
    maxA(1) = params%startA;
    do ixt = 2, Tperiods+1
        maxA(ixt) = (maxA(ixt - 1) + maxInc(ixt-1) ) * (1+params%r)
    end do

    !Create asset grid
    do ixt = 1, Tperiods+1
        span =  (log(1.0+log(1.0+log(1+maxA(ixt)))) - log(1.0+log(1.0+log(1.0))) )/ (numPointsA-1)
        loggrid = log(1.0+log(1.0+log(1.0))) + span*(/(i,i=0,numPointsA-1)/)
        Agrid(ixt, :) = (/(exp(exp(exp(loggrid(i))-1.0)-1.0)-1.0,i=1,numPointsA)/) !exp(exp(exp(loggrid)-1)-1)-1
    end do
    test = sum(Agrid(1, :))/size(Agrid(ixt, :))


    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!solve period
    subroutine solvePeriod(params, grids, Yin, A, AIMEin ,ixt, lbA1, EV1, benefit,policyA1,policyC,policyL,V)
    implicit none
    !input
    real(kind=rk), intent(in) :: Yin, A, lba1, EV1(:,:), benefit, AIMEin
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt
    !output
    real(kind=rk), intent(out) :: policyA1, policyC, V
    integer, intent(out) :: policyL
    !local
    integer :: ixl
    real(kind=rk) :: Y, negV, negVtemp, ubA1, policyA1temp, AIME

    negV = -huge(negv) !-log(0.0) !inf
    AIME = AIMEin
    do ixL = 0,(numPointsL-1),1           ! points of income choice
        ! Value of income and information for optimisation
        Y    = ixL*Yin+ (1-ixL)*benefit;
        !AIME only depends on earned income so add spousal
        !and pensions after calculating it
        call gross2net(params,Y,ixt,ixl,AIME)
        ! Next peridos AIME
        if (ixL==1 .AND.  ixt < spouseretire) then
                AIME =  Yin/ixt + AIME * (ixt-1)/ixt
        else if ((ixL==0 .AND.  ixt < spouseretire)) then
                AIME = AIME * (ixt-1)/ixt
        end if

        ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow

        ! Compute solution
        if (ubA1 - lbA1 < params%minCons) then         ! if liquidity constrained
            negVtemp = objectivefunc(params, grids,lbA1, A, Y,ixL,ixt, AIME,EV1);
            policyA1temp = lbA1;
        else                               ! if interior solution
            ![policyA1temp, negVtemp] = ...
            !    fminbnd(@(A1) objectivefunc(A1, A, Y,ixL,ixt, AIME), lbA1, ubA1, optimset('TolX',tol));
            negVtemp = golden_generic(lbA1, ubA1, policyA1temp, func,params%tol,.FALSE.)
        end if! if (ubA1 - lbA1 < minCons)
        if (negVtemp > negV) then
            negV = negVtemp
            policyA1=policyA1temp
            policyL=ixL
            policyC = A + Y - policyA1/(1+params%r)
            ! Store solution and its value
        end if
    end do
    !testC = policyC(ixt, ixA, ixY,ixAIME)
    V  = negV;
    contains
    function func(x)
    real (kind = rk), intent(in) :: x
    real (kind = rk) :: func
    func = objectivefunc(params, grids,x, A, Y,ixL,ixt, AIME,EV1)
    end function

    end subroutine

    end module routines
