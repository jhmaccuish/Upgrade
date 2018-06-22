    module routines
    use Header
    use routines_generic

    implicit none

    interface unpackArray
    module procedure unpackArrayInt
    module procedure unpackArrayReal
    end interface unpackArray

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
    integer :: t, i, workAge

    sig_inc = params%sigma/((1-params%rho**2)**0.5)

    !Why 3
    call tauchen(numPointsY,params%mu,params%rho,params%sigma,3,ly,incTransitionMrx)

    !Ygrid = exp(repmat(ly', T, 1)+repmat(polyval(delta,1:T)'-fc(1:T),1,numPointsY));
    !minInc = exp(repmat((-normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((-normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));
    !maxInc = exp(repmat((normBnd * sig_inc)', T, 1)+polyval(delta,1:T)');%exp(repmat((normBnd * sig_inc)', T, 1)+repmat(polyval(delta,1:T),1,numPointsY));

    params%pension = 107.45*52
    upper(1) = 0
    do t=1 , Tperiods
        workAge = startAge - 20 + t
        Ygrid(t,:)= exp(ly+params%delta(1)*workAge**2+params%delta(2)*workAge+params%delta(3)-fc(t))
        minInc(t) = exp((-normBnd * sig_inc)+params%delta(1)*workAge**2+params%delta(2)*workAge+params%delta(3)-fc(t))
        maxInc(t) = exp((normBnd * sig_inc)+params%delta(1)*workAge**2+params%delta(2)*workAge+params%delta(3)-fc(t))
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
    subroutine solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, show )
    implicit none

#ifdef mpi  !mpi
    include 'mpif.h'
#endif     
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    logical :: show

    !outputs
    real (kind=rk), intent(out) :: V(Tperiods+1, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(out) :: policyA1(Tperiods, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(out) :: policyC(Tperiods, numPointsA, numAIME, numPointsY)
    integer, intent(out) :: policyL(Tperiods, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(out) :: EV(Tperiods+1, numPointsA, numAIME, numPointsY);
    !local
    integer :: ixt, ixAIME, ixA, ixY, ixL
    real (kind=rk) :: negV, A, Y, lbA1, ubA1, AIME, EV1(numPointsA ,numAIME) ! Agrid1(numPointsA)
    real (kind=rk) :: AIME1grid(numAIME), policyA1temp, negVtemp, realisedV(numPointsY)
    integer:: indexBigN(2), indexSmalln(2), singleIndex, thisCoreStartStore, thisCoreEndStore
    integer :: numTasks, tasksPerCore, leftOverTasks, thisCoreStart, thisCoreEnd, count
    real (kind=rk) :: VecV(mpiDim*numPointsY), VecpolicyA1(mpiDim*numPointsY), VecpolicyC(mpiDim*numPointsY), VecEV(mpiDim*numPointsY);
    integer :: VecpolicyL(mpiDim*numPointsY)
    real (kind=rk) :: tempV(mpiDim*numPointsY), temppolicyA1(mpiDim*numPointsY), temppolicyC(mpiDim*numPointsY), tempEV(mpiDim*numPointsY);
    integer :: temppolicyL(mpiDim*numPointsY)
    real (kind=rk), allocatable :: LocV(:), LocpolicyA1(:), LocpolicyC(:), LocEV(:)
    integer, allocatable :: locpolicyL(:)
    integer :: locVecSize, otherDimP, testRank, thisCoreStartTest, thisCoreEndTest, start, finish

    !Test
    real (kind=rk) :: testC

    !----------------------------------------------------------------------------------!
    ! Initialize MPI
    !----------------------------------------------------------------------------------!

    indexBigN(1) = numPointsA
    indexBigN(2) = numAIME
    numTasks = numPointsA*numAIME

    tasksPerCore = int(numTasks/procSize)
    leftOverTasks = numTasks - tasksPerCore*procSize
    if (leftOverTasks.gt.0) tasksPerCore = tasksPerCore + 1
    thisCoreStart = rank*tasksPerCore + 1
    thisCoreEnd = min((rank+1)*tasksPerCore, numTasks)
#ifdef mpi
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif 

    !Set the terminal value function and expected value function to 0
    EV(Tperiods + 1, :,:,:)  = 0;          ! continuation value at T-1
    V(Tperiods + 1,:,:,:) = 0;

    !Initialise everything
    EV(:, :,:,:)  = 0          ! continuation value at T-1
    V(:,:,:,:) = 0
    policyL(:,:,:,:) = 0
    policyC(:,:,:,:) = 0
    policyA1(:,:,:,:) = 0


    do ixt=Tperiods,1, -1                               ! Loop from time T-1 to 1
        !Agrid1 = grids%Agrid(ixt + 1, :);               ! The grid on assets tomorrow
        AIME1grid = grids%AIMEgrid(ixt + 1, :);
        !!$omp parallel do default(shared) private(ixAIME,ixA,negV,lba1,EV1,A,AIME,y,negVtemp,policyA1temp,ixL,ubA1)
        do ixAIME = 1, numAIME
            do ixA = 1, numPointsA                   ! points on asset grid
                indexSmalln(1) = ixa
                indexSmalln(2) = ixAIME
                singleIndex = sumindex(indexBigN, indexSmalln, .FALSE.)

                if ((singleIndex .ge. thisCoreStart) .and. (singleIndex .le. thisCoreEnd)) then
                    if (ixt < stopwrok) then
                        !Although doesn't recieve income still need to loop around
                        !hypothetical income because participation doesn't effect
                        !earning potential
                        ! STEP 1. solve problem at grid points in assets, income + labour choices
                        ! ---------------------------------------------------------
                        do ixY = 1, numPointsY               ! points on income grid
                            !negV = -huge(negv) !-log(0.0) !inf
                            lbA1 = grids%Agrid(ixt + 1, 1);          ! lower bound: assets tomorrow
                            EV1  = EV(ixt + 1,:,:,ixY);  ! relevant section of EV matrix (in assets tomorrow)
                            call solvePeriod(params, grids, grids%Ygrid(ixt, ixY),grids%Agrid(ixt, ixA), grids%AIMEgrid(ixt,ixAIME), &
                                ixt, lbA1, EV1, grids%benefit(ixt), policyA1(ixt,ixA,ixAIME, ixY), &
                                policyC(ixt, ixA, ixAIME, ixY), policyL(ixt,ixA,ixAIME,ixY), V(ixt, ixA, ixAIME,ixY))
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
                        EV1  = EV(ixt + 1,:,:,1);  ! relevant section of EV matrix (in assets tomorrow)
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
                            policyA1(ixt,ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/))=policyA1temp;
                            policyL(ixt,ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/))=ixL;
                            ! Store solution and its value
                            policyC(ixt,ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/)) = A + Y - policyA1(ixt, ixA,ixAIME,1)/(1+params%r);
                        end if
                        !end
                        !testC = policyC(ixt, ixA, 1,ixAIME)
                        V(ixt, ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/))       = negV;
                        !dU(ixt, ixA, ixY)      = getmargutility(policyC(ixt, ixA, ixY),policyL(ixt,ixA,ixY));
                    end if

                    ! STEP 2. integrate out income today conditional on income
                    ! yesterday to get EV and EdU
                    ! --------------------------------------------------------
                    realisedV(:) = V(ixt, ixA,ixAIME,:);
                    !realiseddU(:,:) = dU(ixt, ixA, :);
                    do ixY = 1,numPointsY,1
                        EV(ixt, ixA,ixAIME,ixY)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV);
                        !EdU(ixt, ixA, ixY) = incTransitionMrx(ixY,:)*realiseddU;
                    end do !ixY
                end if
            end do!ixA
        end do!ixAIME
        !!$omp end parallel do
#ifdef mpi 
        !write(*,*) 2
        locVecSize = thisCoreEnd - thisCoreStart + 1
        otherDimP = numPointsY
        allocate(LocV(locVecSize*otherDimP),LocpolicyA1(locVecSize*otherDimP),LocpolicyC(locVecSize*otherDimP),LocEV(locVecSize*otherDimP),LocpolicyL(locVecSize*otherDimP))


        call unpackArrays(policyL(ixt, :, :, :), policyA1(ixt, :, :, :), policyC(ixt, :, :, :), V(ixt, :, :, :), EV(ixt, :, :, :), &
            locpolicyL, locpolicyA1, locpolicyC, locV, locEV, mpiDim,otherDimP,thisCoreStart,thisCoreEnd)

        !!write(*,*) rank, locVecSize, otherDimP
        !write(*,*) "PL", sum(policyL(ixt, :, :, :))
        !write(*,*) "VecL", sum(LocpolicyL)
        !call mpi_barrier(mpi_comm_world, ierror)
        !!if (ierror.ne.0) stop 'mpi problem180
        !
        !write(*,*) "PA", sum(policyA1(ixt, :, :, :))
        !write(*,*) "VecA", sum(LocpolicyA1)
        !call mpi_barrier(mpi_comm_world, ierror)
        !!if (ierror.ne.0) stop 'mpi problem180
        !
        !write(*,*) "PC", sum(policyC(ixt, :, :, :))
        !write(*,*) "VecC", sum(LocpolicyC)
        !call mpi_barrier(mpi_comm_world, ierror)
        !!if (ierror.ne.0) stop 'mpi problem180
        !
        !write(*,*) "V", sum(V(ixt, :, :, :))
        !write(*,*) "VecV", sum(LocV)
        !call mpi_barrier(mpi_comm_world, ierror)
        !!if (ierror.ne.0) stop 'mpi problem180
        !
        !write(*,*) "EV", sum(EV(ixt, :, :, :))
        !write(*,*) "VecEV", sum(LocEV)
        !



        call mpi_allgather(LocV(1), locVecSize*otherDimP, mpi_double_precision, VecV(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
        call mpi_allgather(LocpolicyA1(1), locVecSize*otherDimP, mpi_double_precision, VecpolicyA1(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
        call mpi_allgather(LocpolicyC(1), locVecSize*otherDimP, mpi_double_precision, VecpolicyC(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
        call mpi_allgather(LocEV(1), locVecSize*otherDimP, mpi_double_precision, VecEV(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
        !write(*,*) 3
        call mpi_allgather(LocpolicyL(1), locVecSize*otherDimP, mpi_integer, VecpolicyL(1), locVecSize*otherDimP, mpi_integer, mpi_comm_world, ierror)
        !       write(*,*) 8

        deallocate(LocV,LocpolicyA1,LocpolicyC,LocEV,LocpolicyL)

        tempV = VecV
        temppolicyA1 = VecpolicyA1
        temppolicyC = VecpolicyC
        tempEV = VecEV
        temppolicyL = VecpolicyL

        do testRank = 0, (procSize-1)
            thisCoreStartTest = testrank*tasksPerCore + 1
            thisCoreEndTest = min((testrank+1)*tasksPerCore, numTasks)
            locVecSize = thisCoreEndTest - thisCoreStartTest + 1
            allocate(LocV(locVecSize*otherDimP),LocpolicyA1(locVecSize*otherDimP),LocpolicyC(locVecSize*otherDimP),LocEV(locVecSize*otherDimP),LocpolicyL(locVecSize*otherDimP))

            start = (thisCoreStartTest-1)*otherDimP+1
            finish = thisCoreEndTest*otherDimP

            locV= tempV(start:finish)
            locpolicyA1 = temppolicyA1(start:finish)
            locpolicyC = temppolicyC(start:finish)
            locEV = tempEV(start:finish)
            locpolicyL = temppolicyL(start:finish)

            !Distribute the contigous section between the contigous section in each column corresponding to the non-splitting dimesnions of the array
            do count=1,otherDimP
                start = (count-1)*mpiDim+thisCoreStartTest
                finish = (count-1)*mpiDim+thisCoreEndTest
                VecV(start:finish) = locV((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
                VecpolicyA1(start:finish) = locpolicyA1((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
                VecpolicyC(start:finish) = locpolicyC((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
                VecEV(start:finish) = locEV((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
                VecpolicyL(start:finish) = locpolicyL((count-1)*locVecSize+1:(count-1)*locVecSize+locVecSize)
            end do
            deallocate(LocV,LocpolicyA1,LocpolicyC,LocEV,LocpolicyL)
        end do!


        V(ixt, :, :, :) = reshape(VecV, (/numPointsA, numAIME, numPointsY/))
        policyA1(ixt, :, :, :) = reshape(VecpolicyA1, (/numPointsA, numAIME, numPointsY/))
        policyC(ixt, :, :, :) = reshape(VecpolicyC, (/numPointsA, numAIME, numPointsY/))
        EV(ixt, :, :, :) = reshape(VecEV, (/numPointsA, numAIME, numPointsY/))
        policyL(ixt, :, :, :) =  reshape(VecpolicyL, (/numPointsA, numAIME, numPointsY/))

        call mpi_barrier(mpi_comm_world, ierror)
        if (ierror.ne.0) stop 'mpi problem180'

#endif   
        if (show .AND. rank==0) WRITE(*,*)  'Passed period ', ixt, ' of ',  Tperiods
        if (show .AND. rank==0 .AND. ixt<20) WRITE(*,*)  "Sum Labour Supply Policy is:", sum( policyL(ixt, :, :, :) )
        !print '( "Avergage Labour Supply Policy is:",f6.3)', sum( policyL(ixt, :, :, :) )/size( policyL(ixt, :, :, :) )
    end do!ixt
    !write(*,*) 4

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
    real (kind=rk) :: cons, VA1, mortal(Tperiods+1), temp(86), VB1
    ! Declare global we need this file have access to
    !mortal = (/1,2,&
    !        3,4/)

    temp = (/0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011071, 0.011907, 0.012807, &
        0.013676, 0.01475, 0.015818, 0.016846, 0.018079, 0.019343, 0.020659, 0.0218, 0.023505, 0.025202, 0.02696, &
        0.028831, 0.031017, 0.033496, 0.036024, 0.038912, 0.042054, 0.045689, 0.049653, 0.054036, 0.05886, 0.064093, &
        0.069636, 0.07533, 0.081069, 0.086912, 0.093067, 0.099807, 0.107483, 0.116125, 0.125196, 0.134361, 0.143881, &
        0.1542, 0.165675, 0.17842, 0.192363, 0.2117, 0.1672,0.1565, 0.1485,0.1459, 1.0/);
    mortal = temp(startAge-20+1:86)
    !Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    !total value (u(c) + b * VA1
    cons = A0  + Y - (A1)/(1+params%r);
    !VA1 = interp1(Agrid1,EV1,A1, interpMethod, 'extrap');
    ![X,Y] = meshgrid(Agrid1,AIME1grid);
    !VA1 = interp2(Agrid1, AIME1grid , EV1, A1, AIME);

    call linearinterp2_withextrap(grids%Agrid(ixP + 1, :), grids%AIMEgrid(ixP + 1, :), &
        numPointsA, numAIME, A1, AIME, VA1, 1, 1, EV1)
    !interp2(EV1, A1/Agrid1(20), AIME/AIME1grid(10));
    VB1 = params%thetab*((A1+params%K)**(1-params%gamma))/(1-params%gamma);
    objectivefunc = utility(params,cons,L) + params%beta * ((1- mortal(ixP))* VA1+mortal(ixP)*VB1);

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
    les=(L)*(1-params%hrsWrk )+(1-L);
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
    real (kind=rk), intent(in) :: policyA1(Tperiods, numPointsA, numAIME, numPointsY)
    !real (kind=rk), intent(in) :: policyC(Tperiods, numPointsA, numPointsY, numAIME)
    integer, intent(in) :: policyL(Tperiods, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(in) :: EV(Tperiods+1, numPointsA, numAIME, numPointsY);

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
    integer :: s, t, seed1, seed2, ios, idxA(1), idxAIME(1), idxY(1), workAge
    !real (kind=rk) :: temp2(numPointsY)
    !CHARACTER(len=255) :: cwd
    !CALL getcwd(cwd)

    integer :: seedIn

    !local
    real (kind=rk) ::  uniformRand(numSims), ltemp, lbA1, EV1(numPointsA ,numAIME)
    INTEGER :: n, i, uniformInt(numSims)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    if (fullLifeCycle) then
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
            uniformInt(i) = floor(uniformRand(i)*numPointsA)+1
            startingA(i) =grids%Agrid(1,uniformInt(i))
        end do
    else
        startingA = grids%initialAssets(1:size(startingA))
        do i =1, size(startingA)
            if (startingA(i) < 0.0) startingA(i) = 0.0
        end do
    end if
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

    !!$omp do private(s,t)
    do s = 1, numSims, 1                           ! loop through individuals
        ! Get all the incomes, recursively
        ly(1, s) = truncate(logy1(1, s), -normBnd*sig_inc,normBnd*sig_inc )
        y(1, s) = exp(ly(1, s)+params%delta(1)*(startAge - 20 + 1)**2+params%delta(2)*(startAge - 20 + 1)+params%delta(3)-grids%fc(1))
        do t = 1,Tperiods-1,1                              ! loop through time periods for a particular individual
            workAge = startAge - 20 + t
            ly(t+1, s) = (1 -params%rho) * params%mu + params%rho * ly(t, s) + e(t + 1, s)
            ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
            y(t+1, s) = exp( ly(t+1, s) + params%delta(1)*(workAge+1)**2+params%delta(2)*(workAge+1)+params%delta(3)-grids%fc(t) )
        end do ! t
    end do! s
    !!$omp end do

    !!$omp do private(s,t,idxY, idxA, idxAIME,lbA1,ev1)
    do s = 1,numSims,1
        a(1, s) = startingA(s)
        AIME(1,s)=startAIME
        do t = 1,Tperiods,1                              ! loop through time periods for a particular individual
            !Should improve this nearest neigbhour
            idxY=minloc(abs(y(t, s)-grids%Ygrid(t,:)))
            idxA=minloc(abs(a(t, s)-grids%Agrid(t,:)))
            idxAIME=minloc(abs(AIME(t, s)-grids%AIMEgrid(t,:)))
            l(t,s)=policyL(t,idxA(1),idxAIME(1),idxY(1))
            yemp(t,s) = y(t, s)

            call linearinterp3_withextrap(grids%Agrid(t,:), grids%AIMEgrid(t,:), grids%Ygrid(t, :), &
                numPointsA, numAIME, numPointsY, a(t, s), AIME(t, s), y(t, s), ltemp,  real(policyL(t, :, :,:),rk))
            if (abs(l(t,s) - ltemp) > 0.01) then
                lbA1 = grids%Agrid(t + 1, 1);          ! lower bound: assets tomorrow
                ev1 = EV(t + 1,:,:,idxY(1))
                call solvePeriod(params, grids, y(t, s), a(t, s), AIME(t, s) ,t, lbA1, ev1, &
                    grids%benefit(t),a(t+1, s),c(t, s),l(t,s),v(t  , s))
            else
                !a(t+1, s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tA1, a(t, s), (y(t, s)));
                !interp3(Agrid(t,:)', Ygrid(t, :)', AIMEgrid(t,:)',tA1, a(t, s), y(t, s),AIME(t, s));
                call linearinterp3_withextrap(grids%Agrid(t,:), grids%AIMEgrid(t,:), grids%Ygrid(t, :), &
                    numPointsA, numAIME, numPointsY, a(t, s), AIME(t, s), y(t, s), a(t+1, s),  policyA1(t, :, :,:))
                !v(t  , s) =  interp2D(Agrid(t,:)', Ygrid(t, :)', tV , a(t, s), (y(t, s)));
                !interp3(Agrid(t,:)', Ygrid(t, :)', AIMEgrid(t,:)',tV, a(t, s), y(t, s),AIME(t, s));
                call linearinterp3_withextrap(grids%Agrid(t,:), grids%AIMEgrid(t,:), grids%Ygrid(t, :),&
                    numPointsA, numAIME, numPointsY,  a(t, s), AIME(t, s), y(t, s), v(t  , s),  EV(t, :, :,:))
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
            end if
            if  (t < Tretire) then ! spouseretire
                if (l(t,s) .EQ. 0) then
                    AIME(t+1, s) =   AIME(t, s) * (t-1)/t
                else
                    AIME(t+1, s) =  y(t,s)/t + AIME(t, s) * (t-1)/t
                end if
            else
                AIME(t+1, s) =   AIME(t, s)
            end if


            call gross2net(params,y(t, s),t,l(t,s), AIME(t, s))
            ! tomorrow's optimal assets
            c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+params%r))

        end do   !t
    end do! s
    !yemp = yemp*l
    !!$omp end do

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
            if (ixL==0 ) then !.AND. ixt >= 50
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
#ifdef mpi  !mpi
    include 'mpif.h'
#endif  

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
    call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, .FALSE. )
    if (rank==0) then
        !simulate
        call simWithUncer(params, grids, policyA1,policyL,EV, ypath, cpath, apath, vpath, lpath, yemp, AIME )
        do n=1,Tperiods
            meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
            meanA(n)=sum(real(Apath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
        end do
        if (fullLifeCycle) then
            gmm = dot_product(abs(meanL(32:32+23)-target(1,:)),abs(meanL(32:32+23)-target(1,:)))! + &
                !dot_product(abs(meanA(32:32+23)-target(2,:)),abs(meanA(32:32+23)-target(2,:)))
            else
        gmm = dot_product(abs(meanL(33 - (StartAge-20):33 - (StartAge-20)+23)-target(1,:)),abs(meanL(33 - (StartAge-20):33 - (StartAge-20)+23)-target(1,:)))! + &
            !dot_product(abs(meanA(33 - (StartAge-20):33 - (StartAge-20)+23)-target(2,:)),abs(meanA(33 - (StartAge-20):33 - (StartAge-20)+23)-target(2,:)))
        end if

    end if

#ifdef mpi 
    call MPI_Bcast( gmm,   1,   mpi_double_precision, 0, mpi_comm_world, ierror)
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif
    !write (*,*) gmm
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
    integer :: n, requiredl , i
    real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanA(Tperiods), meanC(Tperiods), meanY(Tperiods), medianA
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), meanPoor(Tperiods), meanRich(Tperiods)
    integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    !real(kind=rk),DIMENSION(Tperiods:), ALLOCATABLE ::
    integer :: rich, poor
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
        medianA = median(apath(Tretire,:))
        rich = 0
        poor = 0
        do n=1,numSims
            if (apath(Tretire,n) <= medianA) then
                poor = poor + 1
                lpoor(:,poor) = lpath(:,n)
            else
                rich = rich + 1
                lrich(:,rich) = lpath(:,n)
            end if
        end do

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

        inquire (iolength=requiredl)  meanPoor
        open (unit=208, file='./out/PoorPath', status='unknown',recl=requiredl, action='write')
        write (208, * ) 'Header'

        inquire (iolength=requiredl)  meanRich
        open (unit=209, file='./out/RichPath', status='unknown',recl=requiredl, action='write')
        write (209, * ) 'Header'

        inquire (iolength=requiredl)  meanRich
        open (unit=210, file='./out/ldata', status='unknown',recl=requiredl, action='write')


        inquire (iolength=requiredl)  meanRich
        open (unit=211, file='./out/adata', status='unknown',recl=requiredl, action='write')

        inquire (iolength=requiredl)  meanRich
        open (unit=212, file='./out/ydata', status='unknown',recl=requiredl, action='write')

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
            meanRich(n)=sum(real(lrich(n,:),rk))/real(numSims/2,rk) !size(lpath(n,:))
            write (209, * ) meanRich(n)
            meanPoor(n)=sum(real(lpoor(n,:),rk))/real(numSims/2,rk) !size(lpath(n,:))
            write (208, * ) meanPoor(n)
            !            if (n .GE. 32 .AND. n .LE. 55) then
            !            do i=1,numsims
            !                write (210, * ) lpath(n,i)
            !                write (211, * ) apath(n,i)
            !                write (212, * ) lpath(n,i)*ypath(n,i)
            !            end do
            !            end if
        end do
        do i=1,numsims
            do n=32,55
                write (210, * ) lpath(n,i)
                write (211, * ) apath(n,i)
                write (212, * ) lpath(n,i)*ypath(n,i)
            end do
        end do

        close( unit=201)
        close( unit=202)
        close( unit=203)
        close( unit=204)
        close( unit=205)
        close( unit=206)
        close( unit=207)
        close( unit=208)
        close( unit=209)
        close( unit=210)
        close( unit=211)
        close( unit=212)
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
        if (ixL==1 .AND.  ixt < Tretire ) then ! spouseretire
            AIME =  Yin/ixt + AIME * (ixt-1)/ixt
        else if ((ixL==0 .AND.  ixt < Tretire )) then !spouseretire
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
    !!-----------------------------------------------------------------------------------------------------------!
    ! Unpack all arrays
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine unpackArrays(pL,pA,pC,V,eV,vecPl,vecPA,vecPC,vecV, vecEV,dim1,dim2,thisCoreStart,thisCoreEnd)
    implicit none

    !inputs
    integer, intent(in) :: pL( :, :, :), dim1,dim2, thisCoreStart,thisCoreEnd
    real(kind=rk), intent(in) :: pA( :, :, :), pC( :, :, :)
    real(kind=rk), intent(in) :: V( :, :, :), eV( :, :, :)
    !outputs
    integer, intent(out) :: vecPl(:)
    real(kind=rk), intent(out) :: vecPA(:), vecPC(:)
    real(kind=rk), intent(out) :: vecV(:), vecEV(:)

    call unpackArray(pl,vecPl,dim1,dim2,thisCoreStart,thisCoreEnd)
    call unpackArray(pA,vecPA,dim1,dim2,thisCoreStart,thisCoreEnd)
    call unpackArray(pC,vecPC,dim1,dim2,thisCoreStart,thisCoreEnd)
    call unpackArray(v,vecV,dim1,dim2,thisCoreStart,thisCoreEnd)
    call unpackArray(eV,vecEV,dim1,dim2,thisCoreStart,thisCoreEnd)

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    ! Unpack integer array
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine unpackArrayInt(Array,vec,dim1,dim2,thisCoreStart,thisCoreEnd)
    implicit none
    !inputs
    integer, intent(in) :: Array( :, :, :), dim1,dim2,thisCoreStart,thisCoreEnd

    !outputs
    integer, intent(out) :: vec(:)

    !local
    integer, allocatable :: mat(:,:)

    allocate( mat(dim1,dim2))

    mat = reshape(Array,(/dim1,dim2/))
    vec = reshape(mat(thisCoreStart:thisCoreEnd,:),(/(thisCoreEnd-thisCoreStart+1)*dim2/))

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    ! Unpack integer array
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine unpackArrayReal(Array,vec,dim1,dim2,thisCoreStart,thisCoreEnd)
    implicit none
    !inputs
    real(kind=rk), intent(in) :: Array( :, :, :)
    integer, intent(in) :: dim1,dim2,thisCoreStart,thisCoreEnd

    !outputs
    real(kind=rk), intent(out) :: vec(:)

    !local
    real(kind=rk), allocatable :: mat(:,:)

    allocate( mat(dim1,dim2))

    mat = reshape(Array,(/dim1,dim2/))
    vec = reshape(mat(thisCoreStart:thisCoreEnd,:),(/(thisCoreEnd-thisCoreStart+1)*dim2/))

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    !Intialise gues for Nedler-Mead Algorithm
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine initialGuess(rank,params,grids,moments,p,y)
    implicit none
    !inputs
    integer, intent(in) :: rank
    real(kind=rk), intent(in) :: moments(2,24)

    !changing
    type (structparamstype), intent(inout) :: params
    type (gridsType), intent(inout) :: grids

    !outpus
    real(kind=rk), intent(out) :: y(dimEstimation+1), p(dimEstimation+1,dimEstimation)

    !local
    integer :: i, seedIn, n
    real (kind=rk) ::  uniformRand(dimEstimation+1)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    if (fullLifeCycle) then
        if (rank==0) then
            print '("Guess 1")'
        end if
        params%nu =      0.38022456150504280
        params%beta =  0.97235322545400193
        params%gamma =  2.092041817380061
        params%db(1) = 0.91387622628345655
        params%db(2) = -4.7393420983952571E-005

        p(1,1) = params%nu
        p(1,2) = params%beta
        p(1,3) = params%gamma
        p(1,4) = params%db(1)
        p(1,5) = params%db(2)
        y(1) = gmm_criteria(p(1,:))

        !!$omp section
        if (rank==0) then
            print '("Guess 2")'
        end if
        p(2,1) = 0.4637
        p(2,2) = 0.970
        P(2,3) = 1
        p(2,4) = params%db(1)*0.9
        p(2,5) = params%db(2)*1.2
        y(2) = gmm_criteria(p(2,:))

        !!$omp section
        if (rank==0) then
            print '("Guess 3")'
        end if
        p(3,1) = 0.322
        p(3,2) = 0.9843
        P(3,3) = 2
        p(3,4) = params%db(1)*1.1
        p(3,5) = params%db(2)*0.7
        y(3) = gmm_criteria(p(3,:))

        !!$omp section
        if (rank==0) then
            print '("Guess 4")'
        end if
        p(4,1) = 0.55
        p(4,2) = 0.96
        P(4,3) = 0.5
        p(4,4) = params%db(1)*1.3
        p(4,5) = params%db(2)*0.95
        y(4) = gmm_criteria(p(4,:))

        !!$omp section
        if (rank==0) then
            print '("Guess 5")'
        end if
        p(5,1) = 0.15
        p(5,2) = 0.9999
        P(5,3) = 4
        p(5,4) = params%db(1)*0.85
        p(5,5) = params%db(2)*1.15
        y(5) = gmm_criteria(p(5,:))

        !!$omp section
        if (rank==0) then
            print '("Guess 6")'
        end if
        p(6,1) = 0.27
        p(6,2) = 0.986
        P(6,3) = 0.9
        p(6,4) = params%db(1)*0.99
        p(6,5) = params%db(2)*0.87
        y(6) = gmm_criteria(p(6,:))

        !!$omp end sections
        !!$omp end parallel
    else

        seedIn = 16101988
        !Set seed
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        seed = seedIn * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)

        !!get uniform random number
        CALL RANDOM_NUMBER(uniformRand)
        uniformRand = -0.5+uniformRand
        do i = 1, dimEstimation+1
            if (rank==0) then
                write (*,*) "Guess ", i
            end if            
            p(i,1) = params%nu*(1+uniformRand(i))
            p(i,2) = (0.95+0.1*uniformRand(i))!params%beta*(1+uniformRand(i))
            p(i,3) = params%gamma*(1+uniformRand(i))
            p(i,4) = params%db(1)*(1+uniformRand(i))
            p(i,5) = params%db(2)*(1+uniformRand(i))
            y(i) = gmm_criteria(p(i,:))
        end do
    end if
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

    end subroutine
    end module routines
