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
    subroutine getIncomeGrid(params, grids)
    implicit none

    !Changing
    type (structparamstype), intent(inout) :: params
    type (gridsType), intent(inout) :: grids

    !outputs
    !real (kind=rk) :: YgridOut(:,:), incTransitionMrxOut(:,:), maxInc(:), AIMEgrid(:,:), benefit(:)

    !local
    real (kind=rk) :: sig_inc, ly(numPointsProd), upper(numPointsType,Tperiods+1), a !Q(numPointsY,numPointsY)
    real (kind=rk) :: Ygrid(numPointsType,Tperiods,numPointsProd), incTransitionMrx(numPointsProd,numPointsProd)
    integer :: t, i, workAge, typeSim

    sig_inc = params%sigma/((1-params%rho**2)**0.5)

    call tauchen(numPointsProd,params%mu,params%rho,params%sigma,normBnd,ly,incTransitionMrx)

    params%pension = 107.45*52
    upper(:,1) = 0
    do typeSim = 1, numPointsType
        do t=1 , Tperiods
            workAge = startAge - 20 + t
            Ygrid(typeSim,t,:)= exp(ly+params%delta(typeSim,1)*workAge**2+params%delta(typeSim,2)*workAge+params%delta(typeSim,3)-grids%fc(t))
            !if (typeSim==numPointsType) then
            grids%maxInc(typeSim,t) = exp((normBnd * sig_inc)+params%delta(typeSim,1)*workAge**2+params%delta(typeSim,2)*workAge+params%delta(typeSim,3)-grids%fc(t))
            !This should be discounted to average not sum
            upper(typeSim,t+1) = upper(typeSim,t) + grids%maxInc(typeSim,t)
            if (t <=spouseretire) then
                a = (upper(typeSim,t+1)/t)/(numAIME-1)
                grids%AIMEgrid(typeSim,t,:) = a*(/(i,i=0,numAIME-1)/)
            else
                grids%AIMEgrid(typeSim,t,:) = grids%AIMEgrid(typeSim,t-1,:)
            end if
            !if (t <=5) then
            !    grids%benefit(t) = 57.90*52
            !else
            grids%benefit(t) = 0!73.10*52
            !end if
            !end if
            !grids%YgridExtrap(typeSim,t,1) = 0.0
            grids%YgridExtrap(typeSim,t,:) = Ygrid(typeSim,t,:)
        end do
        call addUnemploymentshock(Ygrid(typeSim,:,:), incTransitionMrx,grids,typeSim)
        grids%AIMEgrid(typeSim,Tperiods+1,:) = grids%AIMEgrid(typeSim,Tperiods,:)
    end do


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

#ifdef mpi
    include 'mpif.h'
#endif     
    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    logical :: show

    !outputs
    real (kind=rk), intent(out) :: V(Tperiods+1,numPointsType, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(out) :: policyA1(Tperiods,numPointsType, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(out) :: policyC(Tperiods,numPointsType, numPointsA, numAIME, numPointsY)
    integer, intent(out) :: policyL(Tperiods, numPointsType,numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(out) :: EV(Tperiods+1,numPointsType, numPointsA, numAIME, numPointsY);
    !local
    integer :: ixt, ixAIME, ixA, ixY, ixL
    real (kind=rk) :: negV, A, Y, lbA1, ubA1, AIME, EV1(numPointsA ,numAIME) ! Agrid1(numPointsA)
    real (kind=rk) :: AIME1grid(numAIME), policyA1temp, negVtemp, realisedV(numPointsY)
    integer:: indexBigN(2), indexSmalln(2), singleIndex, typeSim
    integer :: numTasks, tasksPerCore, leftOverTasks, thisCoreStart, thisCoreEnd, requiredl
    character(len=1024) :: outFile
#ifdef mpi
    integer:: thisCoreStartStore, thisCoreEndStore, count
    real (kind=rk) :: VecV(mpiDim*numPointsY), VecpolicyA1(mpiDim*numPointsY), VecpolicyC(mpiDim*numPointsY), VecEV(mpiDim*numPointsY);
    integer :: VecpolicyL(mpiDim*numPointsY)
    real (kind=rk) :: tempV(mpiDim*numPointsY), temppolicyA1(mpiDim*numPointsY), temppolicyC(mpiDim*numPointsY), tempEV(mpiDim*numPointsY);
    integer :: temppolicyL(mpiDim*numPointsY)
    real (kind=rk), allocatable :: LocV(:), LocpolicyA1(:), LocpolicyC(:), LocEV(:)
    integer, allocatable :: locpolicyL(:)
    integer :: locVecSize, otherDimP, testRank, thisCoreStartTest, thisCoreEndTest, start, finish

    !Test
    real (kind=rk) :: testC

#endif   

    !----------------------------------------------------------------------------------!
    ! Initialize MPI
    !----------------------------------------------------------------------------------!

    indexBigN(1) = numPointsA
    indexBigN(2) = numAIME
    numTasks = numPointsA*numAIME !*numPointsType

    tasksPerCore = int(numTasks/procSize)
    leftOverTasks = numTasks - tasksPerCore*procSize
    if (leftOverTasks.gt.0) tasksPerCore = tasksPerCore + 1
    thisCoreStart = rank*tasksPerCore + 1
    thisCoreEnd = min((rank+1)*tasksPerCore, numTasks)
#ifdef mpi
    !write (*,*) 1
    !write (*,*) numTasks, tasksPerCore, leftOverTasks
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif 

    !Set the terminal value function and expected value function to 0
    EV(Tperiods + 1,:,:,:,:)  = 0;          ! continuation value at T-1
    V(Tperiods + 1,:,:,:,:) = 0;

    !Initialise everything
    EV(:,:, :,:,:)  = 0          ! continuation value at T-1
    V(:,:,:,:,:) = 0
    policyL(:,:,:,:,:) = 0
    policyC(:,:,:,:,:) = 0
    policyA1(:,:,:,:,:) = 0

    do typeSim = 1, numPointsType
        if (show .AND. rank==0) WRITE(*,*)  'Solving for type ', typeSim, ' of ',  numPointsType
        do ixt=Tperiods,1, -1                               ! Loop from time T-1 to 1
            AIME1grid = grids%AIMEgrid(typeSim,ixt + 1, :);
            do ixAIME = 1, numAIME
                do ixA = 1, numPointsA                   ! points on asset grid
                    indexSmalln(1) = ixa
                    indexSmalln(2) = ixAIME
                    !indexSmalln(1) = typeSim
                    singleIndex = sumindex(indexBigN, indexSmalln, .FALSE.)

                    if ((singleIndex .ge. thisCoreStart) .and. (singleIndex .le. thisCoreEnd)) then
                        if (ixt < stopwrok) then
                            !Although doesn't recieve income still need to loop around
                            !hypothetical income because participation doesn't effect
                            !earning potential
                            ! STEP 1. solve problem at grid points in assets, income + labour choices
                            ! ---------------------------------------------------------
                            do ixY = 1, numPointsY               ! points on income grid
                                lbA1 = grids%Agrid(typeSim,ixt + 1, 1);          ! lower bound: assets tomorrow
                                EV1  = EV(ixt + 1,typeSim,:,:,ixY);  ! relevant section of EV matrix (in assets tomorrow)
                                call solvePeriod(params, grids, grids%Ygrid(typeSim,ixt, ixY),grids%Agrid(typeSim,ixt, ixA), grids%AIMEgrid(typeSim,ixt,ixAIME), &
                                    ixt, typeSim, lbA1, EV1, grids%benefit(ixt), policyA1(ixt,typeSim,ixA,ixAIME, ixY), &
                                    policyC(ixt,typeSim, ixA, ixAIME, ixY), policyL(ixt,typeSim,ixA,ixAIME,ixY), V(ixt,typeSim, ixA, ixAIME,ixY))
                            end do
                        else
                            negV = -huge(negv)
                            ixL = 0
                            ! Value of income and information for optimisation
                            A    = grids%Agrid(typeSim,ixt, ixA)            ! assets today
                            Y    = grids%benefit(ixt)
                            AIME = grids%AIMEgrid(typeSim,ixt,ixAIME)
                            call gross2net(params,Y,ixt,ixl,typeSim,AIME)
                            if (y < 0) then
                                print *, 'Error Y< 0'
                            end if

                            lbA1 = grids%Agrid(typeSim,ixt + 1, 1);          ! lower bound: assets tomorrow
                            ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow
                            EV1  = EV(ixt + 1,typeSim,:,:,1);  ! relevant section of EV matrix (in assets tomorrow)
                            ! Compute solution
                            if (ubA1 - lbA1 < params%minCons) then         ! if liquidity constrained
                                negVtemp = objectivefunc(params, grids,lbA1, A, Y,ixL,ixt,typeSim, AIME,EV1);
                                policyA1temp = lbA1;
                            else                               ! if interior solution
                                negVtemp = golden_generic(lbA1, ubA1, policyA1temp, func,params%tol,.FALSE.)
                            end if! if (ubA1 - lbA1 < minCons)
                            if (negVtemp > negV) then
                                negV = negVtemp;
                                policyA1(ixt,typeSim,ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/))=policyA1temp;
                                policyL(ixt,typeSim,ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/))=ixL;
                                ! Store solution and its value
                                policyC(ixt,typeSim,ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/)) = A + Y - policyA1(ixt,typeSim, ixA,ixAIME,1)/(1+params%r);
                            end if
                            V(ixt,typeSim, ixA,ixAIME,(/(ixY, ixY=1, numPointsY)/))       = negV;
                        end if

                        ! STEP 2. integrate out income today conditional on income
                        ! yesterday to get EV and EdU
                        ! --------------------------------------------------------
                        realisedV(:) = V(ixt,typeSim, ixA,ixAIME,:);
                        do ixY = 1,numPointsY,1
                            EV(ixt,typeSim, ixA,ixAIME,ixY)  = dot_product(grids%incTransitionMrx(ixY,:),realisedV);
                        end do !ixY
                    end if
                end do!ixA
            end do!ixAIME

#ifdef mpi 
            !write (*,*) 2
            locVecSize = thisCoreEnd - thisCoreStart + 1
            otherDimP = numPointsY
            allocate(LocV(locVecSize*otherDimP),LocpolicyA1(locVecSize*otherDimP),LocpolicyC(locVecSize*otherDimP),LocEV(locVecSize*otherDimP),LocpolicyL(locVecSize*otherDimP))
            !write (*,*) 3
            call unpackArrays(policyL(ixt,typeSim, :, :, :), policyA1(ixt,typeSim, :, :, :), policyC(ixt,typeSim, :, :, :), V(ixt,typeSim, :, :, :), EV(ixt,typeSim, :, :, :), &
                locpolicyL, locpolicyA1, locpolicyC, locV, locEV, mpiDim,otherDimP,thisCoreStart,thisCoreEnd)
            !write (*,*) 4
            !write (*,*)  locVecSize, otherDimP, mpiDim, numPointsY
            call mpi_allgather(LocV(1), locVecSize*otherDimP, mpi_double_precision, VecV(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
            call mpi_allgather(LocpolicyA1(1), locVecSize*otherDimP, mpi_double_precision, VecpolicyA1(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
            call mpi_allgather(LocpolicyC(1), locVecSize*otherDimP, mpi_double_precision, VecpolicyC(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
            call mpi_allgather(LocEV(1), locVecSize*otherDimP, mpi_double_precision, VecEV(1), locVecSize*otherDimP, mpi_double_precision, mpi_comm_world, ierror)
            call mpi_allgather(LocpolicyL(1), locVecSize*otherDimP, mpi_integer, VecpolicyL(1), locVecSize*otherDimP, mpi_integer, mpi_comm_world, ierror)
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
            V(ixt,typeSim, :, :, :) = reshape(VecV, (/numPointsA, numAIME, numPointsY/))
            policyA1(ixt,typeSim, :, :, :) = reshape(VecpolicyA1, (/numPointsA, numAIME, numPointsY/))
            policyC(ixt,typeSim, :, :, :) = reshape(VecpolicyC, (/numPointsA, numAIME, numPointsY/))
            EV(ixt,typeSim, :, :, :) = reshape(VecEV, (/numPointsA, numAIME, numPointsY/))
            policyL(ixt,typeSim, :, :, :) =  reshape(VecpolicyL, (/numPointsA, numAIME, numPointsY/))

            call mpi_barrier(mpi_comm_world, ierror)
            if (ierror.ne.0) stop 'mpi problem180'
#endif  

            if (show .AND. rank==0) WRITE(*,*)  'Passed period ', ixt, ' of ',  Tperiods
            if (show .AND. rank==0 .AND. ixt<= spouseretire) WRITE(*,*)  'Sum Policy L', sum(policyL(ixt,typeSim, :, :, :))
        end do!ixt
        if (rank==0 .AND. intermediateToFile) then
            !L
            write (outFile, *), trim("..\\temp\polL"),typeSim,".txt"
            outfile=ADJUSTL(outfile)
            inquire (iolength=requiredl) policyL(1,typeSim, :, 1, 2)
            open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
            write (201, *) policyL(1,typeSim, :, 1, 2)
            close( unit=201)
            !A
            write (outFile, *), trim("..\\temp\polA"),typeSim,".txt"
            outfile=ADJUSTL(outfile)
            inquire (iolength=requiredl) policyA1(1,typeSim, :, 1, 2)
            open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
            write (201, *) policyA1(1,typeSim, :, 1, 2)
            close( unit=201)
            !EV
            write (outFile, *), trim("..\\temp\EV"),typeSim,".txt"
            outfile=ADJUSTL(outfile)
            inquire (iolength=requiredl) EV(1,typeSim, :, 1, 2)
            open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
            write (201, *) EV(1,typeSim, :, 1, 2)
            close( unit=201)
        end if

    end do
    contains
    function func(x)
    real (kind = rk), intent(in) :: x
    real (kind = rk) :: func
    func = objectivefunc(params, grids,x, A, Y,ixL,ixt,typeSim, AIME,EV1)
    end function

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Objective function
    function objectivefunc(params, grids, A1, A0, Y,L,ixP,ixType, AIME, EV1)
    implicit none
    ! This function returns the following quantity:
    ! - (u(c) +  b V( A1))
    ! where c is calculated from today's assets and tomorrow's assets

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    real (kind=rk), intent(in) :: A1, A0, Y, AIME, EV1(:,:)
    integer, intent(in) :: ixP, L, ixType
    !ouptut
    real (kind=rk) :: objectivefunc
    !local
    real (kind=rk) :: cons, VA1, VB1

    !Get tomorrow's consumption (cons), the value of left over assets (VA1) and
    !total value (u(c) + b * VA1
    cons = A0  + Y - (A1)/(1+params%r);

    call linearinterp2_withextrap(grids%Agrid(ixType,ixP + 1, :), grids%AIMEgrid(ixType,ixP + 1, :), &
        numPointsA, numAIME, A1, AIME, VA1, 1, 1, EV1)
    VB1 = params%thetab*((A1+params%K)**(1-params%gamma))/(1-params%gamma);
    objectivefunc = utility(params,cons,L) + params%beta * ((1- grids%mortal(ixP))* VA1+grids%mortal(ixP)*VB1);

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
    subroutine simWithUncer(params, grids, policyA1,policyL,EV,  yex, c, a, v, l, y, AIME )
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    real (kind=rk), intent(in) :: policyA1(Tperiods,  numPointsType, numPointsA, numAIME, numPointsY)
    integer, intent(in) :: policyL(Tperiods,  numPointsType, numPointsA, numAIME, numPointsY)
    real (kind=rk), intent(in) :: EV(Tperiods+1,  numPointsType, numPointsA, numAIME, numPointsY);
    real (kind=rk), intent(in) :: yex(Tperiods, numSims)

    !outputs
    real (kind=rk), intent(out) :: y(Tperiods, numSims) !income
    real (kind=rk), intent(out) :: c(Tperiods, numSims)  !consumption
    integer, intent(out) :: l(Tperiods, numSims) !labour supply
    real (kind=rk), intent(out) :: v(Tperiods, numSims)  !value
    real (kind=rk), intent(out) :: a(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(out) :: AIME(Tperiods + 1,numSims)

    !local
    real (kind=rk) :: startingA(numSims), startAIME, Aval,AIMEval
    integer :: s, t, idxA(1), idxAIME(1), idxY(1), workAge
    integer :: seedIn, Lcube(8), typeSim
    real (kind=rk) ::  uniformRand(numSims), ltemp, lbA1, EV1(numPointsA ,numAIME)
    INTEGER :: n, i, uniformInt(numSims)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    integer :: unemployed !1 employed, 2 unemploeyed

    real (kind=rk):: LocpolicyA1(Tperiods,  numPointsType, numPointsA, numAIME, numPointsProd,2)
    integer :: LocpolicyL(Tperiods,  numPointsType, numPointsA, numAIME, numPointsProd,2)
    real (kind=rk) :: LocEV(Tperiods+1,  numPointsType, numPointsA, numAIME, numPointsProd,2);
    real (kind=rk) :: gridY(numPointsProd), Yval

    LocpolicyA1(:,:,:,:,:,1) =  policyA1(:,:,:,:,(/(i,i=1,numPointsY-1,2)/))
    LocpolicyA1(:,:,:,:,:,2) =  policyA1(:,:,:,:,(/(i,i=2,numPointsY,2)/))
    LocpolicyL(:,:,:,:,:,1) =  policyL(:,:,:,:,(/(i,i=1,numPointsY-1,2)/))
    LocpolicyL(:,:,:,:,:,2) =  policyL(:,:,:,:,(/(i,i=2,numPointsY,2)/))
    LocEV(:,:,:,:,:,1) =  EV(:,:,:,:,(/(i,i=1,numPointsY-1,2)/))
    LocEV(:,:,:,:,:,2) =  EV(:,:,:,:,(/(i,i=2,numPointsY,2)/))
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
            startingA(i) =grids%Agrid(1,1,uniformInt(i))
        end do
    else
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        seed(n) = 16101988
        CALL RANDOM_SEED(PUT = seed)
        !!get uniform random number
        CALL RANDOM_NUMBER(uniformRand)
        do i = 1, numSims
            startingA(i) = grids%initialAssets(floor(uniformRand(i)*numSims)+1)
            if (startingA(i) < 0.0) startingA(i) = 0.0
        end do
        !startingA = grids%initialAssets(1:size(startingA))
        !do i =1, size(startingA)
            !if (startingA(i) < 0.0) startingA(i) = 0.0
        !end do
    end if
    startAIME = 15500

    do s = 1,numSims,1
        if (numPointsType == 1) then
            typeSim = 1
        else
            if (real(s,rk)/numSims <0.16) then
                typeSim = 1
            else if (real(s,rk)/numSims <0.48) then
                typeSim = 2
            else if (real(s,rk)/numSims <0.59) then
                typeSim =3
            else
                typeSim =4
            end if
        end if
        a(1, s) = startingA(s)
        AIME(1,s)=startAIME
        do t = 1,Tperiods,1                              ! loop through time periods for a particular individual
            if (yex(t,s) <= 0.0001 ) then
                unemployed = 2
                if (numPointsProd==5) then
                    !gridY = (/1,2,3,4,5/)
                else if (numPointsProd==10) then
                    gridY = (/1,2,3,4,5,6,7,8,9,10/)
                end if
                yval = idxY(1)
            else
                unemployed = 1
                !Everyone start employed so this is OK
                idxY=minloc(abs(yex(t,s)-grids%YgridExtrap(typeSim,t,:)))
                gridY =grids%Ygrid(typeSim,t,(/(i,i=1,numPointsY-1,2)/))
                yval = yex(t,s)
                y(t, s) =  yex(t,s)
            end if


            idxA=minloc(abs(a(t, s)-grids%Agrid(typeSim,t,:)))
            idxAIME=minloc(abs(AIME(t, s)-grids%AIMEgrid(typeSim,t,:)))
            l(t,s)=locpolicyL(t,typeSim,idxA(1),idxAIME(1),idxY(1),unemployed)

            y(t, s) = yex(t,s)
            !if (maxval(Lcube)==minval(Lcube)) then
            call linearinterp3_withextrap(grids%Agrid(typeSim,t,:), grids%AIMEgrid(typeSim,t,:),gridY , &
                numPointsA, numAIME, numPointsProd, a(t, s), AIME(t, s),yval, ltemp,  real(locpolicyL(t, typeSim, :, :,:,unemployed),rk))
            if (abs(l(t,s) - ltemp) > 0.01) then
                lbA1 = grids%Agrid(typeSim,t + 1, 1);          ! lower bound: assets tomorrow
                ev1 = EV(t + 1,typeSim,:,:,idxY(1))
                !We can only end up here for ages below 80 as for >80 all entries in polL are 0
                !This also takes account of uneployment shock as the wage offer they have is 0 so won't work
                call solvePeriod(params, grids, y(t, s), a(t, s), AIME(t, s) ,t, typeSim, lbA1, ev1, &
                    grids%benefit(t),a(t+1, s),c(t, s),l(t,s),v(t  , s))
            else
                call linearinterp3_withextrap(grids%Agrid(typeSim,t,:), grids%AIMEgrid(typeSim,t,:), gridY, &
                    numPointsA, numAIME, numPointsProd, a(t, s), AIME(t, s), yval, a(t+1, s),  locpolicyA1(t,typeSim, :, :,:,unemployed))
                call linearinterp3_withextrap(grids%Agrid(typeSim,t,:), grids%AIMEgrid(typeSim,t,:), gridY,&
                    numPointsA, numAIME, numPointsProd+1,  a(t, s), AIME(t, s), yval, v(t  , s),  locEV(t,typeSim, :, :,:,unemployed))
            end if

            if (l(t,s) .EQ. 0) then
                y(t, s)=grids%benefit(t)
                !AIME(t+1, s) =   AIME(t, s) * (t-1)/t
            end if
            workAge = startAge - 20 + t
            if  (t < Tretire) then ! spouseretire
                if (l(t,s) .EQ. 0) then
                    AIME(t+1, s) =   AIME(t, s) * (workAge-1)/workAge
                else
                    AIME(t+1, s) =  y(t,s)/workAge + AIME(t, s) * (workAge-1)/workAge
                end if
            else
                AIME(t+1, s) =   AIME(t, s)
            end if

            call gross2net(params,y(t, s),t,l(t,s),typeSim, AIME(t, s))
            ! tomorrow's optimal assets
            !! THIS is WRONG BECAUSE ABOVE IS WRONG!
            c(t, s) = a(t, s)  + y(t, s) - (a(t+1, s)/(1+params%r))

        end do   !t
    end do! s

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    !!truncate obsrvtion
    !---------------------------------------------------------------------------------------------------------!
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
    !!Convert Groos to net income
    !---------------------------------------------------------------------------------------------------------!
    subroutine gross2net(params,Y,ixt,ixl,ixType, AIME)
    implicit none

    !inputs
    type (structparamstype), intent(in) :: params
    integer, intent(in) :: ixt,ixl, ixType
    real (kind=rk), intent(in) :: AIME

    !changing
    real (kind =rk), intent(inout) :: Y

    !Add own and husbands pension
    if (ixt >= Tretire) then
        Y = Y + params%pension
        if  (ixt >= spouseretire) then
            Y = Y + params%pension
            if (ixL==0 .AND. mod(ixType,2)==0) then ! 
                Y =   Y + dbPension(params,AIME);
            end if
        else
            Y = Y +params%spouseInc(ixType)
        end if
    else
        Y = Y +params%spouseInc(ixType)
    end if

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!Returns GMM Cretieria
    function gmm(params,grids,target,weights)
    implicit none

#ifdef mpi
    include 'mpif.h'
#endif  

    !inputs
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(inout) :: grids
    real (kind=rk), intent(in) :: target(:,:), weights(:,:)

    !output
    real (kind=rk) :: gmm

    !local
    real (kind=rk) :: V(numPointsType,Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: policyA1(numPointsType,Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: policyC(numPointsType,Tperiods, numPointsA, numPointsY, numAIME)
    integer :: policyL(numPointsType,Tperiods, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: EV(numPointsType,Tperiods+1, numPointsA, numPointsY, numAIME)
    real (kind=rk) :: ypath(Tperiods, numSims) !income
    real (kind=rk) :: cpath(Tperiods, numSims)  !consumption
    integer :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk) :: vpath(Tperiods, numSims)  !value
    real (kind=rk) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk) ::  meanL(Tperiods), meanA(Tperiods)
    real (kind=rk) :: AIME(Tperiods + 1,numSims)

    integer :: n, typeSim

    !Set asset grid
    do typeSim = 1, numPointsType
        call getassetgrid( params, grids%maxInc(typeSim,:), grids%Agrid(typeSim,:,:))
    end do
    !solve
    call solveValueFunction( params, grids, policyA1, policyC, policyL, V, EV, .FALSE. )
    if (rank==0) then
        !simulate
        call simWithUncer(params, grids, policyA1,policyL,EV, grids%Simy, cpath, apath, vpath, lpath, ypath, AIME )
        do n=1,Tperiods
            meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
            meanA(n)=sum(real(Apath(n,:),rk))/real(numSims,rk) !size(lpath(n,:))
        end do
        if (fullLifeCycle) then
            gmm = dot_product(abs(meanL(32:32+23)-target(1,:)),abs(meanL(32:32+23)-target(1,:)))! + &
                !dot_product(abs(meanA(32:32+23)-target(2,:)),abs(meanA(32:32+23)-target(2,:)))
            else
        gmm = dot_product(abs(meanL(33 - (StartAge-20):33 - (StartAge-20)+23)-target(1,:)),weights(1,:)*abs(meanL(33 - (StartAge-20):33 - (StartAge-20)+23)-target(1,:))) + &
            dot_product(abs(meanA(33 - (StartAge-20):33 - (StartAge-20)+23)-target(2,:)),weights(2,:)*abs(meanA(33 - (StartAge-20):33 - (StartAge-20)+23)-target(2,:)))
        end if
    end if

#ifdef mpi 
    call MPI_Bcast( gmm,   1,   mpi_double_precision, 0, mpi_comm_world, ierror)
    call mpi_barrier(mpi_comm_world, ierror)
    if (ierror.ne.0) stop 'mpi problem180'
#endif

    end function
    ! ---------------------------------------------------------------------------------------------------------!
    ! !Write to file
    !---------------------------------------------------------------------------------------------------------!
    subroutine writetofile(ypath, cpath, apath, vpath, lpath, yemp, AIME)
    implicit none

    !inputs
    !type (structparamstype), intent(in) :: params
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
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), meanPoor(Tperiods), meanRich(Tperiods), numLC(Tperiods)
    integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    integer :: rich, poor

#ifdef win   
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
        numLC(n) = count(apath(n,:).eq.0)
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

    inquire (iolength=requiredl)  numLC
    open (unit=207, file='..\\out\numLC.txt', status='unknown',recl=requiredl, action='write')
    write (207, '(6E15.3)' ) numLC
    close( unit=207)    
    
#else
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
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk)
        write (201, * ) meanL(n)
        meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
        write (202, * ) meanV(n)
        meanA(n)=sum(real(apath(n,:),rk))/real(numSims,rk)
        write (203, * ) meanA(n)
        meanC(n)=sum(real(cpath(n,:),rk))/real(numSims,rk)
        write (204, * ) meanC(n)
        meanY(n)=sum(real(ypath(n,:),rk))/real(numSims,rk)
        write (205, * ) meanY(n)
        meanYemp(n)=sum(real(yemp(n,:),rk))/real(numSims,rk)
        write (206, * ) meanYemp(n)
        meanAIME(n)=sum(real(AIME(n,:),rk))/real(numSims,rk)
        write (207, * ) meanAIME(n)
        meanRich(n)=sum(real(lrich(n,:),rk))/real(numSims/2,rk)
        write (209, * ) meanRich(n)
        meanPoor(n)=sum(real(lpoor(n,:),rk))/real(numSims/2,rk)
        write (208, * ) meanPoor(n)
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
#endif

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! !Write to file
    !---------------------------------------------------------------------------------------------------------!
    subroutine writetofileByType(ypath, cpath, apath, vpath, lpath, yemp, AIME)
    implicit none

    !inputs
    !type (structparamstype), intent(in) :: params
    real (kind=rk), intent(in) :: ypath(Tperiods, numSims) !income
    real (kind=rk), intent(in) :: cpath(Tperiods, numSims)  !consumption
    integer, intent(in) :: lpath(Tperiods, numSims) !labour supply
    real (kind=rk), intent(in) :: vpath(Tperiods, numSims)  !value
    real (kind=rk), intent(in) :: apath(Tperiods + 1,numSims) !this is the path at the start of each period, so we include the 'start' of death
    real (kind=rk), intent(in)  :: yemp(Tperiods, numSims)
    real (kind=rk), intent(in)  :: AIME(Tperiods + 1,numSims)

    !local
    integer :: n, requiredl , i, finish, start
    integer(1) :: k
    real(kind=rk) :: meanL(Tperiods), meanV(Tperiods), meanA(Tperiods), meanC(Tperiods), meanY(Tperiods), medianA
    real(kind=rk) :: meanYemp(Tperiods), meanAIME(Tperiods), meanPoor(Tperiods), meanRich(Tperiods), numLC(Tperiods)
    integer :: lrich(Tperiods,numSims/2), lpoor(Tperiods,numSims/2)
    integer :: rich, poor
    character(len=1024) :: outFile
    integer, allocatable :: rangeSims(:)

#ifdef win
    do k=1,numPointsType
        if (k == 1)  then
            start = 1
            finish = 0.16*numSims
            allocate(rangeSims(finish))
            rangeSims = (/(i,i=1,finish)/)
        else if (k==2) then
            deallocate(rangeSims)
            start = finish +1
            finish = 0.48*numSims
            allocate(rangeSims(finish-start+1))
            rangeSims = (/(i,i=start,finish)/)
            !real(s,rk)/numSims <0.48
        else if (k==3) then
            deallocate(rangeSims)
            start = finish +1
            finish = 0.59*numSims
            allocate(rangeSims(finish-start+1))
            rangeSims = (/(i,i=start,finish)/)
            !real(s,rk)/numSims <0.59
        else
            deallocate(rangeSims)
            start = finish +1
            finish = numSims
            allocate(rangeSims(finish-start+1))
            rangeSims = (/(i,i=start,finish)/)
        end if

        do n=1,Tperiods
            meanL(n)=sum(real(lpath(n,rangeSims),rk))/real(finish-start+1,rk) !size(lpath(n,:))
            !write (201, * ) meanL(n)
            meanV(n)=sum(real(vpath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (202, * ) meanV(n)
            meanA(n)=sum(real(apath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (203, * ) meanA(n)
            meanC(n)=sum(real(cpath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (204, * ) meanC(n)
            meanY(n)=sum(real(ypath(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (205, * ) meanY(n)
            meanYemp(n)=sum(real(yemp(n,rangeSims),rk))/real(finish-start+1,rk)
            !write (205, * ) meanY(n)
            meanAIME(n)=sum(real(AIME(n,rangeSims),rk))/real(finish-start+1,rk)
            numLC(n) = count(apath(n,:).eq.0)
        end do

        !L
        write (outFile, *), trim("..\\out\lpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanL
        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
        write (201, '(6E15.3)' ) meanL
        close( unit=201)

        !V
        write (outFile, *), trim("..\\out\Vpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanV
        open (unit=202, file=outFile, status='unknown',recl=requiredl, action='write')
        write (202, '(6E15.3)' ) meanV
        close( unit=202)

        !A
        write (outFile, *), trim("..\\out\Apath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanA
        open (unit=203, file=outFile, status='unknown',recl=requiredl, action='write')
        write (203, '(6E15.3)' ) meanA
        close( unit=203)

        !C
        write (outFile, *), trim("..\\out\Cpath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanC
        open (unit=204, file=outFile, status='unknown',recl=requiredl, action='write')
        write (204, '(6E15.3)' ) meanC
        close( unit=204)

        !Y
        write (outFile, *), trim("..\\out\Ypath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanY
        open (unit=205, file=outFile, status='unknown',recl=requiredl, action='write')
        write (205, '(6E15.3)' ) meanY
        close( unit=205)

        write (outFile, *), trim("..\\out\YempPath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanYemp
        open (unit=206, file=outFile, status='unknown',recl=requiredl, action='write')
        write (206, '(6E15.3)' ) meanYemp
        close( unit=206)

        write (outFile, *), trim("..\\out\AIMEPath"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  meanAIME
        open (unit=207, file=outFile, status='unknown',recl=requiredl, action='write')
        write (207, '(6E15.3)' ) meanAIME
        close( unit=207)
        
        write (outFile, *), trim("..\\out\numLC"),k,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  numLC
        open (unit=207, file=outFile, status='unknown',recl=requiredl, action='write')
        write (207, '(6E15.3)' ) numLC
        close( unit=207)        
    end do
#else
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
        meanL(n)=sum(real(lpath(n,:),rk))/real(numSims,rk)
        write (201, * ) meanL(n)
        meanV(n)=sum(real(vpath(n,:),rk))/real(numSims,rk)
        write (202, * ) meanV(n)
        meanA(n)=sum(real(apath(n,:),rk))/real(numSims,rk)
        write (203, * ) meanA(n)
        meanC(n)=sum(real(cpath(n,:),rk))/real(numSims,rk)
        write (204, * ) meanC(n)
        meanY(n)=sum(real(ypath(n,:),rk))/real(numSims,rk)
        write (205, * ) meanY(n)
        meanYemp(n)=sum(real(yemp(n,:),rk))/real(numSims,rk)
        write (206, * ) meanYemp(n)
        meanAIME(n)=sum(real(AIME(n,:),rk))/real(numSims,rk)
        write (207, * ) meanAIME(n)
        meanRich(n)=sum(real(lrich(n,:),rk))/real(numSims/2,rk)
        write (209, * ) meanRich(n)
        meanPoor(n)=sum(real(lpoor(n,:),rk))/real(numSims/2,rk)
        write (208, * ) meanPoor(n)
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
#endif

    end subroutine

    ! ---------------------------------------------------------------------------------------------------------!
    ! get asset grid
    !---------------------------------------------------------------------------------------------------------!
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
        span = (log(1.0+maxA(ixt))-log(1.0))/ (numPointsA-1)
        loggrid = log(1.0)+span*(/(i,i=0,numPointsA-1)/)
        Agrid(ixt, :) = (/((exp(loggrid(i))-1.0),i=1,numPointsA)/)
        !span =  (log(1.0+log(1.0+log(1+maxA(ixt)))) - log(1.0+log(1.0+log(1.0))) )/ (numPointsA-1)
        !loggrid = log(1.0+log(1.0+log(1.0))) + span*(/(i,i=0,numPointsA-1)/)
        !Agrid(ixt, :) = (/(exp(exp(exp(loggrid(i))-1.0)-1.0)-1.0,i=1,numPointsA)/) !exp(exp(exp(loggrid)-1)-1)-1
    end do
    test = sum(Agrid(1, :))/size(Agrid(ixt, :))

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !!solve period
    subroutine solvePeriod(params, grids, Yin, A, AIMEin ,ixt,ixType, lbA1, EV1, benefit,policyA1,policyC,policyL,V)
    implicit none

    !input
    real(kind=rk), intent(in) :: Yin, A, lba1, EV1(:,:), benefit, AIMEin
    type (structparamstype), intent(in) :: params
    type (gridsType), intent(in) :: grids
    integer, intent(in) :: ixt, ixType

    !output
    real(kind=rk), intent(out) :: policyA1, policyC, V
    integer, intent(out) :: policyL

    !local
    integer :: ixl, workAge
    real(kind=rk) :: Y, negV, negVtemp, ubA1, policyA1temp, AIME

    negV = -huge(negv) !-log(0.0) !inf
    AIME = AIMEin
    workAge = startAge - 20 + ixt
    do ixL = 0,(numPointsL-1),1           ! points of income choice
        ! Value of income and information for optimisation
        Y    = ixL*Yin+ (1-ixL)*benefit;

        !AIME only depends on earned income so add spousal
        !and pensions after calculating it
        call gross2net(params,Y,ixt,ixl,ixType,AIME)

        ! Next peridos AIME
        if (ixL==1 .AND.  ixt < Tretire ) then ! spouseretire
            AIME =  Yin/workAge + AIME * (workAge-1)/workAge
        else if ((ixL==0 .AND.  ixt < Tretire )) then !spouseretire
            AIME = AIME * (workAge-1)/workAge
        end if

        ubA1 = (A + Y - params%minCons)*(1+params%r);    ! upper bound: assets tomorrow

        ! Compute solution
        if (ubA1 - lbA1 < params%minCons) then         ! if liquidity constrained
            negVtemp = objectivefunc(params, grids,lbA1, A, Y,ixL,ixt,ixType, AIME,EV1);
            policyA1temp = lbA1;
        else                               ! if interior solution
            negVtemp = golden_generic(lbA1, ubA1, policyA1temp, func,params%tol,.FALSE.)
        end if

        if (negVtemp > negV) then
            negV = negVtemp
            policyA1=policyA1temp
            policyL=ixL
            policyC = A + Y - policyA1/(1+params%r)
            ! Store solution and its value
        end if
    end do

    V  = negV;
    contains
    function func(x)
    real (kind = rk), intent(in) :: x
    real (kind = rk) :: func
    func = objectivefunc(params, grids,x, A, Y,ixL,ixt,ixType, AIME,EV1)
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
    subroutine initialGuess(rank,params,grids,moments,weights,p,y)
    implicit none

    !inputs
    integer, intent(in) :: rank
    real(kind=rk), intent(in) :: moments(2,24), weights(:,:)

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

        if (rank==0) then
            print '("Guess 2")'
        end if
        p(2,1) = 0.4637
        p(2,2) = 0.970
        P(2,3) = 1
        p(2,4) = params%db(1)*0.9
        p(2,5) = params%db(2)*1.2
        y(2) = gmm_criteria(p(2,:))

        if (rank==0) then
            print '("Guess 3")'
        end if
        p(3,1) = 0.322
        p(3,2) = 0.9843
        P(3,3) = 2
        p(3,4) = params%db(1)*1.1
        p(3,5) = params%db(2)*0.7
        y(3) = gmm_criteria(p(3,:))

        if (rank==0) then
            print '("Guess 4")'
        end if
        p(4,1) = 0.55
        p(4,2) = 0.96
        P(4,3) = 0.5
        p(4,4) = params%db(1)*1.3
        p(4,5) = params%db(2)*0.95
        y(4) = gmm_criteria(p(4,:))

        if (rank==0) then
            print '("Guess 5")'
        end if
        p(5,1) = 0.15
        p(5,2) = 0.9999
        P(5,3) = 4
        p(5,4) = params%db(1)*0.85
        p(5,5) = params%db(2)*1.15
        y(5) = gmm_criteria(p(5,:))

        if (rank==0) then
            print '("Guess 6")'
        end if
        p(6,1) = 0.27
        p(6,2) = 0.986
        P(6,3) = 0.9
        p(6,4) = params%db(1)*0.99
        p(6,5) = params%db(2)*0.87
        y(6) = gmm_criteria(p(6,:))

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
            p(i,6) = params%thetab*(1+uniformRand(i))
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
    params%thetab = control(6)
    gmm_criteria = gmm(params,grids,moments,weights) !*-1.0

    end function
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Add Unemployment shocks
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine addUnemploymentshock(Ygrid, incTransitionMrx,grids,typeSim)
    implicit none

    !Input
    real (kind=rk), intent(in) :: Ygrid(:,:), incTransitionMrx(:,:)
    integer, intent(in) :: typeSim

    !Changingin
    type (gridsType), intent(inout) :: grids

    !local
    real (kind=rk) :: reemployRaw(4,10), unemployRaw(4,10)
    integer :: i

    unemployRaw(1,:) = (/0.227, 0.049, 0.027, 0.014, 0.004, 0.007, 0.013, 0.017, 0.007, 0.01/)
    unemployRaw(2,:) = (/0.099, 0.026, 0.01, 0.013, 0.008, 0.012, 0.003, 0.003, 0.002, 0.0/)
    unemployRaw(3,:) = (/0.201, 0.024, 0.054, 0.032, 0.029, 0.038, 0.02, 0.005, 0.004, 0.01/)
    unemployRaw(4,:) = (/0.063, 0.016, 0.01, 0.006, 0.003, 0.002, 0.004, 0.001, 0.003, 0.0/)

    if (numPointsProd==5) then
        do i=1,numPointsType
            grids%unemploy(i,1) =unemployRaw(i,1)+unemployRaw(i,2)!(/0.1762, 0.0415,0.0257,0.0165,0.0090/)
            grids%unemploy(i,2) =unemployRaw(i,3)+unemployRaw(i,4)
            grids%unemploy(i,3) =unemployRaw(i,5)+unemployRaw(i,6)
            grids%unemploy(i,4) =unemployRaw(i,7)+unemployRaw(i,8)
            grids%unemploy(i,5) =unemployRaw(i,9)+unemployRaw(i,10)
        end do
    else if (numPointsProd==10) then
        grids%unemploy = unemployRaw
        !grids%unemploy =(/0.1475, 0.02875, 0.02525, 0.01625, 0.011, 0.01475, 0.01, 0.0065, 0.004, 0.005/)
        !else if (numPointsProd==20) then
        !!grids%unemploy =(/0.1475,0.1475, 0.02875, 0.02875, 0.02525,0.02525, 0.01625, 0.01625, 0.011, 0.011, 0.01475, 0.01475, 0.01, 0.01, 0.0065, 0.0065, 0.004, 0.004, 0.005, 0.005/)
    end if

    reemployRaw(1,:) = (/0.13, 0.043, 0.022, 0.005, 0.002, 0.007, 0.01, 0.002, 0.0, 0.007 /)
    reemployRaw(2,:) = (/0.174, 0.005, 0.013, 0.017, 0.013, 0.017, 0.013, 0.003, 0.003, 0.0 /)
    reemployRaw(3,:) = (/0.118, 0.029, 0.037, 0.011, 0.008, 0.013, 0.003, 0.0, 0.005, 0.0 /)
    reemployRaw(4,:) = (/0.238, 0.032, 0.02, 0.008, 0.008, 0.012, 0.008, 0.004, 0.0, 0.004/)


    if (numPointsProd==5) then
        !grids%reemploy = (/ 0.1923, 0.0333, 0.0200, 0.0108,0.0047/)
        do i=1,numPointsType
            grids%reemploy(i,1) =reemployRaw(i,1)+reemployRaw(i,2)!(/0.1762, 0.0415,0.0257,0.0165,0.0090/)
            grids%reemploy(i,2) =reemployRaw(i,3)+reemployRaw(i,4)
            grids%reemploy(i,3) =reemployRaw(i,5)+reemployRaw(i,6)
            grids%reemploy(i,4) =reemployRaw(i,7)+reemployRaw(i,8)
            grids%reemploy(i,5) =reemployRaw(i,9)+reemployRaw(i,10)
        end do
    else if (numPointsProd==10) then
        grids%reemploy = reemployRaw
        !grids%reemploy =(/0.165, 0.02725, 0.023, 0.01025, 0.00775, 0.01225, 0.0085, 0.00225, 0.002, 0.00275/)
        !else if (numPointsProd==20) then
        !grids%reemploy =(/0.165, 0.165, 0.02725, 0.02725, 0.023, 0.023, 0.01025, 0.01025, 0.00775, 0.00775, 0.01225, 0.01225, 0.0085, 0.0085, 0.00225, 0.00225, 0.002, 0.002, 0.00275, 0.00275/)
    end if

    grids%Ygrid(typeSim,:,:) = 0.0
    grids%incTransitionMrx(:,:) = 0.0
    grids%Ygrid(typeSim,:,(/(i,i=1,numPointsY-1,2)/))= Ygrid(:,:)
    do i = 1, numPointsY
        if (mod(i,2)==1) then
            grids%incTransitionMrx(i,(/(i,i=1,numPointsY-1,2)/)) = (1-grids%unemploy(typeSim,i/2+1))*incTransitionMrx(i/2+1,:)
            grids%incTransitionMrx(i,(/(i,i=2,numPointsY,2)/)) = grids%unemploy(typeSim,i/2+1)*incTransitionMrx(i/2+1,:)
        else
            grids%incTransitionMrx(i,i-1) = grids%reemploy(typeSim,i/2)
            grids%incTransitionMrx(i,i) = 1-grids%reemploy(typeSim,i/2)
        end if
    end do

    end subroutine
    !!-----------------------------------------------------------------------------------------------------------!
    !Intialise gues for Nedler-Mead Algorithm
    !!-----------------------------------------------------------------------------------------------------------!
    subroutine setupMisc(params,grids)
    implicit none

    !Changingin
    type (structparamstype), intent(inout) :: params
    type (gridsType), intent(inout) :: grids

    !local
    REAL(kind=rk) :: temp(85),  temp2(86), temp3(numSims, Tperiods)
    integer :: ios, I
    real (kind=rk) ::  sig_inc, sig_initial
    real (kind=rk) :: e(Tperiods, numSims)   ! the innovations to log income
    real (kind=rk) :: logy1(1,numSims)        ! draws for initial income
    real (kind=rk) :: ly(Tperiods, numSims)           ! log income
    integer :: s, t,  workAge
    real (kind=rk) ::  uniformRand(Tperiods,numSims)
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    integer :: seedIn, n, productiv(1), typeSim, requiredl, finish
    logical :: unemployed
    character(len=1024) :: outFile
    integer, allocatable :: rangeSims(:)

#ifdef win
    open (unit = 1001,file='..\\..\\moments\\InitialAssets.txt', action='read', IOSTAT = ios)
#else
    open (unit = 1001,file='../moments/InitialAssets.txt', action='read', IOSTAT = ios)
#endif
    read (1001, *,  IOSTAT = ios) grids%initialAssets
    close (unit=1001)

    temp = (/-0.01374, -0.01069, -0.00767, -0.00467, -0.00169, 0.00126, 0.00418,  0.00708, 0.00996, 0.01281, 0.01563, 0.01843, &
        0.02121, 0.02396, 0.02668,  0.02938,  0.03206, 0.03471, 0.03733,  0.03994, 0.04251, 0.04506, 0.04759, 0.05009, &
        0.05256, 0.05501,  0.05744,  0.05984,  0.06221, 0.06457, 0.06689,  0.06919, 0.07147, 0.07372, 0.07594, 0.07814, &
        0.08032, 0.08247,  0.08460,  0.08670,  0.08877, 0.09082, 0.09285,  0.09485, 0.09683, 0.09878, 0.10070, 0.10261, &
        0.10448, 0.10633, 0.10816,  0.10996,  0.11174, 0.11349, 0.11521,  0.11691, 0.11859, 0.12024, 0.12187, 0.12347, &
        0.12505, (I*0.0,I=1,24)/)
    grids%fc = temp(startAge-20+1:85)

    call getIncomeGrid(params, grids)

    temp2 = (/0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.011071, 0.011907, 0.012807, &
        0.013676, 0.01475, 0.015818, 0.016846, 0.018079, 0.019343, 0.020659, 0.0218, 0.023505, 0.025202, 0.02696, &
        0.028831, 0.031017, 0.033496, 0.036024, 0.038912, 0.042054, 0.045689, 0.049653, 0.054036, 0.05886, 0.064093, &
        0.069636, 0.07533, 0.081069, 0.086912, 0.093067, 0.099807, 0.107483, 0.116125, 0.125196, 0.134361, 0.143881, &
        0.1542, 0.165675, 0.17842, 0.192363, 0.2117, 0.1672,0.1565, 0.1485,0.1459, 1.0/);
    grids%mortal = temp2(startAge-20+1:86)

    ! Obtain time series of incomes for our simulated individuals
    ! Draw random draws for starting income and for innovations
    sig_inc = params%sigma/ ((1-params%rho**2)**0.5)
    sig_initial = 0.2950 !mean([0.073 0.053 0.110 0.112])^0.5;
    ! normally distributed random draws for the innovation
#ifdef win
    open (unit = 1001,file='..\\Errors.txt', action='read', IOSTAT = ios)
#else
    open (unit = 1001,file='./data/Errors.txt', action='read', IOSTAT = ios)
#endif
    read (1001, *) temp3
    close (unit=1001)
    e = transpose(temp3)
#ifdef win
    open (unit = 1001,file='..\\IntProd.txt', action='read', IOSTAT = ios)
#else
    open (unit = 1001,file='./data/IntProd.txt', action='read', IOSTAT = ios)
#endif

    read (1001, *) logy1
    close (unit=1001)

    seedIn = 01042017
    !Set seed
    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))
    seed = seedIn * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)
    DEALLOCATE(seed)

    !!get uniform random number
    CALL RANDOM_NUMBER(uniformRand)

    do s = 1, numSims, 1                           ! loop through individuals
        if (numPointsType == 1) then
            typeSim = 1
        else
            if (real(s,rk)/numSims <0.16) then
                typeSim = 1
            else if (real(s,rk)/numSims <0.48) then
                typeSim = 2
            else if (real(s,rk)/numSims <0.59) then
                typeSim =3
            else
                typeSim =4
            end if
        end if
        unemployed = .FALSE.
        ! Get all the incomes, recursively
        ly(1, s) = truncate(logy1(1, s), -normBnd*sig_inc,normBnd*sig_inc )
        grids%Simy(1, s) = exp(ly(1, s)+params%delta(typeSim,1)*(startAge - 20 + 1)**2+params%delta(typeSim,2)*(startAge - 20 + 1)+params%delta(typeSim,3)-grids%fc(1))
        do t = 1,Tperiods-1,1                              ! loop through time periods for a particular individual
            workAge = startAge - 20 + t
            if (unemployed) then
                if (uniformRand(t,s) < grids%reemploy(typeSim,productiv(1))) then
                    unemployed = .FALSE.
                    ly(t+1, s) = (1 -params%rho) * params%mu + params%rho * ly(t, s) + e(t + 1, s)
                    ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
                    grids%Simy(t+1, s) = exp( ly(t+1, s) + params%delta(typeSim,1)*(workAge+1)**2+params%delta(typeSim,2)*(workAge+1)+params%delta(typeSim,3)-grids%fc(t) )
                else
                    ly(t+1, s)  = ly(t, s)
                    grids%Simy(t+1, s) = 0
                end if
            else
                productiv = minloc((grids%ygrid(typeSim,t,(/(i,i=1,numPointsY-1,2)/))-grids%Simy(t, s)),(grids%ygrid(typeSim,t,(/(i,i=1,numPointsY-1,2)/))-grids%Simy(t, s))>0)
                if (uniformRand(t,s) < grids%unemploy(typeSim,productiv(1))) then
                    unemployed = .TRUE.
                    ly(t+1, s)  = ly(t, s)
                    grids%Simy(t+1, s) = 0
                else
                    ly(t+1, s) = (1 -params%rho) * params%mu + params%rho * ly(t, s) + e(t + 1, s)
                    ly(t+1, s) = truncate(ly(t+1, s), -normBnd*sig_inc,normBnd*sig_inc )
                    grids%Simy(t+1, s) = exp( ly(t+1, s) + params%delta(typeSim,1)*(workAge+1)**2+params%delta(typeSim,2)*(workAge+1)+params%delta(typeSim,3)-grids%fc(t) )
                end if
            end if
        end do ! t

    end do! s

    if (rank==0) then
        !write (*,*) rank
        finish = 0.16*numSims
        allocate(rangeSims(finish))
        rangeSims = (/(i,i=1,finish)/)
        write (outFile, *), trim("..\\temp\ypath"),1,".txt"
        outfile=ADJUSTL(outfile)
        inquire (iolength=requiredl)  grids%Simy(:, rangeSims)
        open (unit=201, file=outfile, status='unknown',recl=requiredl, action='write')
        write (201, '(6E15.3)' ) grids%Simy(:, rangeSims)
        close( unit=201)
    end if

    end subroutine
    end module routines
