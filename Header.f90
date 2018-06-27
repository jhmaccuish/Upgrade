    module Header

    integer, parameter :: rk = selected_real_kind(15)

    integer, parameter :: numPointsA = 30 !30 !20 !30 !45
    integer, parameter :: numPointsProd = 10
    integer, parameter :: numPointsY = 2*numPointsProd !20
    integer, parameter :: numAIME = 30 !10 !5
    integer, parameter :: numPointsL = 2
    integer, parameter :: numSims = 1016!  250! 5000! 10000!
    integer, parameter :: startAge =  52 !20! 
    integer, parameter :: endAge = 105 
    integer, parameter :: Tperiods = endAge -startAge
    integer, parameter :: Tretire =62 -startAge
    integer, parameter :: normBnd = 4
    integer, parameter :: dimEstimation = 5
    integer, parameter :: spouseretire = 65 -startAge
    integer, parameter :: stopwrok = 80 -startAge

    !Holds the structural parameters that will eventually be estimated (at least in some cases)
    type structparamstype
        !Personal
        real (kind=rk) :: gamma
        real (kind=rk) :: r
        real (kind=rk) :: beta
        real (kind=rk) :: sigma
        real (kind=rk) :: mu
        real (kind=rk) :: rho
        real (kind=rk) :: nu
        !Instutional
        real (kind=rk) :: delta(3)
        real (kind=rk) :: pension
        real (kind=rk) :: hrsWrk
        real (kind=rk) :: spouseInc
        real (kind=rk) :: minCons
        real (kind=rk) :: db(2)
        real (kind=rk) :: startA
        real (kind=rk) :: thetab
        real (kind=rk) :: k
        !real (kind=rk) :: tol, minCons
        real (kind=rk) :: tol
        integer :: system
    end type structparamstype

    type gridsType
        real (kind=rk) :: Agrid(Tperiods+1,numPointsA)
        real (kind=rk) :: Ygrid(Tperiods,numPointsY)
        real (kind=rk) :: incTransitionMrx(numPointsY,numPointsY)
        real (kind=rk) :: AIMEgrid(Tperiods+1,numAIME)
        real (kind=rk) :: benefit(Tperiods)
        real (kind=rk) :: fc(Tperiods)
        real (kind=rk) :: maxInc(Tperiods)
        real (kind=rk) :: initialAssets(1016)
    end type gridsType

    !! For mpi
    integer :: rank, ierror, procsize
    integer, parameter :: mpiDim =  numPointsA * numAIME  

    !! Test controls
    logical, parameter :: fullLifeCycle = .FALSE.
    end module Header
