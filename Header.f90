    module Header

    integer, parameter :: rk = selected_real_kind(15)

    integer, parameter :: numPointsA = 20 !45
    integer, parameter :: numPointsY = 20
    integer, parameter :: numAIME = 5 !5
    integer, parameter :: numPointsL = 2
    integer, parameter :: numSims = 10000 !10000
    integer, parameter :: Tperiods = 75
    integer, parameter :: Tretire =40
    integer, parameter :: normBnd = 4
    integer, parameter :: dimEstimation = 3
    integer, parameter :: spouseretire = 45
    integer, parameter :: stopwrok = 60

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
    end type gridsType

    end module Header
