
    ! This contains generic routines
    ! Random number routines here are scheduled for withdrawal. Need to replace with new ones
    module routines_generic
    !use precisions
    use Header

    contains
    subroutine linearinterp1(xgrid, ygrid, length, xval, yval, extrapbot, extraptop)
    ! Subroutine for linear interpolation
    ! Arguments are 1. (input) xvlalues
    ! 2. (input) yvalues
    ! 3. (input) length of xgrid and ygrid,
    ! 4. (input) value of x to interpolate at
    ! 5. (output) interpolated value
    ! 6. (input) indicator for whether to extrapolate below bottom point of grid
    ! 7. (input indicator for whether to extrapolate above top point of grid

    ! If the choose not to extrapolate interpolated value is set equal to either the top or the bottom
    ! point of the grid
    ! It is fully generic

    implicit none
    ! Arguments
    integer, intent (in) :: length
    integer, intent (in) :: extrapbot
    integer, intent (in) :: extraptop
    real (kind=rk), intent (in) :: xgrid(length)
    real (kind=rk), intent (in) :: ygrid(length)
    real (kind=rk), intent (in) :: xval
    real (kind=rk), intent (out) :: yval

    ! For programme
    integer :: xlowerloc
    integer scratch_xnearestloc ! This areguments required by findingridr but I'm not using them here - hence
    ! scratch

    intrinsic size

    ! Checking inputs - have commented this out for speed
    !   if ((size(xgrid).ne.length) .or. (size(ygrid).ne.length)) then
    !     write (*, *) 'An argument in linearinterp1 is not of the specified length'
    !     stop
    !   end if

    !   if ((extrapbot.ne.0 .and. extrapbot.ne.1) .or. (extraptop.ne.0 .and. extraptop.ne.1)) then
    !     write (*, *) 'An argument in linearterp1 that needs to be either 0 or 1 is not so'
    !     stop
    !   end if

    if (abs(xgrid(length)-xgrid(1)).lt.0.1d-15) then ! If the top of the grid is equal to the bottom
        yval = ygrid(1) ! (which would include the scenario that all the entries are the same) then I
        ! assign the bottom value of y to yval
    else

        call findingridr(xgrid, length, xval, scratch_xnearestloc, xlowerloc)
        call linearinterpfromlocations(xgrid, ygrid, xlowerloc, xval, length, extrapbot, extraptop, yval)
    end if

    end subroutine linearinterp1

    ! ---------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------!


    ! Extrapolation is by default here. There is no option to use the bottom value (as I have
    ! in the one-dimensional linear interpolation).

    subroutine linearinterp2_withextrap(x1grid, x2grid, dim1, dim2, x1val, x2val, yval, extrapbot, extraptop, ygrid)
    implicit none
    ! Arguments
    integer, intent (in) :: dim1
    integer, intent (in) :: dim2
    real (kind=rk), intent (in) :: x1grid(dim1)
    real (kind=rk), intent (in) :: x2grid(dim2)
    real (kind=rk), intent (in) :: ygrid(dim1, dim2)
    real (kind=rk), intent (in) :: x1val
    real (kind=rk), intent (in) :: x2val
    real (kind=rk), intent (out) :: yval
    integer, intent (in) :: extrapbot
    integer, intent (in) :: extraptop

    ! For programme
    integer :: x1lowerloc
    integer :: x2lowerloc
    integer :: scratch_nearestloc
    integer :: newxgridlowerloc
    real (kind=rk) :: newygrid(2)
    real (kind=rk) :: newxgrid(2)

    if ((size(x1grid).ne.dim1) .or. (size(x2grid).ne.dim2)) then
        write (*, *) 'An argument in linearinterp2 is not of the specified length'
        stop
    end if

    ! Find the square
    call findingridr(x1grid, dim1, x1val, scratch_nearestloc, x1lowerloc)
    call findingridr(x2grid, dim2, x2val, scratch_nearestloc, x2lowerloc)

    ! Find if any of the values are off the grid
    if (x1lowerloc.eq.0) then
        x1lowerloc = 1
    end if

    if (x2lowerloc.eq.0) then
        x2lowerloc = 1
    end if

    if (x1lowerloc.eq.dim1) then
        x1lowerloc = dim1 - 1
    end if

    if (x2lowerloc.eq.dim2) then
        x2lowerloc = dim2 - 1
    end if

    ! Now do linearinterpolation along the first dimension
    call linearinterpfromlocations(x1grid, ygrid(:,x2lowerloc), x1lowerloc, x1val, dim1, extrapbot, extraptop, &
        newygrid(1))
    call linearinterpfromlocations(x1grid, ygrid(:,(x2lowerloc+1)), x1lowerloc, x1val, dim1, extrapbot, extraptop, &
        newygrid(2))

    ! Now generate a new grid from from the x2 bounds (a new ygrid has already bee generated)
    newxgrid(1) = x2grid(x2lowerloc)
    newxgrid(2) = x2grid(x2lowerloc+1)

    if (x2lowerloc.eq.0) then
        newxgridlowerloc = 0
    else if (x2lowerloc.eq.dim2) then
        newxgridlowerloc = 2
    else
        newxgridlowerloc = 1
    end if
    call linearinterpfromlocations(newxgrid, newygrid, newxgridlowerloc, x2val, 2, extrapbot, extraptop, yval)

    end subroutine linearinterp2_withextrap


    ! ---------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------!
    ! Extrapolation is by default here. There is no option to use the bottom value (as I have
    ! in the one-dimensional linear interpolation).

    subroutine linearinterp2(x1grid, x2grid, dim1, dim2, x1val, x2val, yval, ygrid)
    !use types
    implicit none
    ! Arguments
    integer, intent (in) :: dim1
    integer, intent (in) :: dim2
    real (kind=rk), intent (in) :: x1grid(dim1)
    real (kind=rk), intent (in) :: x2grid(dim2)
    real (kind=rk), intent (in) :: x1val
    real (kind=rk), intent (in) :: x2val
    real (kind=rk), intent (out) :: yval
    real (kind=rk), intent (in) :: ygrid(dim1, dim2)

    ! For programme
    integer :: x1lowerloc
    integer :: x2lowerloc
    integer :: scratch_nearestloc
    integer :: newxgridlowerloc
    real (kind=rk) :: newygrid(2)
    real (kind=rk) :: newxgrid(2)

    !Have commented out for speed
    if ((size(x1grid).ne.dim1) .or. (size(x2grid).ne.dim2)) then
        write (*, *) 'An argument in linearinterp2 is not of the specified length'
        stop
    end if

    call findingridr(x1grid, dim1, x1val, scratch_nearestloc, x1lowerloc)

    if (x1lowerloc.eq.0) then
        x1lowerloc = 1
    end if
    if (x1lowerloc.eq.dim1) then
        x1lowerloc = dim1 - 1
    end if

    ! Find the square
    call findingridr(x2grid, dim2, x2val, scratch_nearestloc, x2lowerloc)

    ! Find if any of the values are off the grid

    if (x2lowerloc.eq.0) then
        x2lowerloc = 1
    end if
    if (x2lowerloc.eq.dim2) then
        x2lowerloc = dim2 - 1
    end if




    !write(*,*) x1lowerloc
    !write(*,*) x2lowerloc
    !write(*,*) dim1
    !write(*,*) ygrid(:,x2lowerloc)
    !write(*,*) x1grid

    ! Now do linearinterpolation along the first dimension
    !        if (check1) then
    !        write(*,*) x1lowerloc, x2lowerloc
    !        write(*,*) dim1, dim2
    !        write(*,*) x1val
    !        write(*,*) ygrid
    !        stop
    !        end if

    call linearinterpfromlocations(x1grid, ygrid(:,x2lowerloc), x1lowerloc, x1val, dim1, 1, 1, newygrid(1))
    call linearinterpfromlocations(x1grid, ygrid(:,(x2lowerloc+1)), x1lowerloc, x1val, dim1, 1, 1, newygrid(2))


    ! Now generate a new grid from from the x2 bounds (a new ygrid has already bee generated)
    newxgrid(1) = x2grid(x2lowerloc)
    newxgrid(2) = x2grid(x2lowerloc+1)



    !Think this is redundant
    if (x2lowerloc.eq.0) then
        newxgridlowerloc = 0
    else if (x2lowerloc.eq.dim2) then
        newxgridlowerloc = 2
    else
        newxgridlowerloc = 1
    end if




    call linearinterpfromlocations(newxgrid, newygrid, newxgridlowerloc, x2val, 2, 1, 1, yval)

    end subroutine linearinterp2


    ! ---------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------!

    ! Subroutine for interpolating in one dimension given the points on the grid that bound the value for
    ! evaluation
    subroutine linearinterpfromlocations(xgrid, ygrid, xlowerloc, xval, length, extrapbot, extraptop, yval)
    implicit none
    integer, intent (in) :: length
    real (kind=rk), intent (in) :: xgrid(length)
    real (kind=rk), intent (in) :: ygrid(length)
    integer, intent (in) :: xlowerloc
    real (kind=rk), intent (in) :: xval
    integer, intent (in) :: extrapbot
    integer, intent (in) :: extraptop
    real (kind=rk), intent (out) :: yval


    ! For programme
    real (kind=rk) xgridmin
    real (kind=rk) xgridmax
    real (kind=rk) xgriddiffinv
    real (kind=rk) weightmin
    real (kind=rk) weightmax
    real (kind=rk) extrapslope
    integer :: xupperloc

    ! First, and easiest, get the interpolated value if xval is interior to the grid
    if (xlowerloc.ne.0 .and. xlowerloc.ne.length) then
        xupperloc = xlowerloc + 1
        if (ygrid(xupperloc) == ygrid(xlowerloc)) then
            yval = ygrid(xlowerloc)
            return
        end if
        xgridmin = xgrid(xlowerloc)
        xgridmax = xgrid(xupperloc)
        xgriddiffinv = 1/(xgridmax-xgridmin)

        weightmax = (xval-xgridmin)*xgriddiffinv
        weightmin = (xgridmax-xval)*xgriddiffinv

        ! Calculating the interpolated function values
        yval = weightmin*ygrid(xlowerloc) + weightmax*ygrid(xupperloc)
    else if (xlowerloc.eq.0) then
        if (extrapbot.eq.0) then
            yval = ygrid(1)
        else
            extrapslope = (ygrid(2)-ygrid(1))/(xgrid(2)-xgrid(1))
            yval = ygrid(1) + extrapslope*(xval-xgrid(1))
        end if

    else ! a is hu-eg so we extrapolate
        if (extraptop.eq.0) then
            yval = ygrid(length)
        else
            extrapslope = (ygrid(length)-ygrid(length-1))/(xgrid(length)-xgrid(length-1))
            yval = ygrid(length) + extrapslope*(xval-xgrid(length))
        end if
    end if


    end subroutine linearinterpfromlocations

    ! ---------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------!

    ! ---------------------------------------------------------------------------------------!
    ! Fully generic - need to write explanation
    subroutine gengrid(lowerbounds, upperbounds, numyears, gridpoints, growthgrid, assgrid)
    implicit none

    ! Arguments
    integer, intent (in) :: numyears
    integer, intent (in) :: gridpoints
    real (kind=rk), intent (in) :: lowerbounds(numyears)
    real (kind=rk), intent (in) :: upperbounds(numyears)
    real (kind=rk), intent (in) :: growthgrid

    real (kind=rk), intent (out) :: assgrid(gridpoints, (numyears+1))

    ! For program
    integer ixassnode, ixtoday
    real (kind=rk) restofgrid(gridpoints-8)

    ! Getting the asset grid
    do ixtoday = numyears, 1, -1
        ! I start the grid from the constraint plus the mininmum shock
        assgrid(1, ixtoday) = lowerbounds(ixtoday)
        assgrid(2, ixtoday) = assgrid(1, ixtoday) + 0.000001_rk
        assgrid(3, ixtoday) = assgrid(2, ixtoday) + 0.000003_rk
        assgrid(4, ixtoday) = assgrid(3, ixtoday) + 0.000005_rk
        assgrid(5, ixtoday) = assgrid(4, ixtoday) + 0.00001_rk
        assgrid(6, ixtoday) = assgrid(5, ixtoday) + 0.0001_rk
        assgrid(7, ixtoday) = assgrid(6, ixtoday) + 0.001_rk
        assgrid(8, ixtoday) = assgrid(7, ixtoday) + 0.005_rk

        call getnodes(restofgrid, assgrid(8,ixtoday)+0.01_rk, (upperbounds(ixtoday)), gridpoints-8, growthgrid)

        do ixassnode = 9, gridpoints ! Temp: At the moment I'm capping grid at a half of previous income -
            ! extrapolation beyond that
            assgrid(ixassnode, ixtoday) = restofgrid(ixassnode-8)
        end do

    end do ! time loop

    end subroutine gengrid
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------!
    ! Fully generic - need to write explanation

    subroutine gengrid_extra(lowerbounds, breaks1, breaks2, upperbounds, numyears, gridpoints1, gridpoints2, &
        gridpoints3, growthgrid1, growthgrid2, growthgrid3, gridpoints, assgrid)
    implicit none

    ! Arguments
    integer, intent (in) :: numyears
    integer, intent (in) :: gridpoints1
    integer, intent (in) :: gridpoints2
    integer, intent (in) :: gridpoints3
    integer, intent (in) :: gridpoints
    real (kind=rk), intent (in) :: lowerbounds(numyears)
    real (kind=rk), intent (in) :: breaks1(numyears)
    real (kind=rk), intent (in) :: breaks2(numyears)
    real (kind=rk), intent (in) :: upperbounds(numyears)
    real (kind=rk), intent (in) :: growthgrid1
    real (kind=rk), intent (in) :: growthgrid2
    real (kind=rk), intent (in) :: growthgrid3

    real (kind=rk), intent (out) :: assgrid(gridpoints, numyears)

    ! For program
    integer ixassnode, ixtoday, index
    real (kind=rk) restofgrid1(gridpoints1-2)
    real (kind=rk) restofgrid2(gridpoints2+1)
    real (kind=rk) restofgrid3(gridpoints3+1)

    ! Check gridpoints
    if ((gridpoints1+gridpoints2+gridpoints3).ne.gridpoints) then
        write (*, *) 'Error in gengrid_extra; gridpoints1 + gridpoints2 + gridpoints3 is not equal gridpoints'
        stop
    end if

    ! Getting the asset grid
    do ixtoday = 1, numyears, 1

        ! I start the grid from the constraint plus the mininmum shock
        assgrid(1, ixtoday) = lowerbounds(ixtoday)
        assgrid(2, ixtoday) = assgrid(1, ixtoday) + 0.000001_rk

        call getnodes(restofgrid1, assgrid(2,ixtoday)+0.01_rk, breaks1(ixtoday), gridpoints1-2, growthgrid1)
        do ixassnode = 3, gridpoints1 - 1, 1 ! Don't want to add in break1 - leave that till next one
            assgrid(ixassnode, ixtoday) = restofgrid1(ixassnode-2)
        end do

        call getnodes(restofgrid2, breaks1(ixtoday), breaks2(ixtoday), gridpoints2+1, growthgrid2)
        index = 1
        do ixassnode = (gridpoints1), (gridpoints1+gridpoints2-1), 1
            assgrid(ixassnode, ixtoday) = restofgrid2(index)
            index = index + 1
        end do

        call getnodes(restofgrid3, breaks2(ixtoday), upperbounds(ixtoday), gridpoints3+1, growthgrid3)
        index = 1
        do ixassnode = (gridpoints1+gridpoints2), gridpoints, 1
            assgrid(ixassnode, ixtoday) = restofgrid3(index)
            index = index + 1
        end do

    end do ! time loop

    end subroutine gengrid_extra
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!

    ! This function returns a grid where:
    ! 1. The grid starts at firstpoint (input)
    ! 2. The grid ends at lastpoint (input)
    ! 3. The number of nodes is numberofnodes (input)
    ! 4. The growth rate of the step size is growthrate (input)

    ! This function is fully generic

    subroutine getnodes(nodes, firstpoint, lastpoint, numberofnodes, growthrate)
    implicit none

    ! Arguments
    real (kind=rk), intent (in) :: firstpoint
    real (kind=rk), intent (in) :: lastpoint
    real (kind=rk), intent (in) :: growthrate
    integer, intent (in) :: numberofnodes
    real (kind=rk), intent (out) :: nodes(numberofnodes)

    ! For programme
    real (kind=rk) step
    real (kind=rk) denom
    integer ix

    ! Check inputs
    if (firstpoint.ge.lastpoint) then
        write (*, *) firstpoint, lastpoint
        write (*, *) 'In getnodes, firstpoint must be strictly less than lastpoint'
        write (*, *) 'The +10._rk is clumsy in gengrid_Extra'
        stop
    end if

    ! First I need, conditional, on the span of the space, the number of gridpoints, and the growth
    ! rate of the space between successive gridpoints, to work out the size of the first step
    denom = 0
    do ix = 1, numberofnodes - 1, 1
        denom = denom + (1+growthrate)**(ix-1)
    end do

    step = (lastpoint-firstpoint)/denom

    nodes(1) = firstpoint
    do ix = 2, numberofnodes, 1
        nodes(ix) = nodes(ix-1) + step
        step = step*(1+growthrate)
    end do

    end subroutine getnodes
    ! ---------------------------------------------------------------------------------------------------------!
    ! Function to get numerical first derivative. Generic other than use of linearinterp1 and I've hardcoded
    ! last two arguments in linearinterp1
    real (kind=rk) function deriv(x, fx, x0, length)
    implicit none

    ! Arguments
    integer, intent (in) :: length
    real (kind=rk), intent (in) :: x(length)
    real (kind=rk), intent (in) :: fx(length)
    real (kind=rk), intent (in) :: x0

    ! For progamme
    real (kind=rk) :: forderiv1, forderiv2

    call linearinterp1(x, fx, length, x0+0.0000001_rk, forderiv1, 0, 1)
    call linearinterp1(x, fx, length, x0-0.0000001_rk, forderiv2, 0, 1)
    deriv = (forderiv1-forderiv2)/(0.0000002_rk)

    end function deriv
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    ! Function to get the numerical derivative of a function when that function is defined on a 2 dimensional
    ! grid (with two arguments). Generic except for the reliance on linearinterp2
    ! NB: The derivative is with respect to the variable in x1grid at x1val
    real (kind=rk) function deriv_2d(x1grid, x2grid, ygrid, dim1, dim2, x1val, x2val)

    implicit none

    ! Arguments
    integer, intent (in) :: dim1
    integer, intent (in) :: dim2
    real (kind=rk), intent (in) :: x1grid(dim1)
    real (kind=rk), intent (in) :: x2grid(dim2)
    real (kind=rk), intent (in) :: ygrid(dim1, dim2)
    real (kind=rk), intent (in) :: x1val
    real (kind=rk), intent (in) :: x2val


    ! For progamme
    real (kind=rk) :: forderiv1, forderiv2

    call linearinterp2(x1grid, x2grid, dim1, dim2, (x1val+0.0000001_rk), x2val, forderiv1, ygrid)
    call linearinterp2(x1grid, x2grid, dim1, dim2, (x1val-0.0000001_rk), x2val, forderiv2, ygrid)
    deriv_2d = (forderiv1-forderiv2)/(0.0000002_rk)

    end function deriv_2d

    ! ---------------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------------!
    ! A subroutine for finding  the nodes and the values on a real grid that bound a value off the grid
    ! It is fully generic
    subroutine findingridr(xgrid, length, xval, nearestloc, lowerloc)
    implicit none
    ! Arguments
    integer, intent (in) :: length
    real (kind=rk), intent (in) :: xgrid(length)
    real (kind=rk), intent (in) :: xval

    integer, intent (out) :: nearestloc
    integer, intent (out) :: lowerloc

    ! For programme
    integer :: jl
    integer :: jm
    integer :: ju
    integer :: upperloc
    real (kind=rk) :: lowerval
    real (kind=rk) :: upperval
    intrinsic size

    ! Checking inputs
    if ((size(xgrid).ne.length)) then
        write (*, *) 'An argument in findingridr is not of the specified length'
        stop
    end if

    ! This block of code finds the two points on the xgrid between which the x value lies
    jl = 0
    ju = length + 1
    do while (ju-jl.gt.1)
        jm = (ju+jl)/2
        if (xval.gt.xgrid(jm)) then
            jl = jm
        else
            ju = jm
        end if
    end do

    lowerloc = jl
    upperloc = jl + 1

    if (lowerloc.eq.0) then
        lowerval = xgrid(1)
        upperval = xgrid(1)
    else if (lowerloc.eq.length) then
        lowerval = xgrid(length)
        upperval = xgrid(length)
    else
        lowerval = xgrid(lowerloc)
        upperval = xgrid(upperloc)
    end if

    if ((((upperval-xval).lt.(xval-lowerval)) .or. (lowerloc.eq.0)) .and. (lowerloc.ne.length)) then
        nearestloc = upperloc
    else
        nearestloc = lowerloc
    end if

    if ((nearestloc.lt.0) .or. (nearestloc.gt.length)) then
        write (*, *) 'error in findingr - nearestloc is outside the grid'
        stop
    end if

    end subroutine findingridr


    ! ---------------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------------!
    ! This subroutine takes a vector of shock values that represent a stationary Markov process
    ! and generates profiles of shock values using the associated cdf
    ! Fully generic
    !      subroutine simmarkovshocks(zvect, cdf, length, numperiods, nosims, startatmiddle, shockslist, shocksval)
    !        implicit none
    !
    !! Arguments
    !        integer, intent (in) :: length ! These are the number the number of discrete points in the shock vector
    !        real (kind=rk), intent (in) :: zvect(length)
    !        real (kind=rk), intent (in) :: cdf(length, length)
    !        integer, intent (in) :: numperiods
    !        integer, intent (in) :: nosims
    !        logical, intent (in) :: startatmiddle
    !
    !        real (kind=rk), intent (out) :: shocksval(numperiods, nosims)
    !        integer, intent (out) :: shockslist(numperiods, nosims)
    !
    !! For use by subroutine
    !        integer, parameter :: ntforstat = 100 ! No of periods in simulation to findi stationary distribution
    !        real (kind=rk) :: uniformvar(nosims*numperiods)
    !        real (kind=rk) :: uniformvarforstat(nosims*ntforstat)
    !        integer :: zeroshockloc
    !        integer :: todayshockloc
    !        integer :: yesterdayshockloc
    !        integer :: statshockslist(nosims) ! Using 100 periods to find the stationary distribution (i.e. to get the
    !! initial shock for each person)
    !        integer igen, ifail, seed(4)
    !        integer :: lowerboundcdf
    !! Indexes
    !        integer ixtime, ixss, index
    !
    !! Scratch
    !        integer scratch1
    !
    !! external
    !        external g05kbf
    !        external g05lgf
    !
    !! Getting random (uniformly distributed) numbers
    !! This is for the burn in
    !        seed(1) = 559987891
    !        seed(2) = 357865655
    !        seed(3) = 357876222
    !        seed(4) = 981656871
    !        igen = 0
    !        call g05kbf(igen, seed)
    !        igen = 0
    !        ifail = 0
    !        call g05lgf(0.0_rk, 1.0_rk, nosims*ntforstat, uniformvarforstat, igen, seed, ifail)
    !
    !        zeroshockloc = length/2 + 1 ! This will give the zero shock as long as length is odd
    !         !Check
    !        if ((zeroshockloc.lt.1).or.(zeroshockloc.gt.length)) then
    !          write (*,*) 'error in simmarkovshocks'
    !          write (*,*) 'zeroshockloc shock is not in range'
    !          write (*,*) 'Is length of shock vector 1'
    !          stop
    !        end if
    !
    !
    !        index = 0 ! Will run through the random number vector
    !
    !! First finding the stationary distribution of shocks
    !        do ixss = 1, nosims
    !          do ixtime = 1, ntforstat
    !            if (ixtime.eq.1) then
    !              yesterdayshockloc = zeroshockloc
    !            else
    !              yesterdayshockloc = todayshockloc
    !            end if
    !            index = index + 1
    !
    !            call findingridr(cdf(:,yesterdayshockloc), length, uniformvarforstat(index), scratch1, lowerboundcdf)
    !            todayshockloc = lowerboundcdf + 1
    !          end do ! ixtime
    !          statshockslist(ixss) = todayshockloc
    !        end do ! ixss
    !
    !! Burn in over, lets now find the shocks for our simulatiion
    !        seed(1) = 465465465
    !        seed(2) = 894632231
    !        seed(3) = 123654982
    !        seed(4) = 888888888
    !        igen = 0
    !        call g05kbf(igen, seed)
    !        igen = 0
    !        ifail = 0
    !        call g05lgf(0.0_rk, 1.0_rk, nosims*numperiods, uniformvar, igen, seed, ifail)
    !
    !        index = 0
    !        do ixss = 1, nosims
    !          do ixtime = 1, numperiods, 1
    !            if (ixtime.eq.1) then
    !              if (startatmiddle) then ! We want to start everyone from the zero shock
    !                yesterdayshockloc = zeroshockloc
    !              else
    !                yesterdayshockloc = statshockslist(ixss) ! Take the last shock for that individual in the last of
    !! period of the
    !              end if
    !            else ! of the synthethic period used above to find the stationary distribution
    !              yesterdayshockloc = todayshockloc
    !            end if
    !            index = index + 1
    !
    !            call findingridr(cdf(:,yesterdayshockloc), length, uniformvar(index), scratch1, lowerboundcdf)
    !            todayshockloc = lowerboundcdf + 1
    !            shockslist(ixtime, ixss) = todayshockloc
    !            shocksval(ixtime, ixss) = zvect(todayshockloc) ! To look at the actual shocks
    !          end do ! ixtime
    !        end do ! ixss
    !      end subroutine simmarkovshocks
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!

    ! This subroutine takes a vector of shock values that represent a stationary Markov process
    ! and generates profiles of shock values using the associated cdf
    ! Fully generic
    !      subroutine simmarkovshocks_nonstat(zvect, cdf, length, numperiods, startingshock,nosims, shockslist, shocksval, iterateForStarting)
    !        implicit none
    !
    !! Arguments
    !        integer, intent (in) :: length ! These are the number the number of discrete points in the shock vector
    !        integer, intent (in) :: numperiods
    !        integer, intent (in) :: startingshock
    !        real (kind=rk), intent (in) :: zvect(length, numperiods)
    !        real (kind=rk), intent (in) :: cdf(length, length, numperiods)
    !        integer, intent (in) :: nosims
    !
    !        real (kind=rk), intent (out) :: shocksval(numperiods, nosims)
    !        integer, intent (out) :: shockslist(numperiods, nosims)
    !
    !        logical, intent(in), optional :: iterateForStarting ! Do we start everyone in period zero from startingshock or do we iterate on period 1
    !
    !! For use by subroutine
    !        real (kind=rk) :: uniformvar(nosims*numperiods)
    !        integer :: zeroshockloc
    !        integer :: todayshockloc
    !        integer :: yesterdayshockloc
    !        integer igen, ifail, seed(4)
    !        integer :: lowerboundcdf
    !        integer, parameter :: ntforstat = 100 ! No of periods in simulation to findi stationary distribution
    !        real (kind=rk) :: uniformvarforstat(nosims*ntforstat)
    !        integer :: statshockslist(nosims) ! Using 100 periods to find the stationary distribution (i.e. to get the
    !
    !        logical :: doIterate ! Do we start everyone in period zero from startingshock or do we iterate on period 1 - need new variable as iterateForStarting is optional
    !
    !! Indexes
    !        integer ixtime, ixss, index
    !
    !! Scratch
    !        integer scratch1
    !
    !! external
    !        external g05kbf
    !        external g05lgf
    !
    !        !Check
    !        if ((startingshock.lt.1).or.(startingshock.gt.length)) then
    !          write (*,*) 'error in simmarkovshocks_nonstat'
    !          write (*,*) 'starting shock is not in range'
    !          stop
    !        end if
    !
    !        zeroshockloc =  startingshock
    !
    !! Maybe I need to  iterate on the first year's observations - to, for example, get the inititial distribution of
    !! women working
    !        if (present(iterateForStarting)) then
    !            if (iterateForStarting) then
    !                doIterate = .true.
    !            else
    !                doIterate = .false.
    !            end if
    !         else
    !                doIterate = .false.
    !         end if
    !
    !if (doIterate) then
    !! Getting random (uniformly distributed) numbers
    !! This is for the burn in
    !        seed(1) = 559987891
    !        seed(2) = 357865655
    !        seed(3) = 357876222
    !        seed(4) = 981656871
    !        igen = 0
    !        call g05kbf(igen, seed)
    !        igen = 0
    !        ifail = 0
    !        call g05lgf(0.0_rk, 1.0_rk, nosims*ntforstat, uniformvarforstat, igen, seed, ifail)
    !
    !        index = 0 ! Will run through the random number vector
    !        do ixss = 1, nosims
    !          do ixtime = 1, ntforstat
    !            if (ixtime.eq.1) then
    !              yesterdayshockloc = zeroshockloc
    !            else
    !              yesterdayshockloc = todayshockloc
    !            end if
    !            index = index + 1
    !
    !            call findingridr(cdf(yesterdayshockloc,:,1), length, uniformvarforstat(index), scratch1, lowerboundcdf)
    !            todayshockloc = lowerboundcdf + 1
    !          end do ! ixtime
    !          statshockslist(ixss) = todayshockloc
    !        end do ! ixss
    !end if ! if (doiterate)
    !
    !        index = 0 ! Will run through the random number vector
    !
    !! Lets now find the shocks for our simulatiion
    !        seed(1) = 465465465
    !        seed(2) = 894632231
    !        seed(3) = 123654982
    !        seed(4) = 888888888 !cf. Nick Leeson
    !        igen = 0
    !        call g05kbf(igen, seed)
    !        igen = 0
    !        ifail = 0
    !        call g05lgf(0.0_rk, 1.0_rk, nosims*numperiods, uniformvar, igen, seed, ifail)
    !
    !        index = 0
    !        do ixss = 1, nosims
    !          do ixtime = 1, numperiods, 1
    !            if (ixtime.eq.1) then
    !                if (doIterate) then
    !                  yesterdayshockloc = statshockslist(ixss)
    !                else
    !                    yesterdayshockloc = zeroshockloc
    !                end if
    !            else
    !              yesterdayshockloc = todayshockloc
    !            end if
    !            index = index + 1
    !
    !            call findingridr(cdf(yesterdayshockloc,:,ixtime), length, uniformvar(index), scratch1, lowerboundcdf)
    !            todayshockloc = lowerboundcdf + 1
    !            shockslist(ixtime, ixss) = todayshockloc
    !            shocksval(ixtime, ixss) = zvect(todayshockloc, ixtime) ! To look at the actual shocks
    !          end do ! ixtime
    !        end do ! ixss
    !
    !      end subroutine simmarkovshocks_nonstat
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    ! WOuld be a good idea to write this subroutine at some stage - to test ar1
    ! subroutine testar1(
    ! real(kind=rk) :: shocksval(numyears,howmanysims)
    ! integer :: shockslist(numyears,howmanysims)

    ! call getprocesses(numyears,order,params,ar1rho,ar1stdev, &
    ! ar1points,detincome,ar1shocks,ar1pmat,ar1cdf)

    ! call simmarkovshocks_nonstat(ar1shocks,ar1cdf,ar1points,numyears, howmanysims, shockslist,shocksval)
    ! open (unit=555,file='Y:\Cormac\temp.txt',recl=100000,action='write')
    ! write (555,'(999f15.7)') transpose(shocksval)
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------------!
    !subroutine simContAr1(mean, rho, stdev, stdevinit, years, sims, ageStopFanning, seed, draws)
    !  implicit none
    !  real (kind=rk), intent(in) :: mean, rho, stdev, stdevinit
    !  integer, intent(in) ::years, sims, seed, ageStopFanning
    !  real (kind=rk), intent(out) :: draws(years, sims)
    !
    !  ! For prog
    !  real (kind=rk) :: innov(years, sims), innov1(sims)
    !  logical, parameter :: drawNorm = .true.
    !  integer :: ixY, ixS, seed1
    !
    !  seed1 = (int(seed/18_rk) + 150) * 2
    !
    !  if (stdevinit.gt.0.0_rk) then
    !    call getdraws(1, sims, sims, seed1, innov1, drawNorm, 0.0_rk, stdevinit*stdevinit)
    !  end if
    !
    !  call getdraws(years, sims, years*sims, seed, innov, drawNorm, 0.0_rk, stdev*stdev)
    !
    !
    !  do ixS = 1, sims
    !  if (stdevinit.gt.0.0_rk) then
    !    draws(1, ixS) = innov1(ixS)
    !  else
    !    draws(1, ixS) = innov(1, ixS)
    !  end if
    !
    !    do ixY = 2, years
    !        if ( (ixY+ 19).lt.ageStopFanning) then
    !            draws(ixY, ixS) = rho * draws((ixY-1), ixS) + innov(ixY, ixS)
    !        else
    !            draws(ixY, ixS) = draws(ixY - 1, ixS)
    !        end if
    !    end do ! ixY
    !  end do ! ixS
    !
    !
    !end subroutine simContAr1
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------------!


    function transpose1dp(vector, length)
    implicit none
    integer, intent (in) :: length
    real (kind=rk), intent (in) :: vector(length)

    real (kind=rk) :: transpose1dp(1, length)

    integer :: ix

    do ix = 1, length, 1
        transpose1dp(1, ix) = vector(ix)
    end do

    end function transpose1dp
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    function transpose1int(vector, length)
    implicit none
    integer, intent (in) :: length
    integer, intent (in) :: vector(length)

    integer :: transpose1int(1, length)

    integer :: ix

    do ix = 1, length, 1
        transpose1int(1, ix) = vector(ix)
    end do

    end function transpose1int
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!

    ! Check stochastic discount factor with five state variables
    subroutine checksdf5(sdf, licon, dim1, dim2, dim3, dim4, dim5, fail, failprop)
    implicit none
    ! Arguments
    integer, intent (in) :: dim1
    integer, intent (in) :: dim2
    integer, intent (in) :: dim3
    integer, intent (in) :: dim4
    integer, intent (in) :: dim5
    real (kind=rk), intent (in) :: sdf(dim1, dim2, dim3, dim4, dim5)
    integer, intent (in) :: licon(dim1, dim2, dim3, dim4, dim5)
    integer, intent (out) :: fail(4)
    real (kind=rk), intent (out) :: failprop(4)

    ! For program
    ! indexes
    integer ixdim1, ixdim2, ixdim3, ixdim4, ixdim5, countinplay

    real (kind=rk) sdfinplay

    fail(1) = 0
    fail(2) = 0
    fail(3) = 0
    fail(4) = 0
    countinplay = 0 ! This is the number of stochastic discount factors we're testing (i.e. not
    ! liquidity constrained)

    do ixdim5 = dim5, 1, -1
        do ixdim4 = 1, dim4, 1
            do ixdim3 = 1, dim3, 1
                do ixdim2 = 1, dim2, 1
                    do ixdim1 = 1, dim1, 1
                        if ((licon(ixdim1,ixdim2,ixdim3,ixdim4,ixdim5).ne.1) .and. (licon(ixdim1,ixdim2,ixdim3,ixdim4, &
                            ixdim5).ne.-2)) then
                        sdfinplay = sdf(ixdim1, ixdim2, ixdim3, ixdim4, ixdim5)
                        countinplay = countinplay + 1
                        if ((sdfinplay.lt.0.995_rk) .or. (sdfinplay.gt.1.0005_rk)) then
                            fail(1) = fail(1) + 1
                            if ((sdfinplay.lt.0.99_rk) .or. (sdfinplay.gt.1.01_rk)) then
                                fail(2) = fail(2) + 1
                                if ((sdfinplay.lt.0.975_rk) .or. (sdfinplay.gt.1.025_rk)) then
                                    fail(3) = fail(3) + 1
                                    if ((sdfinplay.lt.0.95_rk) .or. (sdfinplay.gt.1.05_rk)) then
                                        fail(4) = fail(4) + 1
                                    end if
                                end if
                            end if
                        end if
                        end if ! If not liquidity constrained
                    end do ! ixdim1
                end do ! ixdim2
            end do ! ixdim3
        end do ! ixdim4
    end do ! ixdim4

    failprop = real(fail, rk)/real(countinplay, rk)
    end subroutine checksdf5
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    ! Check stochastic discount factor with six state variables
    subroutine checksdf6(sdf, licon, dim1, dim2, dim3, dim4, dim5, dim6, indincnonnull, indincnoncon, fail, failprop)
    implicit none
    ! Arguments
    integer, intent (in) :: dim1
    integer, intent (in) :: dim2
    integer, intent (in) :: dim3
    integer, intent (in) :: dim4
    integer, intent (in) :: dim5
    integer, intent (in) :: dim6
    real (kind=rk), intent (in) :: sdf(dim1, dim2, dim3, dim4, dim5, dim6)
    integer, intent (in) :: licon(dim1, dim2, dim3, dim4, dim5, dim6)
    integer, intent (out) :: fail(4)
    integer, intent (out) :: indincnonnull
    integer, intent (out) :: indincnoncon
    real (kind=rk), intent (out) :: failprop(4)

    ! For program
    ! indexes
    integer ixdim1, ixdim2, ixdim3, ixdim4, ixdim5, ixdim6

    real (kind=rk) sdfinplay

    fail(1) = 0
    fail(2) = 0
    fail(3) = 0
    fail(4) = 0
    indincnoncon = 0 ! This is the number of stochastic discount factors we're testing (i.e. not
    ! liquidity constrained)
    indincnonnull = 0 ! This is the number of nonnull points in the state space
    do ixdim6 = 1, dim6, 1
        do ixdim5 = 1, dim5, 1
            do ixdim4 = 1, dim4, 1
                do ixdim3 = 1, dim3, 1
                    do ixdim2 = 1, dim2, 1
                        do ixdim1 = 1, dim1, 1
                            if (licon(ixdim1,ixdim2,ixdim3,ixdim4,ixdim5,ixdim6).ne.2) then
                                indincnonnull = indincnonnull + 1
                                if (licon(ixdim1,ixdim2,ixdim3,ixdim4,ixdim5,ixdim6).ne.1) then
                                    sdfinplay = sdf(ixdim1, ixdim2, ixdim3, ixdim4, ixdim5, ixdim6)
                                    indincnoncon = indincnoncon + 1
                                    if ((sdfinplay.lt.0.995_rk) .or. (sdfinplay.gt.1.0005_rk)) then
                                        fail(1) = fail(1) + 1
                                        if ((sdfinplay.lt.0.99_rk) .or. (sdfinplay.gt.1.01_rk)) then
                                            fail(2) = fail(2) + 1
                                            if ((sdfinplay.lt.0.975_rk) .or. (sdfinplay.gt.1.025_rk)) then
                                                fail(3) = fail(3) + 1
                                                if ((sdfinplay.lt.0.95_rk) .or. (sdfinplay.gt.1.05_rk)) then
                                                    fail(4) = fail(4) + 1
                                                end if
                                            end if
                                        end if
                                    end if
                                end if
                            end if ! If not liquidity constrained
                        end do ! ixdim1
                    end do ! ixdim2
                end do ! ixdim3
            end do ! ixdim4
        end do ! ixdim5
    end do ! ixdim6

    failprop = real(fail, rk)/real(indincnoncon, rk)
    end subroutine checksdf6
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    ! Check stochastic discount factor with with one-dimensional vector
    subroutine checksdf(sdf, licon, dim1, numtests, testvalues, indincnonnull, indincnoncon, countfail, numfails, &
        failprop)
    implicit none

    integer, intent (in) :: dim1
    integer, intent (in) :: numtests
    real (kind=rk), intent (in) :: sdf(dim1)
    integer, intent (in) :: licon(dim1)
    real (kind=rk), intent (in) :: testvalues(numtests)
    integer, intent (out) :: indincnonnull(dim1)
    integer, intent (out) :: indincnoncon(dim1)
    integer, intent (out) :: countfail(dim1, numtests)
    integer, intent (out) :: numfails(numtests)
    real (kind=rk), intent (out) :: failprop(numtests)

    integer :: ixtest, numnoncons

    where (licon.ne.2)
        indincnonnull = 1
    elsewhere
        indincnonnull = 0
    end where

    where ((licon.ne.2) .and. (licon.ne.1))
        indincnoncon = 1
    elsewhere
        indincnonnull = 0
    end where

    do ixtest = 1, numtests, 1
        where ((indincnoncon.eq.1) .and. ((sdf.lt.(1.0_rk-testvalues(ixtest))) .or. (sdf.gt. &
            (1.0_rk+testvalues(ixtest)))))
        countfail(:, ixtest) = 1
        elsewhere
            countfail(:, ixtest) = 0
        end where
    end do

    !Getting sums
    numfails = sum(countfail, 1)
    numnoncons = sum(indincnoncon, 1)

    failprop = numfails/numnoncons

    end subroutine checksdf
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!


    ! A subroutine that takes a (dp) number and assess where the number is far away from another number
    ! Far away is taken as a vector of four (dp) distances. The subroutine returns an integer (dim of 4) vector
    ! with ones where the supplied value is far away and zeros otherwise
    ! subroutine testhowfar(numbertotest, centralval, vecofepsilons, fail)
    ! implicit none



    ! if ((numbertotest<(centralval - vecofepsilons(1) ) .or. (centralval + vecofepsilons(1))) then
    ! fail(1) = fail(1) + 1
    ! if ((numbertotest<(centralval - vecofepsilons(2) ) .or. (centralval + vecofepsilons(2))) then
    ! fail(2) = fail(2) + 1
    ! if ((numbertotest<(centralval - vecofepsilons(3) ) .or. (centralval + vecofepsilons(3))) then
    ! fail(3) = fail(3) + 1
    ! if ((numbertotest<(centralval - vecofepsilons(4) ) .or. (centralval + vecofepsilons(4))) then
    ! fail(4) = fail(4) + 1
    ! end if
    ! end if
    ! end if
    ! end if


    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    ! A routine that takes a vector representing a PDF and returns a vector representing a CDF
    subroutine pdftocdf(length, pdf, cdf)
    implicit none
    integer, intent (in) :: length
    real (kind=rk), intent (in) :: pdf(length)
    real (kind=rk), intent (out) :: cdf(length)

    ! For programme
    integer :: ixlen
    real (kind=rk) :: sumforcdf

    if ((sum(pdf).lt.0.999) .or. (sum(pdf).gt.1.001)) then
        write (*, *) 'error in pdftocdf: probabilities in pdf do not sum to 1'
        write (*, *) 'sum is', sum(pdf)
        stop
    end if

    ! Get CDF
    sumforcdf = 0
    do ixlen = 1, length
        sumforcdf = sumforcdf + pdf(ixlen)
        cdf(ixlen) = sumforcdf
    end do


    if ((cdf(length).lt.1.0_rk).and.(cdf(length).gt.0.999)) then
        cdf(length) = 1.0_rk
    end if

    end subroutine pdftocdf
    !---------------------------------------------------------------------------------------------------------!
    !---------------------------------------------------------------------------------------------------------!
    !Get uniform draws using the NAG libraries
    !      subroutine getdraws(length1, length2, totalrandom, seed, draws, normal, mu, var)
    !        implicit none
    !
    !        integer, intent (in) :: length1, length2
    !        integer, intent (in) :: seed
    !        integer, intent (in) :: totalrandom
    !        logical, intent (in), optional :: normal
    !        real (kind=rk), intent (in), optional :: mu
    !        real (kind=rk), intent (in), optional :: var
    !        real (kind=rk), intent (out) :: draws(length1, length2)
    !
    !!For programme
    !        integer :: genid
    !        integer :: subid !Says use the NAG basic generator
    !        integer :: lseed
    !        integer, parameter :: lstate = 633
    !        real (kind=rk) :: drawsvec(totalrandom)
    !        integer :: state(lstate)
    !        integer :: ifail
    !        external :: g05kff, g05saf
    !
    !!Check that if we have requested normal variables, mu and var are provided
    !        if (present(normal)) then
    !          if ((normal) .and. (.not. (present(mu) .and. present(var)))) then
    !            write (*, *) 'error in getdraws. If normal is .true.,  mean and variance must be supplied'
    !            stop
    !          end if
    !        end if
    !
    !!Check length1 * length2 = totalrandom
    !        if ((length1*length2).ne.totalrandom) then
    !          write (*, *) 'error in getunifdraws: length1 * length2 = totalrandom'
    !          stop
    !        end if
    !
    !        genid = 1 !Says use the NAG basic generator; do not change - will not work otherwise; need to write new subroutine if want to use more genera
    ! !
    !!
    !!
    !!l
    !        subid = 1 !This is irrelavant when we're using the NAG basic generator (i.e. genid = 1)
    !        lseed = 1 !size of seed
    !        ifail = 0
    !
    !! Initialize the generator to a repeatable sequence
    !        call g05kff(genid, subid, seed, lseed, state, lstate, ifail)
    !
    !        ifail = 0
    !
    !        if (present(normal) .and. (normal)) then
    !          call g05skf(totalrandom, mu, var, state, drawsvec, ifail)
    !        else
    !          call g05saf(totalrandom, state, drawsvec, ifail)
    !        end if
    !
    !        draws = reshape(drawsvec, (/length1,length2/) )
    !
    !
    !      end subroutine getdraws

    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !A function that takes an array and a probability distribution over it and works out the expected value
    !Could just matmau, but this avoids having to transpose
    !NBB: Not yet tested - not currently using
    function getexpec(array, probdist, length)
    implicit none
    integer, intent (in) :: length
    real (kind=rk), intent (in) :: array(length)
    real (kind=rk), intent (in) :: probdist(length)
    real (kind=rk) :: getexpec

    ! For programme
    integer ixlen

    getexpec = 0.0_rk
    do ixlen = 1, length, 1
        getexpec = getexpec + array(ixlen)*probdist(ixlen)
    end do

    end function getexpec
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !Function that gets the squared distance between two points
    function dist2(x, y)
    implicit none
    real (kind=rk), intent (in) :: x, y
    real (kind=rk) :: dist2

    dist2 = (x-y)*(x-y)

    end function dist2
    !---------------------------------------------------------------------------------------------------------!
    !---------------------------------------------------------------------------------------------------------!
    ! This is based on code given in Miranda and Fackler:
    ! "Applied Computaiontal Economics and Finance"
    ! I have corrected an error contained in the code that they produce

    real (kind=rk) function golden_generic(a, b, xmax, func,tol,show)
    implicit none

    ! Generic stuff
    ! Arguments
    real (kind=rk), intent (in) :: a
    real (kind=rk), intent (in) :: b
    real (kind=rk) :: func
    real (kind=rk), intent (out) :: xmax
    logical, intent(in) :: show

    ! For program
    real (kind=rk) :: x1
    real (kind=rk) :: x2
    real (kind=rk) :: f1
    real (kind=rk) :: f2
    real (kind=rk) :: diff
    real (kind=rk) :: goldenalpha1
    real (kind=rk) :: goldenalpha2
    real (kind=rk) :: tol
    parameter (goldenalpha2=0.611003399_rk)
    parameter (goldenalpha1=1.0_rk-goldenalpha2)
    !parameter (tol=1.d-10)

    x1 = a + goldenalpha1*(b-a)
    x2 = a + goldenalpha2*(b-a)

    f1 = func(x1)
    f2 = func(x2)
    diff = goldenalpha1*goldenalpha2*(b-a)
    do while (diff.gt.tol)
        diff = diff*goldenalpha2
        if (show) WRITE(*,*)  'Error is: ', diff
        if (f2.lt.f1) then
            x2 = x1
            x1 = x1 - diff

            f2 = f1
            f1 = func(x1)
        else
            x1 = x2
            x2 = x2 + diff

            f1 = f2
            f2 = func(x2)
        end if
    end do

    if (f2.gt.f1) then
        xmax = x2
        golden_generic = f2
    else
        xmax = x1
        golden_generic = f1
    end if

    end function golden_generic

    !---------------------------------------------------------------------------------------------------------!
    !---------------------------------------------------------------------------------------------------------!
    !A subroutine for trilinear interpolation
    !A nice method here (that I'm not using) method at: http://paulbourke.net/miscellaneous/interpolation/
    !This has been tested with the subroutine test below
    subroutine linearinterp3(x1grid, x2grid, x3grid, dim1, dim2, dim3, x1val, x2val, x3val, yval, ygrid)
    !use types
    implicit none

    ! Arguments
    integer, intent (in) :: dim1
    integer, intent (in) :: dim2
    integer, intent (in) :: dim3
    real (kind=rk), intent (in) :: x1grid(dim1)
    real (kind=rk), intent (in) :: x2grid(dim2)
    real (kind=rk), intent (in) :: x3grid(dim3)
    real (kind=rk), intent (in) :: ygrid(dim1, dim2, dim3)
    real (kind=rk), intent (in) :: x1val
    real (kind=rk), intent (in) :: x2val
    real (kind=rk), intent (in) :: x3val
    real (kind=rk), intent (out) :: yval

    ! For programme
    integer :: x3lowerloc
    integer :: scratch_nearestloc
    integer :: newxgridlowerloc
    real (kind=rk) :: newygrid(2)
    real (kind=rk) :: newxgrid(2)


    if ((size(x1grid).ne.dim1) .or. (size(x2grid).ne.dim2) .or. (size(x3grid).ne.dim3)) then
        write (*, *) 'An argument in linearinterp3 is not of the specified length'
        stop
    end if

    ! Locate x3val on the grid
    call findingridr(x3grid, dim3, x3val, scratch_nearestloc, x3lowerloc)

    ! Find if any of the values are off the grid
    if (x3lowerloc.eq.0) then
        x3lowerloc = 1
    end if

    if (x3lowerloc.eq.dim3) then
        x3lowerloc = dim3 - 1
    end if




    !Do two two bilinear interpolations over x1grid and x2grid keeping  x3 at x3lowerloc and x3lowerloc + 1
    call linearinterp2(x1grid, x2grid, dim1, dim2, x1val, x2val, newygrid(1), ygrid(:,:,x3lowerloc))
    call linearinterp2(x1grid, x2grid, dim1, dim2, x1val, x2val, newygrid(2), ygrid(:,:,(x3lowerloc+1)))




    ! Now to a linear interpolation on these found vlues
    ! First generate a new grid from from the x2 bounds (a new ygrid has already bee generated)
    newxgrid(1) = x3grid(x3lowerloc)
    newxgrid(2) = x3grid(x3lowerloc+1)

    !Think this is redundant
    if (x3lowerloc.eq.0) then
        newxgridlowerloc = 0
    else if (x3lowerloc.eq.dim3) then
        newxgridlowerloc = 2
    else
        newxgridlowerloc = 1
    end if

    call linearinterpfromlocations(newxgrid, newygrid, newxgridlowerloc, x3val, 2, 1, 1, yval)


    end subroutine linearinterp3
    !---------------------------------------------------------------------------------------------------------!
    !---------------------------------------------------------------------------------------------------------!
    !A subroutine for testing the routine above
    !      subroutine testlinearinterp3
    !        implicit none
    !
    !!Number of test values
    !        integer, parameter :: numtests = 500
    !
    !!How far outside the grid to check for extrapolation		(=1 for only check inside the grid)
    !        real (kind=rk), parameter :: outside = 1.0
    !
    !!Dimensions
    !        integer, parameter :: x1dim = 50
    !        integer, parameter :: x2dim = 51
    !        integer, parameter :: x3dim = 52
    !
    !!Bounds on the grids
    !        real (kind=rk), parameter :: low1 = 0.0, upp1 = 100.0
    !        real (kind=rk), parameter :: low2 = -100.0, upp2 = 300.0
    !        real (kind=rk), parameter :: low3 = -0.5, upp3 = 0.5
    !
    !!Growth rate of the distance between gridpoints
    !        real (kind=rk), parameter :: growth1 = 0.05
    !        real (kind=rk), parameter :: growth2 = 0.00
    !        real (kind=rk), parameter :: growth3 = 0.00
    !
    !!Grids
    !        real (kind=rk) :: x1grid(x1dim)
    !        real (kind=rk) :: x2grid(x2dim)
    !        real (kind=rk) :: x3grid(x3dim)
    !        real (kind=rk) :: ygrid(x1dim, x2dim, x3dim)
    !
    !!Other stuff
    !        real (kind=rk), dimension (numtests, 1) :: x1val, x2val, x3val
    !        real (kind=rk), dimension (numtests) :: interpval, actval, absinterperror, relinterperror
    !        real (kind=rk) :: maxabserror, maxrelerror
    !
    !!Indexes
    !        integer :: ix1, ix2, ix3, ixtest
    !
    !!Fill in the x grids
    !        call getnodes(x1grid, low1, upp1, x1dim, growth1)
    !        call getnodes(x2grid, low2, upp2, x2dim, growth2)
    !        call getnodes(x3grid, low3, upp3, x3dim, growth3)
    !
    !!Fill in the y grids
    !        do ix3 = 1, x3dim, 1
    !          do ix2 = 1, x2dim, 1
    !            do ix1 = 1, x1dim, 1
    !              ygrid(ix1, ix2, ix3) = funcfortestinterp(x1grid(ix1), x2grid(ix2), x3grid(ix3))
    !            end do
    !          end do
    !        end do
    !
    !!Get random numbers
    !        call getdraws(numtests, 1, numtests, 156898, x1val)
    !        x1val = (x1val*(upp1-low1)+low1)*outside
    !
    !        call getdraws(numtests, 1, numtests, 256898, x2val)
    !        x2val = (x2val*(upp2-low2)+low2)*outside
    !
    !        call getdraws(numtests, 1, numtests, 356898, x3val)
    !        x3val = (x3val*(upp3-low3)+low3)*outside
    !
    !!Now do a test: get the exact value and the interpolated value for some points off the grid
    !        do ixtest = 1, numtests, 1
    !          call linearinterp3(x1grid, x2grid, x3grid, x1dim, x2dim, x3dim, x1val(ixtest,1), x2val(ixtest,1), &
    !            x3val(ixtest,1), interpval(ixtest), ygrid)
    !          actval(ixtest) = funcfortestinterp(x1val(ixtest,1), x2val(ixtest,1), x3val(ixtest,1))
    !        end do
    !
    !        absinterperror = abs(interpval-actval)
    !        relinterperror = abs(interpval/actval-1)
    !
    !        maxabserror = maxval(absinterperror)
    !        maxrelerror = maxval(relinterperror)
    !
    !        write (*, *) 'Max absolute interpolation error is', maxabserror
    !        write (*, *) 'Max relative interpolation error is', maxrelerror
    !      contains
    !
    !!The function over which I'll be interpolating
    !        real (kind=rk) function funcfortestinterp(x1, x2, x3)
    !          implicit none
    !          real (kind=rk) :: x1, x2, x3
    !
    !          funcfortestinterp = x1*x2*(x3)
    !
    !        end function funcfortestinterp
    !!---------------------------------------------------------------------------------------------------------!
    !      end subroutine testlinearinterp3
    !---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    subroutine getannuityrates(allyears,survprobs,intrate,annuityrates)

    implicit none

    integer, intent(in) :: allyears
    real (kind=rk), intent(in) :: intrate
    real (kind=rk), intent(in) :: survprobs(allyears)
    real (kind=rk), intent(out) :: annuityrates(allyears)


    integer :: yearsforpension
    integer :: ixforprob, ixforann
    real (kind=rk), allocatable :: probpay(:)

    real (kind=rk) :: recannuityrates(allyears)

    do ixforann = 1, allyears-1, 1
        yearsforpension = allyears - ixforann

        allocate (probpay(yearsforpension))
        do ixforprob = 1, yearsforpension
            !     probpay(ixforprob) = product(survprobs(ixforann + 1:ixforann+ixforprob))
            probpay(ixforprob) = product(survprobs(ixforann:ixforann+ixforprob - 1))

        end do

        call getfairvalue(yearsforpension,probpay,intrate,annuityrates(ixforann))

        deallocate (probpay)
    end do !ixforann

    annuityrates(allyears) = (huge(annuityrates(allyears)))/10.0_rk



    end subroutine getannuityrates
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    subroutine getdiscountfrom1(allyears,beta,discountfrom1)

    implicit none

    integer, intent(in) :: allyears

    real (kind=rk), intent(in) :: beta
    real (kind=rk), intent(out) :: discountfrom1(allyears)

    integer :: ixyear


    discountfrom1(1) = 1.0_rk
    do ixyear = 2, allyears, 1
        discountfrom1(ixyear) = (beta**(ixyear-1))
    end do !ixforprob

    end subroutine getdiscountfrom1
    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    subroutine getfairvalue(maxyears,probpay,intrate,annuityrate)

    implicit none

    integer, intent(in) :: maxyears
    real (kind=rk) :: probpay(maxyears)
    real (kind=rk) :: intrate
    real (kind=rk) :: annuityrate

    !For program
    integer :: ixyear
    real (kind=rk) :: forsum

    forsum = 0.0_rk
    do ixyear = 1, maxyears, 1
        forsum = forsum + probpay(ixyear)/((1+intrate)**ixyear)
    end do

    annuityrate = 1.0_rk/forsum

    end subroutine
    ! ---------------------------------------------------------------------------------------------------------!

    subroutine getannuityratesnew(maxyears,survmale,survfemale,intrate,annuityrate, penpropinherit, annAdminLoad)

    implicit none

    integer, intent(in) :: maxyears
    real (kind=rk), intent(in) :: survmale(maxyears)
    real (kind=rk), intent(in) :: survfemale(maxyears)
    real (kind=rk), intent(in) :: intrate
    real (kind=rk), intent(out) :: annuityrate(maxyears, 3)
    real (kind=rk), intent(in) :: penpropinherit
    real (kind=rk), intent(in) :: annAdminLoad

    real (kind=rk) :: capitalf(maxyears, 3) !This and the following follow the notation in my notes


    !For program
    integer :: ixtoday
    integer :: tau


    capitalf(maxyears, 1) = 0.0_rk
    capitalf(maxyears, 2) = 0.0_rk
    capitalf(maxyears, 3) = 0.0_rk


    !Elsewhere in the programme - I interpret survprob(ixtoday) as the probability as surviving
    !until ixtoday + 1. So, notwithstanding the fact that in my notes I have written this object
    !as s_{ixtoday + 1) I follow the convention here.

    do ixtoday = maxyears-1, 1, -1
        tau = maxyears - ixtoday !years left

        !Couple
        capitalf(ixtoday, 1)= (survmale(ixtoday)*survfemale(ixtoday)*(1.0_rk +  capitalf(ixtoday+1, 1)) &
            + penpropinherit * survmale(ixtoday) * (1.0_rk - survfemale(ixtoday)) * (1.0_rk + capitalf(ixtoday+1, 2)) &
            + penpropinherit * survfemale(ixtoday) * (1.0_rk - survmale(ixtoday)) * (1.0_rk + capitalf(ixtoday+1, 3)) )/(1+intrate)

        annuityrate(ixtoday, 1) = 1.0_rk/( ((1.0_rk)* annAdminLoad) * capitalf(ixtoday, 1))


        !Male
        capitalf(ixtoday, 2)  =   (survmale(ixtoday) * (1.0_rk + capitalf(ixtoday+1, 2)) ) / (1 + intrate)
        annuityrate(ixtoday, 2) = 1.0_rk/( ((1.0_rk)* annAdminLoad) *  capitalf(ixtoday, 2))

        !Female
        capitalf(ixtoday, 3) =   (survfemale(ixtoday) * (1.0_rk + capitalf(ixtoday+1, 3)) )/(1+intrate)
        annuityrate(ixtoday, 3) = 1.0_rk/(((1.0_rk)* annAdminLoad) * capitalf(ixtoday, 3))

    end do

    annuityrate(maxyears, :) = 999.0_rk


    end subroutine getannuityratesnew
    ! ---------------------------------------------------------------------------------------------------------!


    ! ---------------------------------------------------------------------------------------------------------!
    ! ---------------------------------------------------------------------------------------------------------!
    !********************************************************************
    ! FUNCTION: sumindex
    !I took this from Hamish's code
    ! Instructions:
    ! (1) size of nmindex = size of mindex ( = ndim in the code)
    !********************************************************************
    function sumindex(nmindex, mindex)
    implicit none

    !******* input variables *******
    integer, intent (in), dimension (:) :: nmindex, mindex

    !******* output variables *******
    integer :: sumindex

    !******* local variables *******
    integer :: ndim, dim0

    !******* check the dimension of index *******
    ndim = size(nmindex)

    !******* one dimension *******
    if (ndim.eq.1) then
        sumindex = mindex(1)
        return
    end if

    !******* compute the summarized index *******
    sumindex = 0
    do dim0 = 1, ndim - 1
        sumindex = sumindex + (mindex(dim0)-1)*product(nmindex(dim0+1:ndim))
    end do
    sumindex = sumindex + mindex(ndim)

    !******* end of function *******
    return

    end function sumindex

    !********************************************************************
    ! FUNCTION: inv_sumindex
    !I took this from Hamish's code
    ! Instructions:
    ! (1) size of nmindex >= 2 (=ndim in the code)
    ! (2) sindex >= 1 and <= product(nmindex(1:ndim))
    !********************************************************************
    function inv_sumindex(nmindex, sindex)
    implicit none

    !******* input variables *******
    integer, intent (in), dimension (:) :: nmindex
    integer, intent (in) :: sindex

    !******* output variable *******
    integer, dimension (size(nmindex)) :: inv_sumindex

    !******* local variables *******
    integer :: ndim, dim0, indextemp

    !******* check the dimension of index *******
    ndim = size(nmindex)

    !******* one dimension *******
    if (ndim.eq.1) then
        inv_sumindex(1) = sindex
        return
    end if

    !******* set indextemp *******
    indextemp = sindex

    !******* last dimension is a bit different *******
    inv_sumindex(ndim) = mod(indextemp, nmindex(ndim))
    if (inv_sumindex(ndim).eq.0) inv_sumindex(ndim) = nmindex(ndim)
    indextemp = (indextemp-inv_sumindex(ndim))/nmindex(ndim)

    !******* between the 2nd and the last dimension *******
    if (ndim.gt.2) then
        do dim0 = ndim - 1, 2, -1
            inv_sumindex(dim0) = mod(indextemp, nmindex(dim0)) + 1
            indextemp = (indextemp-(inv_sumindex(dim0)-1))/nmindex(dim0)
        end do
    end if

    !******* 1st dimension is a bit different, too *******
    inv_sumindex(1) = indextemp + 1

    !******* end of function *******
    return

    end function inv_sumindex






    !---------------------------------------------------------------------------------------------------------!
    !---------------------------------------------------------------------------------------------------------!
    ! INTEGER FUNCTION  FindMinimum():
    !    This function returns the location of the minimum in the section
    ! between Start and End.
    ! --------------------------------------------------------------------
    !I took this function from (but have adapted)
    ! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
    INTEGER FUNCTION  FindMinimum(x, Start, End)
    IMPLICIT  NONE
    real (kind=rk), DIMENSION(1:), INTENT(IN) :: x
    INTEGER, INTENT(IN)                :: Start, End
    real (kind=rk)                     :: Minimum
    INTEGER                            :: Location
    INTEGER                            :: i

    Minimum  = x(Start)		! assume the first is the min
    Location = Start			! record its position
    DO i = Start+1, End		! start with next elements
        IF (x(i) < Minimum) THEN	!   if x(i) less than the min?
            Minimum  = x(i)		!      Yes, a new minimum found
            Location = i                !      record its position
        END IF
    END DO
    FindMinimum = Location        	! return the position
    END FUNCTION  FindMinimum

    ! --------------------------------------------------------------------
    ! SUBROUTINE  Swap():
    !    This subroutine swaps the values of its two formal arguments.
    ! --------------------------------------------------------------------
    !I took this function from
    ! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90

    SUBROUTINE  Swap(a, b)
    IMPLICIT  NONE
    real (kind=rk), INTENT(INOUT) :: a, b
    real (kind=rk)                :: Temp

    Temp = a
    a    = b
    b    = Temp
    END SUBROUTINE  Swap

    ! --------------------------------------------------------------------
    ! SUBROUTINE  Sort():
    !    This subroutine receives an array x() and sorts it into ascending
    ! order.
    ! --------------------------------------------------------------------
    !I took this function from
    ! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90

    SUBROUTINE  Sort(x, Size)
    IMPLICIT  NONE
    real (kind=rk), DIMENSION(1:), INTENT(INOUT) :: x
    INTEGER, INTENT(IN)                   :: Size
    INTEGER                               :: i
    INTEGER                               :: Location

    DO i = 1, Size-1			! except for the last
        Location = FindMinimum(x, i, Size)	! find min from this to last
        CALL  Swap(x(i), x(Location))	! swap this and the minimum
    END DO
    END SUBROUTINE  Sort
    ! --------------------------------------------------------------------
    ! --------------------------------------------------------------------
    !subroutine getmatrixinverse(matrix, dim, dimplus1, inverse)
    !!Calculate a matrix inverse using the NAG routine. NAG routine gives answer in an odd format - this converts it to a usable format
    !    !use  globalvalues ! only here so that I can access globalrank for an error message
    !
    !   implicit none
    !        integer, intent(in) :: dim
    !        integer, intent(in) :: dimplus1
    !        real (kind=rk), intent(in) :: matrix(dim, dim)
    !        real (kind=rk), intent(out) :: inverse(dim, dim)
    !
    !        !For program
    !        real (kind=rk) :: forinverse(dimplus1, dim)
    !        real (kind=rk) :: copyofinverse(dim, dim)
    !        integer :: ixi, ixj
    !
    !        !For NAG
    !        integer :: ifail
    !
    !        !For checking
    !        real (kind=rk) :: forcheck(dim, dim)
    !        real (kind=rk) :: sumabs
    !
    !        forinverse = 0.0_rk
    !        forinverse(1:dim, :) = matrix
    !
    !       ifail = 1
    !       call F01ADF(dim, forinverse, dimplus1, ifail)
    !       if (ifail.eq.0) then
    !             do ixi = 1, dim
    !               do ixj = 1, ixi
    !                 inverse(ixi, ixj) = forinverse(ixi + 1, ixj)
    !               end do
    !           end do
    !
    !          do ixi = 1, dim
    !               do ixj = ixi, dim
    !                 inverse(ixi, ixj) = inverse(ixj, ixi)
    !               end do
    !          end do
    !
    !
    !          !Check have I got the inverse right
    !          forcheck = matmul(matrix, inverse) !This should be the identity
    !          do ixi = 1, dim
    !            forcheck(ixi, ixi) = forcheck(ixi, ixi) - 1.0_rk
    !          end do
    !
    !          !Get the sum of absolute values
    !          sumabs = 0.0_rk
    !          do ixi = 1, dim
    !              do ixj = 1, dim
    !                sumabs = sumabs + abs(forcheck(ixi, ixj))
    !              end do
    !          end do
    !
    !
    !          if (sumabs.gt.0.001) then
    !                write(*,*) sumabs
    !              write(*,*) 'error in forinverse - inverse either wrong or innacurate'
    !              stop
    !          end if
    !        else !        if (ifail.eq.0) then Inverse didn't computer
    !            inverse = 0
    !            if (globalrank.eq.0) write(*,*) 'error in getmatrixinverse, NAG couldnt compute inverse'
    !        end if
    !
    !
    !  end subroutine getmatrixinverse

    integer function myMinloc(array, dim) ! The Fortan intrinsic has returned an error
    implicit none
    integer, intent(in) :: dim
    real (kind = rk), intent(in) ::  array(dim)

    integer :: ix

    myMinloc = 1
    do ix = 2, dim
        if (array(ix).lt.(array(myMinloc))) myMinloc = ix
    end do

    end function
    !Last-modified: 15 Aug 2011 06:37:29 PM

    SUBROUTINE amoeba(p,y,ftol,func,iter)
    IMPLICIT NONE
    INTEGER(kind=4), INTENT(OUT) :: iter
    REAL(kind=8), INTENT(IN) :: ftol
    REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: y
    REAL(kind=8), DIMENSION(:,:), INTENT(INOUT) :: p
    INTERFACE
    FUNCTION func(x)
    IMPLICIT NONE
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: x
    REAL(kind=8) :: func
    END FUNCTION func
    END INTERFACE
    INTEGER(kind=4), PARAMETER :: ITMAX=5000
    REAL(kind=8), PARAMETER :: TINY=1.0D-10
    INTEGER(kind=4) :: ihi,ndim
    REAL(kind=8), DIMENSION(size(p,2)) :: psum
    call amoeba_private
    CONTAINS

    SUBROUTINE amoeba_private
    IMPLICIT NONE
    INTEGER(kind=4) :: i,ilo,inhi
    REAL(kind=8) :: rtol,ysave,ytry,ytmp
    IF((size(p,2) == size(p,1) - 1) .and. (size(p,2) == size(y) -1))THEN
        ndim = size(y) - 1
    ELSE
        STOP 'ERROR: terminated in amoeba for inconsistent arr dimensions'
    ENDIF
    iter=0
    psum(:)=sum(p(:,:),dim=1)
    do
        ilo=iminloc(y(:))
        ihi=imaxloc(y(:))
        ytmp=y(ihi)
        y(ihi)=y(ilo)
        inhi=imaxloc(y(:))
        y(ihi)=ytmp
        rtol=2.0D0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
        if (rtol < ftol) then
            call swap_scalar(y(1),y(ilo))
            call swap_vector(p(1,:),p(ilo,:))
            RETURN
        end if
        if (iter >= ITMAX) STOP 'ERROR: ITMAX exceeded in amoeba'
        ytry=amotry(-1.0D0)
        iter=iter+1
        if (ytry <= y(ilo)) then
            ytry=amotry(2.0D0)
            iter=iter+1
        else if (ytry >= y(inhi)) then
            ysave=y(ihi)
            ytry=amotry(0.5D0)
            iter=iter+1
            if (ytry >= ysave) then
                p(:,:)=0.5D0*(p(:,:)+spread(p(ilo,:),1,size(p,1)))
                do i=1,ndim+1
                    if (i /= ilo) y(i)=func(p(i,:))
                end do
                iter=iter+ndim
                psum(:)=sum(p(:,:),dim=1)
            end if
        end if
    end do
    END SUBROUTINE amoeba_private

    FUNCTION amotry(fac)
    IMPLICIT NONE
    REAL(kind=8), INTENT(IN) :: fac
    REAL(kind=8) :: amotry
    REAL(kind=8) :: fac1,fac2,ytry
    REAL(kind=8), DIMENSION(size(p,2)) :: ptry
    fac1=(1.0D0-fac)/ndim
    fac2=fac1-fac
    ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
    ytry=func(ptry)
    if (ytry < y(ihi)) then
        y(ihi)=ytry
        psum(:)=psum(:)-p(ihi,:)+ptry(:)
        p(ihi,:)=ptry(:)
    end if
    amotry=ytry
    END FUNCTION amotry

    FUNCTION imaxloc(arr)
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(kind=4) :: imaxloc
    INTEGER(kind=4), DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc=imax(1)
    END FUNCTION imaxloc
    FUNCTION iminloc(arr)
    REAL(kind=8), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(kind=4) :: iminloc
    INTEGER(kind=4), DIMENSION(1) :: imax
    imax=minloc(arr(:))
    iminloc=imax(1)
    END FUNCTION iminloc
    SUBROUTINE swap_scalar(a,b)
    REAL(kind=8), INTENT(INOUT) :: a,b
    REAL(kind=8) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_scalar
    SUBROUTINE swap_vector(a,b)
    REAL(kind=8), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(kind=8), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
    END SUBROUTINE swap_vector
    END SUBROUTINE amoeba
    !---------------------------------------------------------------------------------------------------------!
    !---------------------------------------------------------------------------------------------------------!
    !A subroutine for trilinear interpolation
    !A nice method here (that I'm not using) method at: http://paulbourke.net/miscellaneous/interpolation/
    !This has been tested with the subroutine test below
    subroutine linearinterp3_withextrap(x1grid, x2grid, x3grid, dim1, dim2, dim3, x1val, x2val, x3val, yval, ygrid)
    !use types
    implicit none

    ! Arguments
    integer, intent (in) :: dim1
    integer, intent (in) :: dim2
    integer, intent (in) :: dim3
    real (kind=rk), intent (in) :: x1grid(dim1)
    real (kind=rk), intent (in) :: x2grid(dim2)
    real (kind=rk), intent (in) :: x3grid(dim3)
    real (kind=rk), intent (in) :: ygrid(dim1, dim2, dim3)
    real (kind=rk), intent (in) :: x1val
    real (kind=rk), intent (in) :: x2val
    real (kind=rk), intent (in) :: x3val
    real (kind=rk), intent (out) :: yval

    ! For programme
    integer :: x3lowerloc
    integer :: scratch_nearestloc
    integer :: newxgridlowerloc
    real (kind=rk) :: newygrid(2)
    real (kind=rk) :: newxgrid(2)


    if ((size(x1grid).ne.dim1) .or. (size(x2grid).ne.dim2) .or. (size(x3grid).ne.dim3)) then
        write (*, *) 'An argument in linearinterp3 is not of the specified length'
        stop
    end if

    ! Locate x3val on the grid
    call findingridr(x3grid, dim3, x3val, scratch_nearestloc, x3lowerloc)

    ! Find if any of the values are off the grid
    if (x3lowerloc.eq.0) then
        x3lowerloc = 1
    end if

    if (x3lowerloc.eq.dim3) then
        x3lowerloc = dim3 - 1
    end if

    !Do two two bilinear interpolations over x1grid and x2grid keeping  x3 at x3lowerloc and x3lowerloc + 1
    call linearinterp2_withextrap(x1grid, x2grid, dim1, dim2, x1val, x2val, newygrid(1), 1, 1, ygrid(:,:,x3lowerloc))
    call linearinterp2_withextrap(x1grid, x2grid, dim1, dim2, x1val, x2val, newygrid(2), 1, 1, ygrid(:,:,(x3lowerloc+1)))

    ! Now to a linear interpolation on these found vlues
    ! First generate a new grid from from the x2 bounds (a new ygrid has already bee generated)
    newxgrid(1) = x3grid(x3lowerloc)
    newxgrid(2) = x3grid(x3lowerloc+1)

    !Think this is redundant
    if (x3lowerloc.eq.0) then
        newxgridlowerloc = 0
    else if (x3lowerloc.eq.dim3) then
        newxgridlowerloc = 2
    else
        newxgridlowerloc = 1
    end if

    call linearinterpfromlocations(newxgrid, newygrid, newxgridlowerloc, x3val, 2, 1, 1, yval)


    end subroutine linearinterp3_withextrap

    
    end module routines_generic
