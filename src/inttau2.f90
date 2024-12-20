module inttau2
!! inttau2 is the heart of the MCRT simulation. It moves the photons though the simulated media.
!! tauint2 is the only public function here and is the main function that moves the photon.
!! Changes should only be made here if bugs are discovered or new methods of tracking photons (i.e phase tracking) or moving photons (i.e new geometry method) is needed.

    use constants, only : wp

    implicit none
    
    private
    public :: tauint2, update_voxels

    contains

    subroutine tauint2(grid, packet, sdfs_array, dects, history)
    !! optical depth integration subroutine
    !! Moves photons to interaction location
    !! Calculated is any reflection or refraction happens whilst moving
    !
        use gridMod,      only : cart_grid
        use photonMod,    only : photon
        use random,       only : ran2
        use sdfs,         only : sdf, calcNormal
        use surfaces,     only : reflect_refract
        use vector_class, only : vector
        use detector_mod,  only : hit_t
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t

        type(cart_grid),       intent(in)    :: grid
        type(photon),          intent(inout) :: packet
        type(sdf),             intent(in)    :: sdfs_array(:)
        type(dect_array),      intent(inout) :: dects(:)
        type(history_stack_t), intent(inout) :: history

        real(kind=wp) :: tau, d_sdf, t_sdf, taurun, ds(size(sdfs_array)), dstmp(size(sdfs_array))
        real(kind=wp) :: eps, dtot, old(size(sdfs_array)), new(size(sdfs_array)), n1, n2, Ri
        integer       :: i, oldlayer, new_layer, smallStepLayer
        type(vector)  :: pos, dir, oldpos, N, smallStepPos
        logical       :: rflag

        real(kind=wp) :: pointSep
        type(vector)  :: startPos
        type(hit_t)   :: hpoint

        !setup temp variables
        pos = packet%pos
        oldpos = pos
        startPos = pos
        dir = vector(packet%nxp, packet%nyp, packet%nzp)


        !round off distance
        eps = 1e-8_wp
        !get random tau
        tau = -log(ran2())
        taurun = 0._wp
        dtot = 0._wp
        do while (taurun <= tau)
            !setup sdf distance
            ds = 0._wp
            do i = 1, size(ds)
                ds(i) = sdfs_array(i)%evaluate(pos)
            end do
            packet%cnts = packet%cnts + size(ds)
            d_sdf = minval(abs(ds), dim=1) ! what is the minimum distance to a SDF
            dstmp = ds



            if(d_sdf <eps)then
                ! packet is on a boundary, step backwards or forwards by a small 
                ! amount to stay in the same medium
                
                d_sdf = minval(abs(ds), dim=1) + 2._wp*eps
                smallStepPos = pos + d_sdf*dir
                ds = 0._wp
                do i = 1, size(ds)
                    ds(i) = sdfs_array(i)%evaluate(smallStepPos)
                end do
                packet%cnts = packet%cnts + size(ds)
                smallStepLayer=maxloc(ds,dim=1, mask=(ds<=0._wp))

                !print*, ""
                !print*, packet%layer
                !print*, pos
                !print*, smallStepLayer
                !print*, smallStepPos
                !print*, dir

                
                if (smallStepLayer == packet%layer) then
                    !print*, "forward"
                    !move the packet forwards slightly to stay in the same layer

                    oldpos = pos
                    t_sdf = d_sdf * sdfs_array(packet%layer)%getkappa()

                    if(taurun + t_sdf < tau)then
                    !comment out for phase screen
                        pos = pos + d_sdf*dir
                        taurun = taurun + t_sdf
                        call update_grids(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%getmua())
                    else 
                        d_sdf = (tau - taurun) / sdfs_array(packet%layer)%getkappa()
                        taurun = taurun + t_sdf
                        call update_grids(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%getmua())
                        packet%tflag = .true.
                    end if
                
                else
                    !move the packet backwards slightly to stay in the same layer
                    !print*, "backward"

                    oldpos = pos
                    t_sdf = d_sdf * sdfs_array(packet%layer)%getkappa()
                    
                    !can we move the full amount
                    if(taurun + t_sdf < tau)then
                        pos = pos - d_sdf*dir
                        taurun = taurun + t_sdf
                        !comment out for phase screen
                        call update_grids(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%getmua())
                    else
                        d_sdf = (tau - taurun) / sdfs_array(packet%layer)%getkappa()
                        pos = pos - d_sdf*dir
                        call update_grids(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%getmua())
                        packet%tflag = .true.
                    end if

                end if
            
                !Have we moved through a detector?
                pointSep = sqrt((pos%x - startPos%x)**2 + (pos%y - startPos%y)**2 + (pos%z - startPos%z)**2)
                hpoint = hit_t(startPos, dir, pointSep, sqrt(pos%x**2+pos%y**2), packet%weight)
                do i = 1, size(dects)
                    call dects(i)%p%record_hit(hpoint, history)
                end do
                startPos = pos
                
                !setup sdf distance and current layer
                ds = 0._wp
                do i = 1, size(ds)
                    ds(i) = sdfs_array(i)%evaluate(pos)
                end do
                packet%cnts = packet%cnts + size(ds)
                d_sdf = minval(abs(ds), dim=1)
                dstmp = ds

                !check if outside all sdfs
                if(minval(ds) > -eps)then
                    packet%tflag = .true.
                end if
            end if

            !exit early if conditions met
            if(taurun >= tau .or. packet%tflag)then
                exit
            end if
            
            !print*, ""
            !print*, pos
            !print*, packet%layer
            !print*, sdfs_array(packet%layer)%getkappa()
            !move to the edge or till the packet has moved the full optical depth
            do while(d_sdf >= eps)
                t_sdf = d_sdf * sdfs_array(packet%layer)%getkappa()

                if(taurun + t_sdf < tau)then
                    !move full distance to sdf surface
                    taurun = taurun + t_sdf
                    oldpos = pos
                    !comment out for phase screen
                    call update_grids(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%getmua())
                    pos = pos + d_sdf * dir
                    !dtot = dtot + d_sdf
                else
                    !run out of tau so move remaining tau and exit
                    d_sdf = (tau - taurun) / sdfs_array(packet%layer)%getkappa()
                    dtot = dtot + d_sdf
                    taurun = tau
                    oldpos = pos
                    pos = pos + d_sdf * dir
                    !comment out for phase screen
                    call update_grids(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%getmua())
                    exit
                end if
                ! get distance to nearest sdf
                ds = 0._wp
                do i = 1, size(ds)
                    ds(i) = sdfs_array(i)%evaluate(pos)
                end do
                d_sdf = minval(abs(ds), dim=1)
                packet%cnts = packet%cnts + size(ds)
                dstmp = ds

                !check if outside all sdfs
                if(minval(ds) > -eps)then
                    packet%tflag = .true.
                    exit
                end if
            end do


            ! From the the last startPos has the packet interacted with a detector?
            ! This check will only occur if we need to move the particle across a boundary
            pointSep = sqrt((pos%x - startPos%x)**2 + (pos%y - startPos%y)**2 + (pos%z - startPos%z)**2)
            hpoint = hit_t(startPos, dir, pointSep, sqrt(pos%x**2+pos%y**2), packet%weight)
            do i = 1, size(dects)
                call dects(i)%p%record_hit(hpoint, history)
            end do
            startPos = pos

            !exit early if conditions met
            if(taurun >= tau .or. packet%tflag)then
                exit
            end if

            ! Otherwise the particle must cross a boundary,
            ! we must now consider refraction

            !print*, "across boundary"
            
            !step a small bit into next sdf to get n2
            d_sdf = minval(abs(ds), dim=1) + 2.0_wp*eps
            smallStepPos = pos + d_sdf*dir
            ds = 0._wp
            do i = 1, size(ds)
                ds(i) = sdfs_array(i)%evaluate(smallStepPos)
            end do
            packet%cnts = packet%cnts + size(ds)
            new_layer = maxloc(ds,dim=1, mask=(ds<=0._wp))

            !print*, pos
            !print*, smallStepPos

            !Get n1 and n2
            n1 = sdfs_array(packet%layer)%getn()
            n2 = sdfs_array(new_layer)%getn()           

            !carry out refelction/refraction
            if (n1 /= n2) then !check for fresnel reflection
                print*, "n are different"
                !Need to test and understand
                !get correct sdf normal
                if(ds(packet%layer) < 0._wp .and. ds(new_layer) < 0._wp)then
                    oldlayer = minloc(abs([ds(packet%layer), ds(new_layer)]), dim=1)

                elseif(dstmp(packet%layer) < 0._wp .and. dstmp(new_layer) < 0._wp)then
                    oldlayer = maxloc([dstmp(packet%layer), dstmp(new_layer)], dim=1)

                elseif(ds(packet%layer) > 0._wp .and. ds(new_layer) < 0._wp)then
                    oldlayer = packet%layer

                elseif(ds(packet%layer) > 0._wp .and. ds(new_layer) > 0._wp)then
                    packet%tflag = .true.
                    exit
                else
                    error stop "This should not be reached!"
                end if
                if(oldlayer == 1)then
                    oldlayer = packet%layer
                else
                    oldlayer = new_layer
                end if

                !print*, oldlayer
                !print*, new_layer
                !print*, pos
                !print*, smallStepPos
                N = calcNormal(pos, sdfs_array(oldlayer))
                !print*, N

                rflag = .false.
                call reflect_refract(dir, N, n1, n2, rflag, Ri)
                packet%weight = packet%weight * Ri
                ! Why on earth are we doing this
                !tau = -log(ran2())
                !taurun = 0._wp
                if(.not.rflag)then
                    !update layer and step across the boundary
                    packet%layer = new_layer
                    pos = smallStepPos

                    pointSep = sqrt((pos%x - startPos%x)**2 + (pos%y - startPos%y)**2 + (pos%z - startPos%z)**2)
                    hpoint = hit_t(startPos, dir, pointSep, sqrt(pos%x**2+pos%y**2), packet%weight)
                    do i = 1, size(dects)
                        call dects(i)%p%record_hit(hpoint, history)
                    end do
                    startPos = pos
                else
                    !step back inside original sdf
                    pos = oldpos
                    startPos = oldpos
                    
                    !reflect so incrment bounce counter
                    packet%bounces = packet%bounces + 1
                    ! Also why on earth are we doing this
                    !if(packet%bounces > 1000)then
                    !    packet%tflag=.true.
                    !    exit
                    !end if
                end if
            else
                ! n are equal, no change to the direction
                ! update layer and step across the boundary by a small amount
                packet%layer = new_layer

                oldpos = pos
                !comment out for phase screen
                call update_grids(grid, oldpos, dir, d_sdf, packet, sdfs_array(packet%layer)%getmua())
                t_sdf = d_sdf * sdfs_array(packet%layer)%getkappa()
                taurun = taurun + t_sdf
                pos = smallStepPos

                pointSep = sqrt((pos%x - startPos%x)**2 + (pos%y - startPos%y)**2 + (pos%z - startPos%z)**2)
                hpoint = hit_t(startPos, dir, pointSep, sqrt(pos%x**2+pos%y**2), packet%weight)
                do i = 1, size(dects)
                    call dects(i)%p%record_hit(hpoint, history)
                end do
                startPos = pos

            end if
            if(packet%tflag)exit
        end do

        packet%pos = pos
        packet%nxp = dir%x
        packet%nyp = dir%y
        packet%nzp = dir%z

        packet%phi = atan2(dir%y, dir%x)
        packet%sinp = sin(packet%phi)
        packet%cosp = cos(packet%phi)

        packet%cost = dir%z
        packet%sint = sqrt(1._wp - packet%cost**2)

        ! packet%step = dtot
        if(abs(packet%pos%x) > grid%xmax)then
            packet%tflag = .true.
        end if
        if(abs(packet%pos%y) > grid%ymax)then
            packet%tflag = .true.
        end if
        if(abs(packet%pos%z) > grid%zmax)then
            packet%tflag = .true.
        end if

    end subroutine tauint2


    subroutine update_grids(grid, pos, dir, d_sdf, packet, mua)
    !! record fluence using path length estimators. Uses voxel grid

        use vector_class
        use photonMod
        use gridMod
        use iarray,     only: phasor, jmean, emission, absorb
        use constants , only : sp
        
        !> grid stores voxel grid information (voxel walls and etc)
        type(cart_grid), intent(IN)    :: grid
        !> dir is the current direction (0,0,1) is up
        type(vector),    intent(IN)    :: dir
        !> d_sdf is the distance to travel in voxel grid
        real(kind=wp),   intent(IN)    :: d_sdf
        !> absoprtion coefficent
        real(kind=wp), optional, intent(IN) :: mua
        !> pos is current position with origin in centre of medium (0,0,0)
        type(vector),    intent(INOUT) :: pos
        !> packet stores the photon related variables
        type(photon),    intent(INOUT) :: packet
        
        complex(kind=sp) :: phasec
        type(vector)  :: old_pos
        logical       :: ldir(3)
        integer       :: celli, cellj, cellk
        real(kind=wp) :: dcell, delta=1e-8_wp, d, mua_real

        if(present(mua))then
            mua_real = mua
        else
            mua_real = 1._wp
        end if

        !convert to different coordinate system. Origin is at lower left corner of fluence grid
        old_pos = vector(pos%x+grid%xmax, pos%y+grid%ymax, pos%z+grid%zmax)
        call update_voxels(grid, old_pos, celli, cellj, cellk)
        packet%xcell = celli
        packet%ycell = cellj
        packet%zcell = cellk

        d = 0._wp
        !if packet outside grid return
        if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
            packet%tflag = .true.
            pos = vector(old_pos%x-grid%xmax, old_pos%y-grid%ymax, old_pos%z-grid%zmax)
            return
        end if
        !move photon through grid updating path length estimators
        do
            ldir = (/.FALSE., .FALSE., .FALSE./)

            dcell = wall_dist(grid, celli, cellj, cellk, old_pos, dir, ldir)
            if(d + dcell > d_sdf)then
                dcell = d_sdf - d
                d = d_sdf
! needs to be atomic so dont write to same array address with more than 1 thread at a time
                packet%phase = packet%phase + dcell
!$omp atomic
                    jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + real(dcell, kind=sp)
                call update_pos(grid, old_pos, celli, cellj, cellk, dcell, .false., dir, ldir, delta)
                exit
            else
                d = d + dcell
                packet%phase = packet%phase + dcell
!$omp atomic
                    jmean(celli,cellj,cellk) = jmean(celli,cellj,cellk) + real(dcell, kind=sp)
                call update_pos(grid, old_pos, celli, cellj, cellk, dcell, .true., dir, ldir, delta)
            end if
            if(celli == -1 .or. cellj == -1 .or. cellk == -1)then
                packet%tflag = .true.
                exit
            end if
        end do
        pos = vector(old_pos%x-grid%xmax, old_pos%y-grid%ymax, old_pos%z-grid%zmax)
        packet%xcell = celli
        packet%ycell = cellj
        packet%zcell = cellk

    end subroutine update_grids

    function wall_dist(grid, celli, cellj, cellk, pos, dir, ldir) result(res)
    !! funtion that returns distant to nearest wall and which wall that is (x, y, or z)
    
        use vector_class
        use gridMod

        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: pos, dir
        logical,         intent(INOUT) :: ldir(:)
        integer,         intent(INOUT) :: celli, cellj, cellk
        real(kind=wp) :: res

        real(kind=wp) :: dx, dy, dz

        dx = -999._wp
        dy = -999._wp
        dz = -999._wp

        if(dir%x > 0._wp)then
            dx = (grid%xface(celli+1) - pos%x)/dir%x
        elseif(dir%x < 0._wp)then
            dx = (grid%xface(celli) - pos%x)/dir%x
        elseif(dir%x == 0._wp)then
            dx = 100000._wp
        end if

        if(dir%y > 0._wp)then
            dy = (grid%yface(cellj+1) - pos%y)/dir%y
        elseif(dir%y < 0._wp)then
            dy = (grid%yface(cellj) - pos%y)/dir%y
        elseif(dir%y == 0._wp)then
            dy = 100000._wp
        end if

        if(dir%z > 0._wp)then
            dz = (grid%zface(cellk+1) - pos%z)/dir%z
        elseif(dir%z < 0._wp)then
            dz = (grid%zface(cellk) - pos%z)/dir%z
        elseif(dir%z == 0._wp)then
            dz = 100000._wp
        end if

        res = min(dx, dy, dz)
        if(res < 0._wp)then
            print*,'dcell < 0.0 warning! ',res
            print*,dx,dy,dz
            print*,dir
            print*,celli,cellj,cellk
            error stop 1
        end if

        ldir = [res == dx, res==dy, res==dz]
        if(.not.ldir(1) .and. .not.ldir(2) .and. .not.ldir(3))print*,'Error in dir flag'
      
   end function wall_dist


    subroutine update_pos(grid, pos, celli, cellj, cellk, dcell, wall_flag, dir, ldir, delta)
    !! routine that updates positions of photon and calls Fresnel routines if photon leaves current voxel
    
        use vector_class
        use gridMod
        use utils, only : str
      
        type(cart_grid), intent(IN)    :: grid
        type(vector),    intent(IN)    :: dir
        logical,         intent(IN)    :: wall_flag, ldir(:)
        real(kind=wp),   intent(IN)    :: dcell, delta
        type(vector),    intent(INOUT) :: pos
        integer,         intent(INOUT) :: celli, cellj, cellk

        if(wall_flag)then

            if(ldir(1))then
                if(dir%x > 0._wp)then
                    pos%x = grid%xface(celli+1) + delta
                elseif(dir%x < 0._wp)then
                    pos%x = grid%xface(celli) - delta
                else
                    print*,'Error in x ldir in update_pos', ldir, dir
                end if
                pos%y = pos%y + dir%y*dcell 
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(2))then
                if(dir%y > 0._wp)then
                    pos%y = grid%yface(cellj+1) + delta
                elseif(dir%y < 0._wp)then
                    pos%y = grid%yface(cellj) - delta
                else
                    print*,'Error in y ldir in update_pos', ldir, dir
                end if
                pos%x = pos%x + dir%x*dcell
                pos%z = pos%z + dir%z*dcell
            elseif(ldir(3))then
                if(dir%z > 0._wp)then
                    pos%z = grid%zface(cellk+1) + delta
                elseif(dir%z < 0._wp)then
                    pos%z = grid%zface(cellk) - delta
                else
                    print*,'Error in z ldir in update_pos', ldir, dir
                end if
                pos%x = pos%x + dir%x*dcell
                pos%y = pos%y + dir%y*dcell 
            else
                print*,'Error in update_pos... '//str(ldir)
                error stop 1
            end if
        else
            pos%x = pos%x + dir%x*dcell
            pos%y = pos%y + dir%y*dcell 
            pos%z = pos%z + dir%z*dcell
        end if

        if(wall_flag)then
            call update_voxels(grid, pos, celli, cellj, cellk)
        end if

    end subroutine update_pos


    subroutine update_voxels(grid, pos, celli, cellj, cellk)
    !! updates the current voxel based upon position

        use vector_class
        use gridmod
        
        !> grid
        type(cart_grid), intent(IN)    :: grid
        !> current photon packet position
        type(vector),    intent(IN)    :: pos
        !> position of photon packet in grid
        integer,         intent(INOUT) :: celli, cellj, cellk

        !accurate but slow
        ! celli = find(pos%x, grid%xface) 
        ! cellj = find(pos%y, grid%yface)
        ! cellk = find(pos%z, grid%zface) 

        !fast but can be inaccurate in some cases...
        celli = floor(grid%nxg * (pos%x) / (2. * grid%xmax)) + 1
        cellj = floor(grid%nyg * (pos%y) / (2. * grid%ymax)) + 1
        cellk = floor(grid%nzg * (pos%z) / (2. * grid%zmax)) + 1

        if(celli > grid%nxg .or. celli < 1)celli = -1
        if(cellj > grid%nyg .or. cellj < 1)cellj = -1
        if(cellk > grid%nzg .or. cellk < 1)cellk = -1

    end subroutine update_voxels

    integer function find(val, a)
    !! searches for bracketing indices for a value value in an array a

        !> value to find in array
        real(kind=wp), intent(in) :: val
        !> array to find val in
        real(kind=wp), intent(in) :: a(:)
        
        integer :: n, lo, mid, hi

        n = size(a)
        lo = 0
        hi = n + 1

        if (val == a(1)) then
            find = 1
        else if (val == a(n)) then
            find = n-1
        else if((val > a(n)) .or. (val < a(1))) then
            find = -1
        else
            do
                if (hi-lo <= 1) exit
                mid = (hi+lo)/2
                if (val >= a(mid)) then
                    lo = mid
                else
                    hi = mid
                end if
            end do
            find = lo
        end if
    end function find
end module inttau2