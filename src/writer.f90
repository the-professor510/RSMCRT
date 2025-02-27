module writer_mod
!! This module defines all functions that write simulation data to the disk or pre-process data before writing.
!! normalise_fluence. Normalises fluence by number of photons run and size of each voxel. **!Does not normalise by power!**
!! write_fluence. Write out fluence in either raw or nrrd format. Default is nrrd.
!! write_detected_photons. Write out photons detected by detectors.

!! Changes should only be made here if there is a bug or new data types need to be written to disk (phase information) or new file format is needed.

    use constants, only : wp

    implicit none

    interface nrrd_write
        module procedure write_3d_r8_nrrd, write_3d_r4_nrrd
    end interface nrrd_write

    interface raw_write
        module procedure write_3d_r8_raw, write_3d_r4_raw
    end interface raw_write

    private
    public :: normalise_fluence, write_data, write_detected_photons, checkpoint, write_escape

    contains
        subroutine normalise_fluence(grid, array, nphotons)
        !! normalise fluence in the Lucy 1999 way
            
            use gridMod
            use constants, only : sp

            !> grid class
            type(cart_grid), intent(in) :: grid
            !> array to normalise
            real(kind=sp),   intent(inout) :: array(:, :, :)
            !> number of photons run
            integer,         intent(in) :: nphotons
            
            real(kind=wp) :: xmax, ymax, zmax
            integer       :: nxg, nyg, nzg

            nxg = grid%nxg
            nyg = grid%nyg
            nzg = grid%nzg
            xmax = grid%xmax
            ymax = grid%ymax
            zmax = grid%zmax

            array  = array * ((2._sp*xmax*2._sp*ymax*2._sp*zmax)&
                            /(nphotons * (2._sp * xmax / nxg) * &
                            (2._sp * ymax / nyg) * (2._sp * zmax / nzg)))

        end subroutine normalise_fluence


        subroutine write_detected_photons(dects)

            use detectors
            use constants, only: fileplace, wp
            use utils, only : str
            use sim_state_mod, only : state

            type(dect_array), intent(in) :: dects(:)

            integer :: i, j, u, test
            character(len=:), allocatable :: hdr

            do i = 1, size(dects)
                
                open(newunit=u,file=trim(fileplace)//"detectors/detector_"//str(i)//".dat",&
                    access='stream',status='REPLACE',form='unformatted')
                    associate(x => dects(i)%p)
                    select type(x)
                    !write out different data for different detector types
                    type is(circle_dect)
                        ! write out the circle detector 
                        write(u)  1.0_wp ! What type of detector is it
                        write(u)  real(len(x%ID), kind=wp) !write out the detector ID
                        do j = 1, len(x%ID) 
                            write(u)  real(ichar(x%ID(j:j)), kind=wp)
                        end do
                        write(u)  real(state%nphotons, kind=wp)
                        write(u)  x%radius
                        write(u)  x%pos
                        write(u)  x%dir
                        do j = 1, x%nbins
                            write(u)(real(j,kind=wp)-0.5_wp) * x%bin_wid, x%data(j)
                        end do
                    type is(fibre_dect)
                        ! write out the fibre detector
                        write(u)  2.0_wp !What type of detector is it
                        write(u)  real(len(x%ID), kind=wp) !write out the detector ID
                        do j = 1, len(x%ID) 
                            write(u)  real(ichar(x%ID(j:j)), kind=wp)
                        end do
                        write(u)  real(state%nphotons, kind=wp)
                        write(u)  x%pos
                        write(u)  x%dir
                        write(u)  x%focalLength1
                        write(u)  x%focalLength2
                        write(u)  x%f1Aperture
                        write(u)  x%f2Aperture
                        write(u)  x%frontOffset
                        write(u)  x%backOffset
                        write(u)  x%frontToPinSep
                        write(u)  x%pinToBackSep
                        write(u)  x%pinAperture
                        write(u)  x%acceptAngle
                        write(u)  x%coreDiameter
                        do j = 1, x%nbins
                            write(u)(real(j,kind=wp)-0.5_wp) * x%bin_wid, x%data(j)
                        end do
                    type is(annulus_dect)
                        ! write out the annulus detector
                        write(u)  3.0_wp ! What type of detector is it
                        write(u)  real(len(x%ID), kind=wp) !write out the detector ID
                        do j = 1, len(x%ID) 
                            write(u)  real(ichar(x%ID(j:j)), kind=wp)
                        end do
                        write(u)  real(state%nphotons, kind=wp)
                        write(u)  x%r1
                        write(u)  x%r2
                        write(u)  x%pos
                        write(u)  x%dir
                        do j = 1, x%nbins
                            write(u)((real(j,kind=wp)-0.5_wp) * x%bin_wid + x%r1), x%data(j)
                        end do
                    type is(camera)
                        print*,"Warning camera detector not yet implmented!"
                    end select
                    end associate
                close(u)
            end do

        end subroutine write_detected_photons

        subroutine write_escape(dects, dict, overwrite)

            use iarray
            use detectors
            use constants, only: fileplace, wp
            use utils, only : str
            use sim_state_mod, only : state
            use tomlf,         only : toml_table


            !> list of detectors
            type(dect_array), intent(in) :: dects(:)
            !> dictionary of metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> overwrite flag
            logical,          optional, intent(IN)    :: overwrite

            !> filename to save array as
            character(len=:), allocatable    :: filename

            integer :: i, j, u
            character(len=:), allocatable :: hdr

            do i = 1, size(dects)
                filename = trim(fileplace)//"escape/dectID_"//trim(dects(i)%p%ID)//"__escape"//trim(str(i))//".nrrd" 
                call write_data(escape(i,:,:,:), filename, state, dict, overwrite, dects(i)%p%ID)

                filename = trim(fileplace)//"escape/dectID_"//trim(dects(i)%p%ID)//"__escapeSym"//trim(str(i))//".nrrd" 
                call write_data(escapeSymmetry(i,:,:,:), filename, state, dict, overwrite, dects(i)%p%ID)
            end do
        end subroutine write_escape


        subroutine write_data(array, filename, state, dict, overwrite, dect_ID)
        !! routine automatically selects which way to write out results based upon file extension
            
            use sim_state_mod, only : settings_t
            use tomlf,         only : toml_table, get_value
            use constants,     only : sp

            !> simulation state
            type(settings_t),           intent(IN)    :: state
            !> array to write out
            real(kind=sp),              intent(IN)    :: array(:,:,:)
            !> filename to save array as
            character(*),               intent(IN)    :: filename
            !> dictionary of metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> overwrite flag
            logical,          optional, intent(IN)    :: overwrite
            !> detector ID, only used writing the escape function
            character(len=:), allocatable, optional, intent(in) :: dect_ID

            Logical :: over_write
            integer :: pos
            
            if(present(overwrite))then
                over_write = overwrite
            else
                over_write = state%overwrite
            end if

            pos = index(filename, ".nrrd")
            if(pos > 0)then
                if(present(dict))then
                    call nrrd_write(array, filename, over_write, dict, dect_ID = dect_ID)
                else
                    call nrrd_write(array, filename, over_write, dect_ID = dect_ID)
                end if
                return
            end if

            pos = index(filename, ".raw")
            if(pos > 0)then
                call raw_write(array, filename, over_write)
                return
            end if

            pos = index(filename, ".dat")
            if(pos > 0)then
                call raw_write(array, filename, over_write)
                return
            end if

            error stop "File type not supported!"

        end subroutine write_data

        subroutine write_3d_r8_raw(array, filename, overwrite)
        !! write 3D array of float64s to disk as raw binary data
            
            !> array to write to disk
            real(kind=wp), intent(IN) :: array(:, :, :)
            !> filename to save array as
            character(*),  intent(IN) :: filename
            !> overwrite flag
            logical,       intent(IN) :: overwrite

            integer :: u
            character(len=:), allocatable :: file

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if
            open(newunit=u,file=file,access='stream',status='REPLACE',form='unformatted')
            write(u) array
            close(u)

        end subroutine write_3d_r8_raw

        subroutine write_3d_r4_raw(array, filename, overwrite)
        !! write 3D array of float32's to disk as raw binary data
            use constants, only : sp

            !> array to write to disk
            real(kind=sp), intent(IN) :: array(:, :, :)
            !> filename to save array as
            character(*),  intent(IN) :: filename
            !> overwrite flag
            logical,       intent(IN) :: overwrite

            integer :: u
            character(len=:), allocatable :: file

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if
            open(newunit=u,file=file,access='stream',status='REPLACE',form='unformatted')
            write(u) array
            close(u)

        end subroutine write_3d_r4_raw

        function get_new_file_name(file) result(res)
            !! If file exits, get numeral to append to filename 

            use utils, only : str

            !> file to be checked
            character(len=*), intent(IN) :: file

            character(len=:), allocatable :: res
            integer :: pos, i

            i = 1
            do
                pos = scan(trim(file), ".", back=.true.)
                res = file(1:pos-1)//" ("//str(i)//")"//file(pos:)
                if(.not. check_file(res))exit
                i = i + 1
            end do

        end function get_new_file_name

        logical function check_file(file) result(res)
            !! Functional wrapper around inquire to check if file exits

            !> file to be checked
            character(len=*), intent(IN) :: file

            inquire(file=trim(file), exist=res)
        
        end function check_file

        subroutine write_hdr(u, sizes, type, dect_ID)
            !! write out header information for .nrrd file format
            use utils, only : str

            !> data dtype
            character(*), intent(IN) :: type
            !> file handle
            integer,      intent(IN) :: u
            !> dimensions of data
            integer,      intent(IN) :: sizes(:)
            !> detector ID, only used writing the escape function
            character(len=:), allocatable, optional, intent(in) :: dect_ID
            
            character(len=100) :: string
            integer :: i

            string = ""
            string = str(sizes(3))
            string = trim(string) // " " // str(sizes(2))
            string = trim(string) // " " // str(sizes(1))

            write(u,"(A)")"NRRD0004"
            write(u,"(A)")"type: "//type
            write(u,"(A)")"dimension: "//str(size(sizes))
            write(u,"(A)")"sizes: "//trim(string)
            write(u,"(A)")"space dimension: "//str(size(sizes))
            write(u,"(A)")"encoding: raw"
            write(u,"(A)")"endian: little"

            if(present(dect_ID)) then
                write(u, "(A)")"dector: "//(trim(dect_ID))
            end if

        end subroutine write_hdr

        subroutine write_3d_r8_nrrd(array, filename, overwrite, dict, dect_ID)
            !! write 3D array of float64's to .nrrd fileformat

            use tomlf,           only : toml_table, toml_dump, toml_error
            use iso_fortran_env, only : int32, int64, real32, real64
            use utils,    only : str
            
            !> filename
            character(*),               intent(IN)    :: filename
            !> array to be written to disk
            real(kind=wp),              intent(IN)    :: array(:, :, :)
            !> dictionary of metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> overwrite flag
            logical,                    intent(IN)    :: overwrite
            !> detector ID, only used writing the escape function
            character(len=:), allocatable, optional, intent(in) :: dect_ID

            type(toml_error), allocatable :: error
            character(len=:), allocatable :: file
            integer :: u

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file,form="formatted")
            !to do fix precision
            call write_hdr(u, [size(array, 1), size(array, 2), size(array, 3)], "double", dect_ID)

            if(present(dict))then
                call toml_dump(dict, u, error)
            end if
            write(u,"(A)")new_line("C")
            close(u)
            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)array
            close(u)
        
        end subroutine write_3d_r8_nrrd

        subroutine write_3d_r4_nrrd(array, filename, overwrite, dict, dect_ID)
            !! write 3D array of float32's to .nrrd fileformat

            use tomlf,           only : toml_table, toml_dump, toml_error
            use iso_fortran_env, only : int32, int64, real32, real64
            use utils,           only : str
            use constants,       only : sp
            
            !> filename
            character(*),               intent(IN)    :: filename
            !> array to be written to disk
            real(kind=sp),              intent(IN)    :: array(:, :, :)
            !> dictionary of metadata
            type(toml_table), optional, intent(INOUT) :: dict
            !> overwrite flag
            logical,                    intent(IN)    :: overwrite
            !> detector ID, only used writing the escape function
            character(len=:), allocatable, optional, intent(in) :: dect_ID

            type(toml_error), allocatable :: error
            character(len=:), allocatable :: file
            integer :: u

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file,form="formatted")
            !to do fix precision
            call write_hdr(u, [size(array, 1), size(array, 2), size(array, 3)], "float", dect_ID)

            if(present(dict))then
                call toml_dump(dict, u, error)
            end if
            write(u,"(A)")new_line("C")
            close(u)
            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)array
            close(u)
        
        end subroutine write_3d_r4_nrrd

        subroutine checkpoint(toml_filename, filename, nphotons_run, overwrite)

            use iarray, only : jmean

            !> filename of toml file used in simulation
            character(*), intent(IN)    :: toml_filename
            !> name of checkpoint file to be saved
            character(*), intent(IN)    :: filename
            !> flag which determines if file is to be overwritten or adjusted
            logical,      intent(IN) :: overwrite
            !> number of photons run up to checkpoint
            integer,      intent(IN) :: nphotons_run

            character(len=:), allocatable :: file
            integer :: u

            if(check_file(filename) .and. .not. overwrite)then
                file = get_new_file_name(filename)
            else
                file = filename
            end if

            open(newunit=u,file=file)
            write(u,"(a,a)")"tomlfile=",toml_filename
            write(u,"(a,i0)")"photons_run=",nphotons_run
            close(u)

            open(newunit=u,file=file,access="stream",form="unformatted",position="append")
            write(u)jmean
            close(u)

        end subroutine checkpoint
end module writer_mod
