module parse_inverseMod
    !! routine to parse the source table from the input Toml file.
    use constants, only : wp
    use parse_HelpersMod
    use parse_SpectrumMod
    use vector_class

    use tomlf

    implicit none
    
    private
    public :: parse_inverse

contains

    subroutine parse_inverse(table, dict, error)
        !! parse symmetry information, only used in the computation of the escape function
        use sim_state_mod, only : state
        use gridMod,       only : init_grid_cart, init_grid_cyl 
        use vector_class,  only : vector, magnitude
        use tomlf
        use tomlf_error

        !> Input Toml table 
        type(toml_table),              intent(inout) :: table
        !> Error message
        type(toml_error), allocatable, intent(out)   :: error
        !> Dictonary used to store metadata
        type(toml_table),               intent(inout) :: dict

        type(toml_table), pointer :: child
        type(toml_array), pointer :: children

        integer :: maxNumSteps, layer
        real(kind = wp) :: maxStepSize, gradStepSize, accuracy
        logical :: findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing
        character(len=:), allocatable :: outputFile
        real(kind=wp) :: muaUpper, muaLower
        real(kind=wp) :: musUpper, musLower
        real(kind=wp) :: hggUpper, hggLower
        real(kind=wp) :: nUpper, nLower
        real(kind=wp) :: reducedmusUpper, reducedmusLower
        real(kind=wp) :: mueffUpper, mueffLower

        ! Global Optimization Method
        integer :: optimizationMethod

        !used by AdaLIPO
        real(kind=wp) :: alpha

        !used by Bayesian optimizatoin
        real(kind=wp) :: observationNoise
        integer :: trainDataSize
        integer :: fittingDataSize

        !used by AdaLIPO with Bayesian Trust Regieme
        integer :: LIPOtrainingDataSize
        real(kind=wp) :: bayesianMinDist
        real(kind=wp) :: tune
        integer :: numAdaLIPOtoBay

        call get_value(table, "inverse", child)

        if(associated(child))then
            
            call get_value(child, "maxStepSize", maxStepSize, 1.0_wp)
            call set_value(dict, "maxStepSize", maxStepSize)

            call get_value(child, "gradStepSize", gradStepSize, 0.0001_wp)
            call set_value(dict, "gradStepSize", gradStepSize)

            call get_value(child, "accuracy", accuracy, 0.01_wp)
            call set_value(dict, "accuracy", accuracy)

            call get_value(child, "maxNumSteps", maxNumSteps, 1000)
            call set_value(dict, "maxNumSteps", maxNumSteps)
            
            call get_value(child, "Findmua", findmua, .false.)
            call set_value(dict, "Findmua", findmua)

            call get_value(child, "Findmus", findmus, .false.)
            call set_value(dict, "Findmus", findmus)

            call get_value(child, "Findg", findg, .false.)
            call set_value(dict, "Findg", findg)

            call get_value(child, "Findn", findn, .false.)
            call set_value(dict, "Findn", findn)

            call get_value(child, "ReducedmusGuessing", reducedmusGuessing, .false.)
            call set_value(dict, "ReducedmusGuessing", reducedmusGuessing)

            call get_value(child, "MueffGuessing", mueffGuessing, .false.)
            call set_value(dict, "MueffGuessing", mueffGuessing)

            if (reducedmusGuessing) then
                if (.not. (findmus .or. findg)) then
                    call make_error(error, "Must set findmus or findg to be true if using Reducedmus Guessing")
                    return
                end if
            end if

            if (mueffGuessing) then
                if(.not. (findmus .or. findmua .or. findg)) then
                    call make_error(error, "Must set findmus, or findmua, of findg to be true if using mueff Guessing")
                    return
                end if
            end if

            !get bounds on mua
            call get_value(child, "muaUpper", muaUpper, 100.0_wp)
            if (muaUpper < 0.0_wp) then
                call make_error(error, "Must set muaUpper to be greater than 0.0")
                return
            end if
            call set_value(dict, "muaUpper", muaUpper)
            call get_value(child, "muaLower", muaLower, 0.0_wp)
            if(muaLower < 0.0_wp) then
                muaLower = 0.0_wp
            end if
            if(muaLower >= muaUpper) then
                call make_error(error, "Must set muaLower to be less than muaUpper")
                return
            end if
            call set_value(dict, "muaLower", muaLower)

            !get bounds on mus
            call get_value(child, "musUpper", musUpper, 100.0_wp)
            if (musUpper < 0.0_wp) then
                call make_error(error, "Must set musUpper to be greater than 0.0")
                return
            end if
            call set_value(dict, "musUpper", musUpper)
            call get_value(child, "musLower", musLower, 0.0_wp)
            if(musLower < 0.0_wp) then
                musLower = 0.0_wp
            end if
            if(musLower >= musUpper) then
                call make_error(error, "Must set musLower to be less than musUpper")
                return
            end if
            call set_value(dict, "musLower", musLower)

            !get bounds on hgg
            call get_value(child, "hggUpper", hggUpper, 1.0_wp)
            if (hggUpper < -1.0_wp) then
                call make_error(error, "Must set hggUpper to be greater than or equal to -1.0")
                return
            else if (hggUpper > 1.0_wp) then
                call make_error(error, "Must set hggUpper to be less than or equal to 1.0")
                return
            end if
            call set_value(dict, "hggUpper", hggUpper)
            call get_value(child, "hggLower", hggLower, -1.0_wp)
            if (hggLower < -1.0_wp) then
                call make_error(error, "Must set hggLower to be greater than or equal to -1.0")
                return
            else if (hggLower > 1.0_wp) then
                call make_error(error, "Must set hggLower to be less than or equal to 1.0")
            end if
            if(hggLower >= hggUpper) then
                call make_error(error, "Must set hggLower to be less than hggUpper")
                return
            end if
            call set_value(dict, "hggLower", hggLower)

            !get bounds on n
            call get_value(child, "nUpper", nUpper, 20.0_wp)
            if (nUpper < 1.0_wp) then
                call make_error(error, "Must set nUpper to be greater than 1.0")
                return
            end if
            call set_value(dict, "nUpper", nUpper)
            call get_value(child, "nLower", nLower, 1.0_wp)
            if(nLower < 1.0_wp) then
                nLower = 1.0_wp
            end if
            if(nLower >= nUpper) then
                call make_error(error, "Must set nLower to be less than nUpper")
                return
            end if
            call set_value(dict, "nLower", nLower)

            !get bounds on reducedmus
            call get_value(child, "reducedmusUpper", reducedmusUpper, 200.0_wp)
            if (reducedmusUpper < 0.0_wp) then
                call make_error(error, "Must set reducedmusUpper to be greater than 0.0")
                return
            end if
            call set_value(dict, "reducedmusUpper", reducedmusUpper)
            call get_value(child, "reducedmusLower", reducedmusLower, 0.0_wp)
            if(reducedmusLower < 0.0_wp) then
                reducedmusLower = 0.0_wp
            end if
            if(reducedmusLower > reducedmusUpper) then
                call make_error(error, "Must set reducedmusLower to be less than reducedmusUpper")
                return
            end if
            call set_value(dict, "reducedmusLower", reducedmusLower)

            !get bounds on mueff
            call get_value(child, "mueffUpper", mueffUpper, 300.0_wp)
            if (mueffUpper < 0.0_wp) then
                call make_error(error, "Must set mueffUpper to be greater than 0.0")
                return
            end if
            call set_value(dict, "mueffUpper", mueffUpper)
            call get_value(child, "mueffLower", mueffLower, 0.0_wp)
            if(mueffLower < 0.0_wp) then
                mueffLower = 0.0_wp
            end if
            if(mueffLower > mueffUpper) then
                call make_error(error, "Must set mueffLower to be less than mueffUpper")
                return
            end if
            call set_value(dict, "mueffLower", mueffLower)

            call get_value(child, "layer", layer, -985464082)
            if(layer /= -985464082) then
               call set_value(dict, "inverseLayer", layer)
            else 
                call make_error(error, "Must specifiy a layer in inverse table", -1)
                return
            end if

            call get_value(child, "inverseFileName", outputFile, "inverse")
            call set_value(dict, "inverseOutputFileName", outputFile)


        
            call get_value(child, "optimizationMethod", optimizationMethod, 3)
            call set_value(dict, "optimizationMethod", optimizationMethod)           

            call get_value(child, "alpha", alpha, 0.01_wp)
            call set_value(dict, "alpha", alpha)

            call get_value(child, "observationNoise", observationNoise, 0.01_wp)
            call set_value(dict, "observationNoise", observationNoise)

            call get_value(child, "trainDataSize", trainDataSize, 50)
            call set_value(dict, "trainDataSize", trainDataSize)

            call get_value(child, "fittingDataSize", fittingDataSize, 1000)
            call set_value(dict, "fittingDataSize", fittingDataSize)

            call get_value(child, "LIPOtrainingDataSize", LIPOtrainingDataSize, 50)
            call set_value(dict, "LIPOtrainingDataSize", LIPOtrainingDataSize)

            call get_value(child, "bayesianMinDist", bayesianMinDist, 0.2_wp)
            call set_value(dict, "bayesianMinDist", bayesianMinDist)

            call get_value(child, "tune", tune, 0.0_wp)
            call set_value(dict, "tune", tune)

            call get_value(child, "numAdaLIPOtoBay", numAdaLIPOtoBay, 5)
            call get_value(dict, "numAdaLIPOtoBay", numAdaLIPOtoBay)


            if (alpha <= 0.0) then
                call make_error(error, "alpha must be larger than 0.0")
                return
            end if
            if (observationNoise <= 0.0) then
                call make_error(error, "observationNoise must be larger than 0.0")
                return
            end if
            if (trainingDataSize < 0) then
                call make_error(error, "trainingDataSize must be positive")
                return
            end if
            if (fittingDataSize < 1) then
                call make_error(error, "fittingDataSize must be greater than 0")
                return
            end if
            if (LIPOtrainingDataSize < 0) then
                call make_error(error, "LIPOtrainingdataSize must be positive")
                return
            end if
            if (bayesianMinDist <= 0.0) then
                call make_error(error, "bayesianMinDist must be greater than 0.0")
                return
            end if 
            if (tune < 0.0) then
                call make_error(error, "tune must be greater than or equal to 0.0")
                return
            end if
            if (numAdaLIPOtoBay < 0) then
                call make_error(error, "numAdaLIPOtoBay must be greater than or equal to 0")
                return
            end if

        else
            call make_error(error, "Need inverse table in input param file", -1)
            return
        end if

    end subroutine parse_inverse

end module parse_inverseMod