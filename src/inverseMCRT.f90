module inverseMCRTMod
!! Contains the main program and scattering loop. Calls all other routine to setup, run and break down the simulation.

    implicit none
    
    private
    public :: inverse_MCRT  
    
    
contains  
    subroutine inverse_MCRT(input_file)
        !Shared data
        use iarray
        use constants, only : wp, fileplace, sp, TWOPI

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use opticalProperties, only : opticalProp_t, mono
        use setupMod, only : setup_inverseDirectory
        use writer_mod,    only : write_inverse
        use kernels, only : setup, finalise, reset_detectors

        !interface for stdlib
        use Interfaces, only : sposv

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, get_value
#ifdef _OPENMP
        use omp_lib
#endif

        character(len=*), intent(in) :: input_file
        
        integer                       :: i, j, loopCounter
        type(history_stack_t)         :: history
        type(photon)                  :: packet
        type(toml_table)              :: dict
        real(kind=wp),    allocatable :: distances(:), image(:,:,:)
        type(dect_array), allocatable :: dects(:)
        type(sdf),        allocatable :: array(:)
        real(kind=wp)                 :: nscatt, start
        type(spectrum_t)              :: spectrum
        type(tevipc)                  :: tev

        integer :: nphotons_run,pos
        character(len=128) :: line
        character(len=:), allocatable :: checkpt_input_file
        character(len=:), allocatable :: outputFile
        integer :: maxNumSteps, layer, numGuesses
        real(kind = wp) :: accuracy
        logical :: findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing
        real(kind=wp) :: error
        real(kind=wp), allocatable :: optimizationData(:,:)
        real(kind=wp) :: Tmus, Tmua, Thgg, Tn, Treducedmus, Tmueff !Accpeted true values of mus, mua, hgg, and n
        real(kind=wp) :: mus, mua, hgg, n, reducedmus, mueff !Current guess values
        real(kind=wp) :: domain(6,2) !domain of mus, mua, hgg, n, reducedmus, and museff
        integer :: SDF_array_index


        !optical properties of the simulation
        integer :: NoVariablesToOptimize
        type(opticalProp_t) :: trialOptProp
        real(kind=wp) :: temp !temp variable for changing the optical properties of the simulation


        !random seed for choosing random optical properties
        integer, allocatable :: seed(:)
        integer :: sizeRanSeed
        real(kind=wp) :: ranNum

        integer :: inverseMethod

        !used by AdaLIPO_Method
        real(kind=wp) ::  alpha

        !used by Bayesian_Method
        real(kind=wp) :: observationNoise
        integer :: trainingDataSize
        integer :: fittingDataSize

        !used by combined
        integer :: LIPOtrainingDataSize
        real(kind=wp) :: tune
        real(kind=wp) :: bayesianMinDist
        integer :: numAdaLIPOtoBay


        !setup the simulation
        if(state%loadckpt)then
            call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .false.)
            open(newunit=j,file=state%ckptfile, access="stream", form="formatted")
            read(j,"(a)")line
            pos = scan(line, "=")
            checkpt_input_file = trim(line(pos+1:))

            read(j,"(a)")line
            pos = scan(line, "=")
            read(line(pos+1:),*) nphotons_run

            inquire(j,pos=pos)
            close(j)

            open(newunit=j,file=state%ckptfile, access="stream", form="unformatted")
            read(j,pos=pos)jmean
            close(j)

            call setup(checkpt_input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .true.)
            state%iseed=state%iseed*101
            state%nphotons = state%nphotons - nphotons_run
        else
            call setup(input_file, tev, dects, array, packet, spectrum, dict, distances, image, nscatt, start, .true.)
        end if

        !read in the data used in AdaLIPO
        call get_value(dict, "accuracy", accuracy)
        call get_value(dict, "maxNumSteps", maxNumSteps)
        call get_value(dict, "Findmua", findmua)
        call get_value(dict, "Findmus", findmus)
        call get_value(dict, "Findg", findg)
        call get_value(dict, "Findn", findn)
        call get_value(dict, "ReducedmusGuessing", reducedmusGuessing)
        call get_value(dict, "MueffGuessing", mueffGuessing)
        call get_value(dict, "inverseLayer", layer)
        call get_value(dict, "inverseOutputFileName", outputFile)

        
        !check if the inverse MCRT directory exists
        call setup_inverseDirectory()

        !setup different random number seed for the 
        call random_seed(size=sizeRanSeed)
        allocate(seed(sizeRanSeed))
        seed = 0
        seed = state%iseed
        call random_seed(put=seed)
        call random_seed(get=seed)
        !discard first 100 to ensure that we get random numbers
        do i = 1, 100
            call random_number(ranNum)
        end do
        call random_seed(get=seed)


        !check how many variables we will optimize/search for
        NoVariablesToOptimize = 0
        if (findmua) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if (findmus) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if (findg) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if (findn) then
            NoVariablesToOptimize = NoVariablesToOptimize + 1
        end if
        if(NoVariablesToOptimize == 0) then
            print*, "Please select at least one of mus, mua, hgg, n to find with inverse MCRT"
            return 
        end if

        !loop through the layers finding the index in the SDF array of the selected layer and its initial optical properties
        SDF_array_index = -1
        do i = 1, size(array)
            if (array(i)%getLayer() == layer) then 
                !we have found the layer, store the index of this and its optical properties
                SDF_array_index = i
                Tmua = array(i)%getMua()
                Tmus = array(i)%getKappa() - Tmua
                Thgg = array(i)%gethgg()
                Tn = array(i)%getN()
                exit
            end if
        end do

        !check that the selected layer is found in the SDF array
        if (SDF_array_index == -1) then
            print*, "Selected layer not found in SDF array please choose a layer inside the SDF array"
            return
        end if

        !bounds for AdaLIPO
        call get_value(dict, "musLower", domain(1,1))
        call get_value(dict, "musUpper", domain(1,2))
        call get_value(dict, "muaLower", domain(2,1))
        call get_value(dict, "muaUpper", domain(2,2))
        call get_value(dict, "hggLower", domain(3,1))
        call get_value(dict, "hggUpper", domain(3,2))
        call get_value(dict, "nLower", domain(4,1))
        call get_value(dict, "nUpper", domain(4,2))
        call get_value(dict, "reducedmusLower", domain(5,1))
        call get_value(dict, "reducedmusUpper", domain(5,2))        
        call get_value(dict, "mueffLower", domain(6,1))
        call get_value(dict, "mueffUpper", domain(6,2))

        !set the initial guesses
        allocate(optimizationData(maxNumSteps, 5))


        !set the values for AdaLIPO
        alpha = 0.01

        !set the values for Bayesian
        observationNoise = 0.01
        trainingDataSize = 50
        fittingDataSize = 1000

        !set the values for AdaLIPO with Bayesian Trust Regieme
        alpha = 0.01
        observationNoise = 0.01
        fittingDataSize = 1000
        LIPOtrainingDataSize = 50
        bayesianMinDist = 0.2_wp
        tune = 0.0
        numAdaLIPOtoBay = (5) + 1 !(user value) + 1 for modulus to work and make sense


        !get inverseMethod from the input file
        ! 1 = AdaLIPO only
        ! 2 = Bayesian only
        ! 3 = AdaLIPO with Bayesian trust region
        inverseMethod = 3
        
        if (inverseMethod == 1) then
            call AdaLIPO_Method(input_file, history, packet, dict, distances, image, dects, array, nscatt, start, tev, &
                                spectrum, seed, NoVariablesToOptimize, trialOptProp, SDF_array_index, optimizationData, &
                                Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, &
                                domain, outputFile, layer, accuracy, maxNumSteps, alpha)
                                
        else if (inverseMethod == 2) then
            call Bayesian_Method(input_file, history, packet, dict, distances, image, dects, array, nscatt, start, tev, &
                                spectrum, seed, NoVariablesToOptimize, trialOptProp, SDF_array_index, optimizationData, &
                                Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, &
                                domain, outputFile, layer, accuracy, maxNumSteps, trainingDataSize, fittingDataSize, &
                                observationNoise)

        else if (inverseMethod == 3) then
            call AdaLIPO_with_BayesianTrust_Method(input_file, history, packet, dict, distances, image, dects, array, nscatt, &
                                start, tev, spectrum, seed, NoVariablesToOptimize, trialOptProp, SDF_array_index, &
                                optimizationData, Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, reducedmusGuessing, &
                                mueffGuessing, domain, outputFile, layer, accuracy, maxNumSteps, alpha, fittingDataSize, &
                                observationNoise, LIPOtrainingDataSize, bayesianMinDist, tune, numAdaLIPOtoBay)
        end if       

    end subroutine inverse_MCRT


    subroutine AdaLIPO_Method(input_file, history, packet, dict, distances, image, dects, array, nscatt, start, tev, spectrum, &
        seed, NoVariablesToOptimize, trialOptProp, SDF_array_index, optimizationData, Tmus, Tmua, Thgg, Tn, findmua, findmus, &
        findg, findn, reducedmusGuessing, mueffGuessing, domain, outputFile, layer, accuracy, maxNumSteps, alpha)
        !AdaLIPO algorithm for global optimization

        !Shared data
        use iarray
        use constants, only : wp, fileplace

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use opticalProperties, only : opticalProp_t, mono
        use setupMod, only : setup_inverseDirectory
        use writer_mod,    only : write_inverse
        use default_MCRTMod, only : run_MCRT
        use kernels, only : setup, finalise, reset_detectors

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, get_value
#ifdef _OPENMP
        use omp_lib
#endif

        character(len=*),              intent(in) :: input_file
        type(history_stack_t),         intent(inout) :: history
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt, start
        type(spectrum_t),              intent(inout) :: spectrum
        type(tevipc),                  intent(inout) :: tev

        integer, allocatable, intent(inout) :: seed(:)
        integer, intent(in) :: SDF_array_index
        integer, intent(in) :: NoVariablesToOptimize
        type(opticalProp_t), intent(inout) :: trialOptProp


        real(kind=wp), intent(inout) :: optimizationData(:,:)
        real(kind=wp), intent(in) :: Tmus, Tmua, Thgg, Tn
        logical, intent(in) :: findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing
        real(kind=wp), intent(in) :: domain(:,:)
        character(len=:), allocatable, intent(inout) :: outputFile
        integer, intent(in) :: maxNumSteps
        integer, intent(in) :: layer
        real(kind=wp), intent(in) :: accuracy

        !unique to AdaLIPO_Method
        real(kind=wp), intent(in) :: alpha


        integer :: i, j
        real(kind=wp) :: mus, mua, hgg, n, reducedmus, mueff !Current guess values
        real(kind=wp) :: temp !temp variable for changing the optical properties of the simulation
        real(kind=wp) :: error
        integer :: numGuesses


        !used by AdaLIPO
        real(kind=wp) :: probability, k, it
        integer :: indexOfMinError, indexOfMaxRatio, index, ratioCounter
        real(kind=wp) :: minError, maxRatio
        real(kind=wp), allocatable :: ratios(:)
        real(kind=wp) :: ranNum !used to temporarily store a random number
        real(kind=wp) :: leftMin, rightMax, tempMin !sides of the LIPO condition


        !choose random mus, mua, hgg, n
        call randomOptProp(optimizationData(1,1), optimizationData(1,2), optimizationData(1,3), optimizationData(1,4), & 
                            reducedmus, mueff, Tmus, Tmua, Thgg, Tn, &
                            findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, domain)

        trialOptProp = mono(optimizationData(1,1), optimizationData(1,2), optimizationData(1,3), optimizationData(1,4))
        temp = array(SDF_array_index)%updateOptProp(trialOptProp)

        !evaluate, by running MCRT
        call random_seed(get=seed) !store the random seed for optical properties
        call run_MCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)
        call random_seed(put=seed) !restart the random seed for optical properties
        error = 0._wp
        call inverse_evaluate(dects, error)
        optimizationData(1,5) = error

        ! reset the arrays storing data
        call reset_detectors(dects)  

        !store the position of minimum error
        minError = optimizationData(1,5)
        indexOfMinError = 1

        !store the maximum value
        rightMax = optimizationData(1,5)

        !allocate the ratio array (it will be of the size maxNumSteps triangular number)
        allocate(ratios(int(nint(maxNumSteps * (maxNumSteps+1)/2.0_wp))))
        ratios = 0.0_wp
        ratioCounter = 1
        maxRatio = 0.0_wp
        indexOfMaxRatio = 1
        k = 0.0_wp
        probability = 1.0_wp

        do i = 2, maxNumSteps

            
            ranNum = ran2()
            if( ranNum <= probability) then
                !we are in the explore stage
                !get the new guesses for the mua, mus, n, and g
                
                !choose random mus, mua, hgg, n
                call randomOptProp(mus, mua, hgg, n, reducedmus, mueff, &
                                Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, domain)

                !update optimizationData
                optimizationData(i,1) = mus
                optimizationData(i,2) = mua
                optimizationData(i,3) = hgg
                optimizationData(i,4) = n

                !update the optical properties
                trialOptProp = mono(optimizationData(i,1), optimizationData(i,2), optimizationData(i,3), optimizationData(i,4))
                temp = array(SDF_array_index)%updateOptProp(trialOptProp)

            else
                do while(.true.)
                    !get the new guesses for the mua, mus, n, and g
                    !choose random mus, mua, hgg, n
                    call randomOptProp(mus, mua, hgg, n, reducedmus, mueff, &
                                    Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, &
                                    reducedmusGuessing, mueffGuessing, domain)


                    !find the minimum of left_Min
                    leftMin = optimizationData(1,5) + k * sqrt((mus - optimizationData(1,1))**2 & 
                                                            + (mua - optimizationData(1,2))**2 & 
                                                            + (hgg - optimizationData(1,3))**2 & 
                                                            + (n - optimizationData(1,4))**2)
                    do j = 2,(i-1)
                        tempMin = optimizationData(j,5) + k * sqrt((mus - optimizationData(j,1))**2 &  
                                                                + (mua - optimizationData(j,2))**2 & 
                                                                + (hgg - optimizationData(j,3))**2 & 
                                                                + (n - optimizationData(j,4))**2)
                        if (tempMin < leftMin) then
                            leftMin = tempMin
                        end if
                    end do

                    !LIPO condition
                    if (leftMin >= rightMax) then

                        !update optimizationData
                        optimizationData(i,1) = mus
                        optimizationData(i,2) = mua
                        optimizationData(i,3) = hgg
                        optimizationData(i,4) = n

                        !update the optical properties
                        trialOptProp=mono(optimizationData(i,1), optimizationData(i,2), optimizationData(i,3),&
                                            optimizationData(i,4))
                        temp = array(SDF_array_index)%updateOptProp(trialOptProp)

                        !exit while loop
                        exit

                    end if
                end do
            end if

            !evaluate, by running MCRT
            call random_seed(get=seed) !store the random seed for optical properties
            call run_MCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)
            call random_seed(put=seed) !restart the random seed for optical properties
            error = 0._wp
            call inverse_evaluate(dects, error)
            optimizationData(i,5) = error

            ! reset the arrays storing data
            call reset_detectors(dects) 

            print*, " "
            print*, "Last Guess", i
            print*, "mus", optimizationData(i, 1)
            print*, "mua", optimizationData(i, 2)
            print*, "hgg", optimizationData(i, 3)
            print*, "n", optimizationData(i, 4)
            print*, "mus'", optimizationData(i, 1)*(1-optimizationData(i, 3))
            print*, "mueff", mueff, sqrt(3.0_wp*optimizationData(i,2)*(optimizationData(i,2)+&
                            (optimizationData(i,1)*(1-optimizationData(i,3)))))
            print*, "error", optimizationData(i, 5)
            print*, "k", k
            print*, " "

            !check if this is a new minimum score
            if(optimizationData(i,5) > minError) then
                minError = optimizationData(i,5)
                indexOfMinError = i
                rightMax = optimizationData(i,5)
            end if

            ! find the ratios, and find the new maximum ratio
            probability = 1.0_wp/log(real(i, kind=wp))

            
            do j = 1, i-1
                if (sqrt((optimizationData(i,1) - optimizationData(j,1))**2 & 
                + (optimizationData(i,2) - optimizationData(j,2))**2 & 
                + (optimizationData(i,3) - optimizationData(j,3))**2 & 
                + (optimizationData(i,4) - optimizationData(j,4))**2) == 0.0_wp) then
                    print*, "zero"
                    ratios(ratioCounter) = -99999._wp
                    ratioCounter = ratioCounter + 1
                    cycle
                end if
                ratios(ratioCounter) = abs(optimizationData(i,5) - optimizationData(j,5))/ & 
                                sqrt((optimizationData(i,1) - optimizationData(j,1))**2 & 
                                + (optimizationData(i,2) - optimizationData(j,2))**2 & 
                                + (optimizationData(i,3) - optimizationData(j,3))**2 & 
                                + (optimizationData(i,4) - optimizationData(j,4))**2)
                if (ratios(ratioCounter) > maxRatio) then
                    maxRatio = ratios(ratioCounter)
                    indexOfMaxRatio = ratioCounter
                end if
                ratioCounter = ratioCounter + 1
            end do

            !update the value of k
            it = int(ceiling(log(maxRatio)/log(1+alpha)))
            k = (1+alpha)** it

            print*, " "
            print*, "Best Guess", indexOfMinError
            print*, "mus", optimizationData(indexOfMinError, 1)
            print*, "mua", optimizationData(indexOfMinError, 2)
            print*, "hgg", optimizationData(indexOfMinError, 3)
            print*, "n", optimizationData(indexOfMinError, 4)
            print*, "mus'", optimizationData(indexOfMinError, 1)*(1-optimizationData(indexOfMinError, 3))
            print*, "mueff", sqrt(3.0_wp*optimizationData(indexOfMinError,2)*(optimizationData(indexOfMinError,2)+&
                                        (optimizationData(indexOfMinError,1)*(1-optimizationData(indexOfMinError,3)))))
            print*, "error", optimizationData(indexOfMinError, 5)
            print*, " "
            print*, " "

            !check if we have reached an accuracy value less than the target accuracy, if true then break
            if (abs(minError) > (1.0_wp-accuracy)) then
                print*, "Below min error threshold"

                numGuesses = i
                call write_inverse(optimizationData, outputFile, numGuesses, indexOfMinError)
                return
            end if

        end do        
        
        !we have finished the gradient descent output the optimizationData to a file and write the error
        numGuesses = i - 1
        call write_inverse(optimizationData, outputFile, numGuesses, indexOfMinError)

    end subroutine AdaLIPO_Method

    subroutine Bayesian_Method(input_file, history, packet, dict, distances, image, dects, array, nscatt, start, tev, spectrum, &
        seed, NoVariablesToOptimize, trialOptProp, SDF_array_index, optimizationData, Tmus, Tmua, Thgg, Tn, findmua, findmus, &
        findg, findn, reducedmusGuessing, mueffGuessing, domain, outputFile, layer, accuracy, maxNumSteps, trainingDataSize, &
        fittingDataSize, observationNoise)

        !Shared data
        use iarray
        use constants, only : wp, fileplace, TWOPI

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use opticalProperties, only : opticalProp_t, mono
        use setupMod, only : setup_inverseDirectory
        use writer_mod,    only : write_inverse
        use default_MCRTMod, only : run_MCRT
        use kernels, only : setup, finalise, reset_detectors
        use Interfaces, only : sposv

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, get_value
#ifdef _OPENMP
        use omp_lib
#endif

        character(len=*),              intent(in) :: input_file
        type(history_stack_t),         intent(inout) :: history
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt, start
        type(spectrum_t),              intent(inout) :: spectrum
        type(tevipc),                  intent(inout) :: tev

        integer, allocatable, intent(inout) :: seed(:)
        integer, intent(in) :: SDF_array_index
        integer, intent(in) :: NoVariablesToOptimize
        type(opticalProp_t), intent(inout) :: trialOptProp


        real(kind=wp), intent(inout) :: optimizationData(:,:)
        real(kind=wp), intent(in) :: Tmus, Tmua, Thgg, Tn
        logical, intent(in) :: findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing
        real(kind=wp), intent(in) :: domain(:,:)
        character(len=:), allocatable, intent(inout) :: outputFile
        integer, intent(in) :: maxNumSteps
        integer, intent(in) :: layer
        real(kind=wp), intent(in) :: accuracy

        !unique to Bayesian_Method
        integer, intent(in) :: trainingDataSize
        integer :: fittingDataSize
        real(kind=wp), intent(in) :: observationNoise


        integer :: i, j
        real(kind=wp) :: mus, mua, hgg, n, reducedmus, mueff !Current guess values
        real(kind=wp) :: temp !temp variable for changing the optical properties of the simulation
        real(kind=wp) :: error
        integer :: numGuesses

        !unique to Bayesian_Method
        real(kind=wp) :: tune
        integer :: x, y, count
        real(kind=wp), allocatable :: trainingData(:,:), temptrainingData(:,:)
        real(kind=wp), allocatable :: trainingDataError(:,:)
        real(kind=wp), allocatable :: fittingData(:,:)

        real(kind=sp), allocatable :: kernelObs(:,:)
        real(kind=sp), allocatable :: kernelObstoPred(:,:)
        real(kind=sp), allocatable :: solved(:,:)
        real(kind=sp), allocatable :: tempKernelObs(:,:)
        integer :: INFO
        real(kind=wp), allocatable :: mean(:,:)
        real(kind=sp), allocatable :: kernelPred(:,:)
        real(kind=sp), allocatable :: solvedmatmulkernelObstoPred(:,:)
        real(kind=wp), allocatable :: covariance(:,:)
        real(kind=wp), allocatable :: std(:,:)
        real(kind=wp), allocatable :: expectedImp(:)

        real(kind=wp) :: bestGuess, minError
        integer :: maxExpectedImpIndx, indexOfMinError




        !fill optimizationData with the training data
        count = 0
        do i = 1, trainingDataSize
            call randomOptProp(mus, mua, hgg, n, reducedmus, mueff, Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, &
                                reducedmusGuessing, mueffGuessing, domain)
            
            optimizationData(i,1) = mus
            optimizationData(i,2) = mua
            optimizationData(i,3) = hgg
            optimizationData(i,4) = n
            
            trialOptProp = mono(optimizationData(i,1), optimizationData(i,2), optimizationData(i,3), optimizationData(i,4))
            temp = array(SDF_array_index)%updateOptProp(trialOptProp)

            !evaluate, by running MCRT
            call random_seed(get=seed) !store the random seed for optical properties
            call run_MCRT(input_file, history, packet, dict, & 
                            distances, image, dects, array, nscatt, start, & 
                            tev, spectrum)
            call random_seed(put=seed) !restart the random seed for optical properties
            error = 0._wp
            call inverse_evaluate(dects, error)
            optimizationData(i,5) = error
            ! reset the arrays storing data
            call reset_detectors(dects)

            print*, " "
            print*, "mus: ", mus
            print*, "mua: ", mua
            print*, "hgg: ", hgg
            print*, "n: ", n
            print*, "reducedmus: ", reducedmus
            print*, "mueff: ", mueff
            print*, "error: ", optimizationData(i,5)
        end do
              

        print*, " "
        print*, " "
        print*, "training data over"
        print*, " "
        print*, " "

        

        !store the position of minimum error
        minError = optimizationData((maxloc(optimizationData(1:trainingDataSize,5),dim=1)), 5)
        indexOfMinError = (maxloc(optimizationData(1:trainingDataSize,5),dim=1))
        
        !do BayesianOptimization
        do i = trainingDataSize+1, maxNumSteps
            
            tune = 1.0_wp-((i-trainingDataSize)/(maxNumSteps-trainingDataSize))
            
            if(allocated(trainingData)) deallocate(trainingData)
            allocate(trainingData(i-1, 4))

            if(allocated(trainingDataError)) deallocate(trainingDataError)
            allocate(trainingDataError(i-1, 1))

            ! create training data arrays
            trainingData = 0.0_wp
            trainingData = optimizationData(1:i-1, 1:4)
            trainingDataError = 0.0_wp
            trainingDataError = optimizationData(1:i-1, 5:5)

            
            ! create array of data to be fitted to the model
            if(allocated(fittingData)) deallocate(fittingData)
            allocate(fittingData(fittingDataSize,4))
            do j = 1, fittingDataSize
                call randomOptProp(mus, mua, hgg, n, reducedmus, mueff, Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, &
                                reducedmusGuessing, mueffGuessing, domain)

                fittingData(j,1) = mus
                fittingData(j,2) = mua
                fittingData(j,3) = hgg
                fittingData(j,4) = n
            end do

            ! Kernel of the observations
            call produceKernel(trainingData, trainingData, kernelObs)
            do x = 1, size(kernelObs, dim=1)
                do y = 1, size(kernelObs, dim=2)
                    if (x==y) kernelObs(x,y) = kernelObs(x, y) + real(observationNoise**2)
                end do
            end do 


            ! Kernel of observations to predictions
            call produceKernel(trainingData, fittingData, kernelObstoPred)


            if(allocated(solved)) deallocate(solved)
            allocate(solved(size(kernelObstoPred, dim=1), size(kernelObstoPred, dim=2)))
            solved = kernelObstoPred

            if(allocated(tempKernelObs)) deallocate(tempKernelObs)
            allocate(tempKernelObs(size(kernelObs, dim=1), size(kernelObs, dim=2)))
            tempKernelObs = kernelObs

            

            call sposv("U", size(tempKernelObs, dim = 1), size(solved, dim = 2), &
                    tempKernelObs, size(tempKernelObs, dim=2), solved, size(solved, dim=1), INFO)

            if (INFO /= 0) then
                do x = 1, size(kernelObs, dim=1)
                    do y = 1, size(kernelObs, dim=2)
                        if (kernelObs(x,y) < 0.0_sp) then 
                            print*, " "
                            print*, x, y
                            print*, kernelObs(x,y)
                        end if
                    end do
                end do

                print*, "kernelObstoPred"

                do x = 1, size(kernelObstoPred, dim=1)
                    do y = 1, size(kernelObstoPred, dim=2)
                        if (kernelObstoPred(x,y) < 0.0_sp) then 
                            print*, " "
                            print*, x, y
                            print*, kernelObstoPred(x,y)
                        end if
                    end do
                end do

                print*, "error could not perform solution of linear equations in fitting Bayesian Model"
                print*, "INFO: ", INFO
                print*, "= 0:  successful exit"
                print*, "< 0:  if INFO = -i, the i-th argument had an illegal value"
                print*, "> 0:  if INFO = i, the leading principal minor of order i of &
                        &A is not positive, so the factorization could not be completed, &
                        &and the solution has not been computed."
                exit
            end if

            solved = transpose(solved)
            mean = matmul(solved, trainingDataError)
            
            ! Kernel of predictions to predictions
            call produceKernel(fittingData, fittingData, kernelPred)
            solvedmatmulkernelObstoPred = matmul(solved, kernelObstoPred)

            if(allocated(covariance)) deallocate(covariance)
            allocate(covariance(size(kernelPred, dim =1), size(kernelPred, dim =2)))
            if(allocated(std)) deallocate(std)
            allocate(std(size(kernelPred, dim =1), 1))

            !calculate covariance matrix and standard devaition from sqrt(diag(covariance))

            covariance = 0.0_wp
            std = 0.0_wp
            do x = 1, size(kernelPred, dim =1)
                do y = 1, size(kernelPred, dim =2)
                    covariance(x,y) = kernelPred(x,y) - solvedmatmulkernelObstoPred(x,y)

                    if (x==y) std(x,1) = sqrt(covariance(x,y))
                end do
            end do

            
            !calculate expected improvement
            if(allocated(expectedImp)) deallocate(expectedImp)
            allocate(expectedImp(size(mean, dim=1)))
            bestGuess = minError
            do j = 1, size(mean, dim=1)
                expectedImp(j) = (mean(j,1)-bestGuess-tune)*(0.5*(1+erf((mean(j,1)-bestGuess-tune)/ &
                                                                (sqrt(2.0)*(std(j,1) + 1e-8_wp))))) &
                                + (std(j,1)+1e-8_wp)*(1/sqrt(TWOPI))*exp(-((mean(j,1)-bestGuess-tune)/(std(j,1) + 1e-8_wp))**2/2)
            end do


            !find the maximum of the acquisition function
            maxExpectedImpIndx = maxloc(expectedImp, dim = 1)
            mus = fittingData(maxExpectedImpIndx,1)
            mua = fittingData(maxExpectedImpIndx,2)
            hgg = fittingData(maxExpectedImpIndx,3)
            n = fittingData(maxExpectedImpIndx,4)
            reducedmus = mus * (1-hgg)
            mueff = sqrt(3.0_wp*mua*(mua+reducedmus))

            optimizationData(i,1) = mus
            optimizationData(i,2) = mua
            optimizationData(i,3) = hgg
            optimizationData(i,4) = n

            trialOptProp = mono(optimizationData(i,1), optimizationData(i,2), optimizationData(i,3), optimizationData(i,4))
            temp = array(SDF_array_index)%updateOptProp(trialOptProp)

            !evaluate maximum of acquisition function, by running MCRT
            call random_seed(get=seed) !store the random seed for optical properties
            call run_MCRT(input_file, history, packet, dict, & 
                            distances, image, dects, array, nscatt, start, & 
                            tev, spectrum)
            call random_seed(put=seed) !restart the random seed for optical properties
            error = 0._wp
            call inverse_evaluate(dects, error)
            optimizationData(i,5) = error
            ! reset the arrays storing data
            call reset_detectors(dects)

            print*, "Last Guess", i
            print*, "mus: ", mus
            print*, "mua: ", mua
            print*, "hgg: ", hgg
            print*, "n: ", n
            print*, "reducedmus", reducedmus
            print*, "mueff", mueff, sqrt(3.0_wp*optimizationData(i,2)*(optimizationData(i,2)+&
                                        (optimizationData(i,1)*(1-optimizationData(i,3)))))
            print*, "error: ", optimizationData(i,5)

            indexOfMinError = (maxloc(optimizationData(1:i,5),dim=1))
            minError = optimizationData(indexOfMinError,5)

            print*, " "
            print*, "Best Guess", indexOfMinError
            print*, "mus: ", optimizationData(indexOfMinError,1)
            print*, "mua: ", optimizationData(indexOfMinError,2)
            print*, "hgg: ", optimizationData(indexOfMinError,3)
            print*, "n: ", optimizationData(indexOfMinError,4)
            print*, "reducedmus:", optimizationData(indexOfMinError,1)*(1-optimizationData(indexOfMinError,3))
            print*, "mueff", sqrt(3.0_wp*optimizationData(indexOfMinError,2)*(optimizationData(indexOfMinError,2)+&
                                        (optimizationData(indexOfMinError,1)*(1-optimizationData(indexOfMinError,3)))))
            print*, "error:", optimizationData(indexOfMinError,5)

            !check if we have reached an accuracy value less than the target accuracy, if true then break
            if (abs(minError) > (1.0_wp-accuracy)) then
                print*, "Below min error threshold"

                numGuesses = i
                call write_inverse(optimizationData, outputFile, numGuesses, indexOfMinError)
                return
            end if
        end do

        !we have finished the gradient descent output the optimizationData to a file and write the error
        numGuesses = i - 1
        call write_inverse(optimizationData, outputFile, numGuesses, indexOfMinError)

    end subroutine Bayesian_Method

    subroutine AdaLIPO_with_BayesianTrust_Method(input_file, history, packet, dict, distances, image, dects, array, nscatt, &
        start, tev, spectrum, seed, NoVariablesToOptimize, trialOptProp, SDF_array_index, optimizationData, Tmus, Tmua, Thgg, &
        Tn, findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, domain, outputFile, layer, accuracy, maxNumSteps, &
        alpha, fittingDataSize, observationNoise, LIPOtrainingDataSize, bayesianMinDist, tune, numAdaLIPOtoBay)

        !Shared data
        use iarray
        use constants, only : wp, fileplace, TWOPI

        !subroutines
        use detectors,     only : dect_array
        use historyStack,  only : history_stack_t
        use photonMod,     only : photon
        use piecewiseMod
        use random,        only : ran2, init_rng
        use sdfs,          only : sdf
        use sim_state_mod, only : state
        use opticalProperties, only : opticalProp_t, mono
        use setupMod, only : setup_inverseDirectory
        use writer_mod,    only : write_inverse
        use default_MCRTMod, only : run_MCRT
        use kernels, only : setup, finalise, reset_detectors
        use Interfaces, only : sposv

        !external deps
        use tev_mod, only : tevipc
        use tomlf,   only : toml_table, get_value
#ifdef _OPENMP
        use omp_lib
#endif

        character(len=*),              intent(in) :: input_file
        type(history_stack_t),         intent(inout) :: history
        type(photon),                  intent(inout) :: packet
        type(toml_table),              intent(inout) :: dict
        real(kind=wp),    allocatable, intent(inout) :: distances(:), image(:,:,:)
        type(dect_array), allocatable, intent(inout) :: dects(:)
        type(sdf),        allocatable, intent(inout) :: array(:)
        real(kind=wp),                 intent(inout) :: nscatt, start
        type(spectrum_t),              intent(inout) :: spectrum
        type(tevipc),                  intent(inout) :: tev

        integer, allocatable, intent(inout) :: seed(:)
        integer, intent(in) :: SDF_array_index
        integer, intent(in) :: NoVariablesToOptimize
        type(opticalProp_t), intent(inout) :: trialOptProp


        real(kind=wp), intent(inout) :: optimizationData(:,:)
        real(kind=wp), intent(in) :: Tmus, Tmua, Thgg, Tn
        logical, intent(in) :: findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing
        real(kind=wp), intent(in) :: domain(:,:)
        character(len=:), allocatable, intent(inout) :: outputFile
        integer, intent(in) :: maxNumSteps
        integer, intent(in) :: layer
        real(kind=wp), intent(in) :: accuracy

        !unique to AdaLIPO_Method
        real(kind=wp), intent(in) :: alpha

        !unique to Bayesian_Method
        integer :: fittingDataSize
        real(kind=wp), intent(in) :: observationNoise

        !uniques to combined
        integer, intent(in) :: LIPOtrainingDataSize
        real(kind=wp), intent(in) :: tune
        real(kind=wp), intent(in) :: bayesianMinDist
        integer, intent(in) :: numAdaLIPOtoBay


        integer :: i, j
        real(kind=wp) :: mus, mua, hgg, n, reducedmus, mueff !Current guess values
        real(kind=wp) :: temp !temp variable for changing the optical properties of the simulation
        real(kind=wp) :: error
        integer :: numGuesses

        !used by AdaLIPO
        real(kind=wp) :: probability, k, it
        integer :: indexOfMinError, indexOfMaxRatio, index, ratioCounter
        real(kind=wp) :: minError, maxRatio
        real(kind=wp), allocatable :: ratios(:)
        real(kind=wp) :: ranNum !used to temporarily store a random number
        real(kind=wp) :: leftMin, rightMax, tempMin !sides of the LIPO condition

        !unique to Bayesian_Method
        integer :: x, y, count
        real(kind=wp), allocatable :: trainingData(:,:), temptrainingData(:,:)
        real(kind=wp), allocatable :: trainingDataError(:,:)
        real(kind=wp), allocatable :: fittingData(:,:)

        real(kind=sp), allocatable :: kernelObs(:,:)
        real(kind=sp), allocatable :: kernelObstoPred(:,:)
        real(kind=sp), allocatable :: solved(:,:)
        real(kind=sp), allocatable :: tempKernelObs(:,:)
        integer :: INFO
        real(kind=wp), allocatable :: mean(:,:)
        real(kind=sp), allocatable :: kernelPred(:,:)
        real(kind=sp), allocatable :: solvedmatmulkernelObstoPred(:,:)
        real(kind=wp), allocatable :: covariance(:,:)
        real(kind=wp), allocatable :: std(:,:)
        real(kind=wp), allocatable :: expectedImp(:)

        real(kind=wp) :: bestGuess
        integer :: maxExpectedImpIndx

        !unique to combined
        real(kind=wp) :: bayesianDomain(6,2)
        real(kind=wp) :: normDist
        integer :: sizeTrustTrainingData
        real(kind=wp), allocatable :: upperConf(:)
        integer :: maxUpperConfIndx



        !choose random mus, mua, hgg, n
        call randomOptProp(optimizationData(1,1), optimizationData(1,2), optimizationData(1,3), optimizationData(1,4), & 
                            reducedmus, mueff, Tmus, Tmua, Thgg, Tn, &
                            findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, domain)

        trialOptProp = mono(optimizationData(1,1), optimizationData(1,2), optimizationData(1,3), optimizationData(1,4))
        temp = array(SDF_array_index)%updateOptProp(trialOptProp)

        !evaluate, by running MCRT
        call random_seed(get=seed) !store the random seed for optical properties
        call run_MCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)
        call random_seed(put=seed) !restart the random seed for optical properties
        error = 0._wp
        call inverse_evaluate(dects, error)
        optimizationData(1,5) = error

        ! reset the arrays storing data
        call reset_detectors(dects)  

        !store the position of minimum error
        minError = optimizationData(1,5)
        indexOfMinError = 1

        !store the maximum value
        rightMax = optimizationData(1,5)

        !allocate the ratio array (it will be of the size maxNumSteps triangular number)
        allocate(ratios(int(nint(maxNumSteps * (maxNumSteps+1)/2.0_wp))))
        ratios = 0.0_wp
        ratioCounter = 1
        maxRatio = 0.0_wp
        indexOfMaxRatio = 1

        do i = 2, maxNumSteps

            if (i > LIPOtrainingDataSize .and. mod(i,numAdaLIPOtoBay) == 0) then
                !use bayesian on the trust region around the best guess so far
                
                
                !we have the best, we need to build a list of all the closest 
                if(allocated(temptrainingData)) deallocate(temptrainingData)
                allocate(temptrainingData((i-1),5))

                !constants for the Bayesian Optimization
                bestGuess = optimizationData(indexOfMinError,5)



                !find the sampled points closest to the current best guess to use as training data for a bayesian model
                sizeTrustTrainingData = 0               
                do j = 1, (i-1)
                    normDist = sqrt(((optimizationData(j,1) - optimizationData(indexOfMinError,1))/(domain(1,2)-domain(1,1)))**2 &
                                  + ((optimizationData(j,2) - optimizationData(indexOfMinError,2))/(domain(2,2)-domain(2,1)))**2 &
                                  + ((optimizationData(j,3) - optimizationData(indexOfMinError,3))/(domain(3,2)-domain(3,1)))**2 &
                                  + ((optimizationData(j,4) - optimizationData(indexOfMinError,4))/(domain(4,2)-domain(4,1)))**2) 

                    
                    if (normDist <= bayesianMinDist) then
                        sizeTrustTrainingData = sizeTrustTrainingData + 1

                        temptrainingData(sizeTrustTrainingData,:) = optimizationData(j,:)
                    end if
                end do


                if (sizeTrustTrainingData < 2) then
                    !choose random mus, mua, hgg, n
                    call randomOptProp(mus, mua, hgg, n, reducedmus, mueff, &
                                    Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, domain)                  
                else

                    if(allocated(trainingData)) deallocate(trainingData)
                    if(allocated(trainingDataError)) deallocate(trainingDataError)
                    allocate(trainingData(sizeTrustTrainingData,4))
                    allocate(trainingDataError(sizeTrustTrainingData,1))

                    trainingData = 0.0_sp
                    trainingData = temptrainingData(1:sizeTrustTrainingData,1:4)
                    trainingDataError = 0.0_wp
                    trainingDataError = temptrainingData(1:sizeTrustTrainingData,5:5)

                    
                    !define the domain overwhich to sample the bayesian function
                    bayesianDomain(1,1) = minval(trainingData(:,1)) - (domain(1,2)-domain(1,1))*0.001
                    if (bayesianDomain(1,1) < 0.0_wp) bayesianDomain(1,1) = 0.0_wp !ensure mus always greater than zero
                    bayesianDomain(1,2) = maxval(trainingData(:,1)) + (domain(1,2)-domain(1,1))*0.001
                    bayesianDomain(2,1) = minval(trainingData(:,2)) - (domain(2,2)-domain(2,1))*0.001
                    if (bayesianDomain(2,1) < 0.0_wp) bayesianDomain(2,1) = 0.0_wp !ensure mua always greater than zero
                    bayesianDomain(2,2) = maxval(trainingData(:,2)) + (domain(2,2)-domain(2,1))*0.001
                    bayesianDomain(3,1) = minval(trainingData(:,3)) - (domain(3,2)-domain(3,1))*0.001
                    if (bayesianDomain(3,1) < -1.0_wp) bayesianDomain(3,1) = -1.0_wp !ensure hgg always greater than -1.0
                    bayesianDomain(3,2) = maxval(trainingData(:,3)) + (domain(3,2)-domain(3,1))*0.001
                    if (bayesianDomain(3,2) > 1.0_wp) bayesianDomain(3,2) = 1.0_wp !ensure hgg always less than 1.0
                    bayesianDomain(4,1) = minval(trainingData(:,4)) - (domain(4,2)-domain(4,1))*0.001
                    if (bayesianDomain(4,1) < 1.0_wp) bayesianDomain(4,1) = 1.0_wp !ensure n always greater than 1.0
                    bayesianDomain(4,2) = maxval(trainingData(:,4)) + (domain(4,2)-domain(4,1))*0.001
                    bayesianDomain(5,1) = bayesianDomain(1,1) * (1.0_wp - bayesianDomain(3,2)) ! minMus' = minMus(1-maxHgg)
                    bayesianDomain(5,2) = bayesianDomain(1,2) * (1.0_wp - bayesianDomain(3,1)) ! maxMus' = maxMus(1-minHgg)
                    bayesianDomain(6,1) = sqrt(3.0_wp*bayesianDomain(2,1)*(bayesianDomain(2,1) + bayesianDomain(5,1))) ! minMueff = sqrt(3*minMua(minMua + minMus'))
                    bayesianDomain(6,2) = sqrt(3.0_wp*bayesianDomain(2,2)*(bayesianDomain(2,2) + bayesianDomain(5,2))) ! maxMueff = sqrt(3*maxMua(maxMua + maxMus'))
                  
                    !create array of data to be fitted to the model
                    if (allocated(fittingData)) deallocate(fittingData)
                    allocate(fittingData(fittingDataSize,4))
                    fittingData = 0.0_wp
                    do j = 1, fittingDataSize
                        call randomOptProp(fittingData(j,1), fittingData(j,2), fittingData(j,3), fittingData(j,4), & 
                                reducedmus, mueff, Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, &
                                reducedmusGuessing, mueffGuessing, bayesianDomain)
                    end do 
                    
                    !create training data kernel
                    call produceKernel(trainingData, trainingData, kernelObs)
                    do x = 1, size(kernelObs, dim=1)
                        do y = 1, size(kernelObs, dim=2)
                            if (x==y) kernelObs(x,y) = kernelObs(x, y) + real(observationNoise**2)
                        end do
                    end do


                    !create training to fitting data kernel
                    call produceKernel(trainingData, fittingData, kernelObstoPred)

                    !allocate solved as kernelObstoPred for use in sposv
                    if(allocated(solved)) deallocate(solved)
                    allocate(solved(size(kernelObstoPred, dim=1), size(kernelObstoPred, dim=2)))
                    solved = kernelObstoPred

                    !allocate tempKernelObs as kernelObs for use in sposv
                    if(allocated(tempKernelObs)) deallocate(tempKernelObs)
                    allocate(tempKernelObs(size(kernelObs, dim=1), size(kernelObs, dim=2)))
                    tempKernelObs = kernelObs

                    !solve the linear series of equations
                    call sposv("U", size(tempKernelObs, dim = 1), size(solved, dim = 2), &
                            tempKernelObs, size(tempKernelObs, dim=2), solved, size(solved, dim=1), INFO)

                    solved = transpose(solved)
                    mean = matmul(solved, trainingDataError)

                    ! Kernel of predictions to predictions
                    call produceKernel(fittingData, fittingData, kernelPred)

                    solvedmatmulkernelObstoPred = matmul(solved, kernelObstoPred)

                    !allocate covariance and std vectors
                    if(allocated(covariance)) deallocate(covariance)
                    allocate(covariance(size(kernelPred, dim =1), size(kernelPred, dim =2)))
                    if(allocated(std)) deallocate(std)
                    allocate(std(size(kernelPred, dim =1), 1))

                    !calculate covariance matrix and standard devaition from sqrt(diag(covariance))
                    covariance = 0.0_wp
                    std = 0.0_wp
                    do x = 1, size(kernelPred, dim =1)
                        do y = 1, size(kernelPred, dim =2)
                            covariance(x,y) = kernelPred(x,y) - solvedmatmulkernelObstoPred(x,y)

                            if (x==y) std(x,1) = sqrt(covariance(x,y))
                        end do
                    end do
                    
                    !calculate expected improvement
                    if(allocated(expectedImp)) deallocate(expectedImp)
                    allocate(expectedImp(size(mean, dim=1)))
                    
                    do j = 1, size(mean, dim=1)
                        expectedImp(j) = (mean(j,1)-bestGuess-tune)*(0.5*(1+erf((mean(j,1)-bestGuess-tune)/ &
                                                                        (sqrt(2.0)*(std(j,1) + 1e-8_wp))))) &
                                + (std(j,1)+1e-8_wp)*(1/sqrt(TWOPI))*exp(-((mean(j,1)-bestGuess-tune)/(std(j,1) + 1e-8_wp))**2/2)
                    end do

                    !calculate upper confidence bound
                    if(allocated(upperConf)) deallocate(upperConf)
                    allocate(upperConf(size(mean, dim=1)))

                    do j=1, size(mean, dim=1)
                        upperConf(j) = mean(j,1) + tune*std(j,1)
                    end do

                    !use the maximum expected improvedment to chose the next spot to evaluate
                    maxExpectedImpIndx = maxloc(expectedImp, dim = 1)
                    maxUpperConfIndx = maxloc(upperConf, dim=1)
                    mus = fittingData(maxUpperConfIndx,1)
                    mua = fittingData(maxUpperConfIndx,2)
                    hgg = fittingData(maxUpperConfIndx,3)
                    n   = fittingData(maxUpperConfIndx,4)
                    reducedmus = mus * (1.0_wp-hgg)
                    mueff = sqrt(3.0_wp*mua*(mua+reducedmus))
                end if

                !update optimizationData
                optimizationData(i,1) = mus
                optimizationData(i,2) = mua
                optimizationData(i,3) = hgg
                optimizationData(i,4) = n

                !update the optical properties
                trialOptProp = mono(optimizationData(i,1), optimizationData(i,2), optimizationData(i,3), optimizationData(i,4))
                temp = array(SDF_array_index)%updateOptProp(trialOptProp)

            else
                ranNum = ran2()
                if( ranNum <= probability) then
                    !we are in the explore stage
                    !get the new guesses for the mua, mus, n, and g
                    
                    !choose random mus, mua, hgg, n
                    call randomOptProp(mus, mua, hgg, n, reducedmus, mueff, &
                                    Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, domain)

                    !update optimizationData
                    optimizationData(i,1) = mus
                    optimizationData(i,2) = mua
                    optimizationData(i,3) = hgg
                    optimizationData(i,4) = n

                    !update the optical properties
                    trialOptProp=mono(optimizationData(i,1),optimizationData(i,2),optimizationData(i,3),optimizationData(i,4))
                    temp = array(SDF_array_index)%updateOptProp(trialOptProp)

                else
                    do while(.true.)
                        !get the new guesses for the mua, mus, n, and g
                        !choose random mus, mua, hgg, n
                        call randomOptProp(mus, mua, hgg, n, reducedmus, mueff, &
                                        Tmus, Tmua, Thgg, Tn, findmua, findmus, findg, findn, &
                                        reducedmusGuessing, mueffGuessing, domain)


                        !find the minimum of left_Min
                        leftMin = optimizationData(1,5) + k * sqrt((mus - optimizationData(1,1))**2 & 
                                                                + (mua - optimizationData(1,2))**2 & 
                                                                + (hgg - optimizationData(1,3))**2 & 
                                                                + (n - optimizationData(1,4))**2)
                        do j = 2,(i-1)
                            tempMin = optimizationData(j,5) + k * sqrt((mus - optimizationData(j,1))**2 &  
                                                                    + (mua - optimizationData(j,2))**2 & 
                                                                    + (hgg - optimizationData(j,3))**2 & 
                                                                    + (n - optimizationData(j,4))**2)
                            if (tempMin < leftMin) then
                                leftMin = tempMin
                            end if
                        end do

                        !LIPO condition
                        if (leftMin >= rightMax) then

                            !update optimizationData
                            optimizationData(i,1) = mus
                            optimizationData(i,2) = mua
                            optimizationData(i,3) = hgg
                            optimizationData(i,4) = n

                            !update the optical properties
                            trialOptProp=mono(optimizationData(i,1),optimizationData(i,2),optimizationData(i,3),&
                                                optimizationData(i,4))
                            temp = array(SDF_array_index)%updateOptProp(trialOptProp)

                            !exit while loop
                            exit

                        end if
                    end do
                end if
            end if

            !evaluate, by running MCRT
            call random_seed(get=seed) !store the random seed for optical properties
            call run_MCRT(input_file, history, packet, dict, & 
                        distances, image, dects, array, nscatt, start, & 
                        tev, spectrum)
            call random_seed(put=seed) !restart the random seed for optical properties
            error = 0._wp
            call inverse_evaluate(dects, error)
            optimizationData(i,5) = error

            ! reset the arrays storing data
            call reset_detectors(dects) 

            print*, " "
            print*, "Last Guess", i
            print*, "mus", mus, optimizationData(i, 1)
            print*, "mua", mua, optimizationData(i, 2)
            print*, "hgg", hgg, optimizationData(i, 3)
            print*, "n", n, optimizationData(i, 4)
            print*, "mus'", reducedmus, optimizationData(i, 1)*(1-optimizationData(i, 3))
            print*, "mueff", mueff, sqrt(3.0_wp*optimizationData(i,2)*(optimizationData(i,2)+&
                                        (optimizationData(i,1)*(1-optimizationData(i,3)))))
            print*, "error", optimizationData(i, 5)
            print*, "k", k
            print*, " "

            !check if this is a new minimum score
            if(optimizationData(i,5) > minError) then
                minError = optimizationData(i,5)
                indexOfMinError = i
                rightMax = optimizationData(i,5)
            end if

            ! find the ratios, and find the new maximum ratio
            probability = 1.0_wp/log(real(i, kind=wp))

            
            do j = 1, i-1
                if (sqrt((optimizationData(i,1) - optimizationData(j,1))**2 & 
                + (optimizationData(i,2) - optimizationData(j,2))**2 & 
                + (optimizationData(i,3) - optimizationData(j,3))**2 & 
                + (optimizationData(i,4) - optimizationData(j,4))**2) == 0.0_wp) then
                    print*, "zero"
                    ratios(ratioCounter) = -99999._wp
                    ratioCounter = ratioCounter + 1
                    cycle
                end if
                ratios(ratioCounter) = abs(optimizationData(i,5) - optimizationData(j,5))/ & 
                                sqrt((optimizationData(i,1) - optimizationData(j,1))**2 & 
                                + (optimizationData(i,2) - optimizationData(j,2))**2 & 
                                + (optimizationData(i,3) - optimizationData(j,3))**2 & 
                                + (optimizationData(i,4) - optimizationData(j,4))**2)
                if (ratios(ratioCounter) > maxRatio) then
                    maxRatio = ratios(ratioCounter)
                    indexOfMaxRatio = ratioCounter
                end if
                ratioCounter = ratioCounter + 1
            end do

            !update the value of k
            it = int(ceiling(log(maxRatio)/log(1+alpha)))
            k = (1+alpha)** it

            print*, " "
            print*, "Best Guess", indexOfMinError
            print*, "mus", optimizationData(indexOfMinError, 1)
            print*, "mua", optimizationData(indexOfMinError, 2)
            print*, "hgg", optimizationData(indexOfMinError, 3)
            print*, "n", optimizationData(indexOfMinError, 4)
            print*, "mus'", optimizationData(indexOfMinError, 1)*(1-optimizationData(indexOfMinError, 3))
            print*, "mueff", sqrt(3.0_wp*optimizationData(indexOfMinError,2)*(optimizationData(indexOfMinError,2)+&
                                        (optimizationData(indexOfMinError,1)*(1-optimizationData(indexOfMinError,3)))))
            print*, "error", optimizationData(indexOfMinError, 5)
            print*, " "
            print*, " "

            !check if we have reached an accuracy value less than the target accuracy, if true then break
            if (abs(minError) > (1.0_wp-accuracy)) then
                print*, "Below min error threshold"

                numGuesses = i
                call write_inverse(optimizationData, outputFile, numGuesses, indexOfMinError)
                return
            end if

        end do        
        
        !we have finished the gradient descent output the optimizationData to a file and write the error
        numGuesses = i - 1
        call write_inverse(optimizationData, outputFile, numGuesses, indexOfMinError)

    end subroutine AdaLIPO_with_BayesianTrust_Method



    !produce a covariance matrix from two arrays using the guassian kernel
    subroutine produceKernel(pointsArray1, pointsArray2, kernel)

        use constants, only : wp, sp

        real(kind=wp), intent(in) :: pointsArray1(:,:), pointsArray2(:,:)
        real(kind=sp), allocatable, intent(inout) :: kernel(:,:)

        integer :: i, j
        real(kind=wp) :: temp

        if(allocated(kernel)) deallocate(kernel)
        allocate(kernel(size(pointsArray1, dim = 1), size(pointsArray2, dim = 1)))

        do i = 1, size(pointsArray1, dim = 1)
            do j = 1, size(pointsArray2, dim = 1)
                !take exp(-0.5 * square euclidean distance)
                temp = -0.5 *  ((real(pointsArray1(i,1)) - real(pointsArray2(j,1)))**2 + &
                                (real(pointsArray1(i,2)) - real(pointsArray2(j,2)))**2 + &
                                (real(pointsArray1(i,3)) - real(pointsArray2(j,3)))**2 + &
                                (real(pointsArray1(i,4)) - real(pointsArray2(j,4)))**2)

                ! if kernel(i,j) < 1e-10 round to zero
                if (temp < -23.02585) then
                    kernel(i,j) = 0.0
                else 
                    kernel(i,j) = exp(temp)
                end if

            end do
        end do


    end subroutine produceKernel

    !return a random set of optical properties
    subroutine randomOptProp(mus, mua, hgg, n, reducedmus, mueff, Tmus, Tmua, Thgg, Tn, &
                            findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing, domain)

        use random, only : ran2
        use constants, only : wp

        real(kind=wp), intent(out) :: mus, mua, hgg, n, reducedmus, mueff
        real(kind=wp), intent(in) :: Tmus, Tmua, Thgg, Tn
        logical, intent(in) :: findmua, findmus, findg, findn, reducedmusGuessing, mueffGuessing
        real(kind=wp) :: domain(6,2)

        !choose mus, mua, hgg, and n from a random range
        if (findmus) then
            mus = ran2() * (domain(1,2)-domain(1,1)) + domain(1,1)
        else 
            mus = Tmus
        end if 
        if (findg) then
            hgg = ran2() * (domain(3,2)-domain(3,1)) + domain(3,1)
        else 
            hgg = Thgg
        end if 
        if (findmua) then
            mua = ran2() * (domain(2,2)-domain(2,1)) + domain(2,1)
        else
            mua = Tmua
        end if 
        if (findn) then
            n = ran2() * (domain(4,2)-domain(4,1)) + domain(4,1)
        else 
            n = Tn
        end if

        !check if using reduced mus Guessing mode or mueff Guessing mode or both

        if(reducedmusGuessing .and. mueffGuessing) then
            do while (.true.)
                reducedmus = mus * (1.0_wp - hgg)
                mueff = sqrt(3*mua*(mua + reducedmus))
                if (reducedmus >= domain(5,1) .and. reducedmus <= domain(5,2) .and. &
                    mueff >= domain(6,1) .and. mueff <= domain(6,2)) then
                    !mua and mus' are withing bounds, exit while loop
                    exit
                end if

                !otherwise choose new random values for mus, mua, and hgg
                if (findmus) then
                    mus = ran2() * (domain(1,2)-domain(1,1)) + domain(1,1)
                else 
                    mus = Tmus
                end if 
                if (findg) then
                    hgg = ran2() * (domain(3,2)-domain(3,1)) + domain(3,1)
                else 
                    hgg = Thgg
                end if 
                if (findmua) then
                    mua = ran2() * (domain(2,2)-domain(2,1)) + domain(2,1)
                else
                    mua = Tmua
                end if 
            end do
        else if(reducedmusGuessing) then
            do while (.true.)
                reducedmus = mus * (1.0_wp - hgg)
                mueff = sqrt(3*mua*(mua + reducedmus))
                if (reducedmus >= domain(5,1) .and. reducedmus <= domain(5,2)) then
                    !g and mus are within bounds, exit while loop
                    exit
                end if

                !otherwise choose new mus and hgg
                if (findmus) then
                    mus = ran2() * (domain(1,2)-domain(1,1)) + domain(1,1)
                else 
                    mus = Tmus
                end if 
                if (findg) then
                    hgg = ran2() * (domain(3,2)-domain(3,1)) + domain(3,1)
                else 
                    hgg = Thgg
                end if

            end do
        else if(mueffGuessing) then
            do while (.true.)

                reducedmus = mus * (1.0_wp - hgg)
                mueff = sqrt(3*mua*(mua + reducedmus))
                if (mueff >= domain(6,1) .and. mueff <= domain(6,2)) then
                    !mua and mus' are withing bounds, exit while loop
                    exit
                end if

                !choose new mus, mua, and hgg
                if (findmus) then
                    mus = ran2() * (domain(1,2)-domain(1,1)) + domain(1,1)
                else 
                    mus = Tmus
                end if 
                if (findg) then
                    hgg = ran2() * (domain(3,2)-domain(3,1)) + domain(3,1)
                else 
                    hgg = Thgg
                end if 
                if (findmua) then
                    mua = ran2() * (domain(2,2)-domain(2,1)) + domain(2,1)
                else
                    mua = Tmua
                end if 

            end do 
        else 
            reducedmus = mus * (1.0_wp - hgg)
            mueff = sqrt(3*mua*(mua + reducedmus))
        end if

    end subroutine randomOptProp

    subroutine inverse_evaluate(dects, error)
        !Calculate the error between the detector target values and given detector actual values

        use constants, only : wp
        use detectors,     only : dect_array
        use sim_state_mod, only : state

        type(dect_array), allocatable, intent(inout) :: dects(:)
        real(kind=wp), intent(inout) :: error

        integer :: counter, loopCounter
        real(kind = wp) :: total, targetVal

        error = 0._wp
        counter = 0

        do loopCounter = 1, size(dects)
            targetVal = dects(loopCounter)%p%targetValue
            if (targetVal /= -1) then

                !get the total value of the detector
                total = 0._wp
                call dects(loopCounter)%p%total_dect(total)
                total = total / real(state%nphotons, kind=wp)

                !error is defined as average absolute difference between all detectors
                error = error + abs((total-targetVal)/(targetVal+1e-16_wp))             !relative difference
                !error = error + (log(total+1e-16_wp) - log(targetVal+1e-16_wp))**2     !log difference
                !error = error + abs(total-targetVal)                                   !absolute difference

                counter = counter + 1
            end if
        end do

        !! current tested method of error
        error = 1.0_wp-error/counter
        
        !limit error between 1.0_wp and 0.0_wp
        if (error < 0.0_wp) then
            error = 0.0_wp
        end if
    end subroutine inverse_evaluate

end module inverseMCRTMod