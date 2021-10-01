PROGRAM neuralModel
    integer, parameter :: dp = selected_real_kind(15, 500)
    integer :: eventSize = 0, systemSize, placeIndex = 0, maxK,placeIndex2=0
    real(dp) :: driving, propogation, decay, refractory, neuronPlaceholder, probFlip
    integer(kind=8), allocatable :: theBrain(:), avalanches(:,:), theActives(:),outputBrain_Array(:,:)
    integer(kind=8), allocatable :: frequencies_Sizes(:,:),frequencies_Durations(:,:)
    integer(kind=8), allocatable :: networkConnected(:,:),degreeK(:),currentNeuron(:),currentNeuronIA(:)
    integer(kind=8) :: timer,currentTime,kSum
    LOGICAL :: coin
    character(len=1024) :: filename1,filename2,filename3,filename4,filename5,filename6, h1,h2,h3    ! f1 - Raw Output Data; f2 - Histrogram_List of Sizes; f3 - Histrogram_List of Durations; f4 - Unfinished Avalanchews

    integer, parameter :: LONG=SELECTED_INT_KIND(18)
    integer(long) :: seed,seedRecord
    real rnd48, r
    integer(kind=8) :: i,k=1,t,tt,p,q,l,j
    integer(kind=8) :: sizeCounter, a1,a2,a3,a4,a5,activityTotal, clusterTotal

    read *,timer, systemSize, driving, propogation, decay, refractory, seed, probFlip               ! Read in system parameters
    seedRecord = seed

    write (filename1, fmt='(a,I0,a)') "brain",seedRecord,"_rawData.txt"                             ! Create and open files to store rawDZata and timeData
    write (filename4, fmt='(a,I0,a)') "brain",seedRecord,"_timeData.txt"                            ! " 
    write (filename5, fmt='(a,I0,a)') "brain",seedRecord,"_network.txt"                             ! " 
    write (filename6, fmt='(a,I0,a)') "brain",seedRecord,"_netCon.txt"                              ! " 
    open(unit=5,file=filename1)                                                                     ! " 
    open(unit=4,file=filename4)                                                                     ! " 
    open(unit=24,file=filename5)                                                                    ! " 
    open(unit=25,file=filename6)                                                                    ! " 

    write(5,*) "Number of loads: ", timer, "; seed: ", seedRecord                                   ! Write parameters to lines 1-3 of textfile
    write(5,*) "SystemSize, driving, propogation, decay, refractory, probability"                   ! " 
    write(5,*) systemSize, driving, propogation, decay, refractory, probFlip                        ! " 

    allocate(theBrain(systemSize), theActives(systemSize),avalanches(5,systemSize))                 ! Initialize necessary arrays for siteIndex and avalancheNumberIndex
                                                                                                    ! For avalanches, each column is {avalancheIndex, avalancheSize, activated-deactivated/time, end-time, start-time}

    write(4,*) "Number of loads: ", timer, "; seed: ", seedRecord
    write(4,*) "SystemSize, driving, propogation, decay, refractory"
    write(4,*) systemSize, driving, propogation, decay, refractory

! ==================================================================================================
                    
! Initialize the Brain                                                                              !!!!!!!!!!!! Initialize the Brain !!!!!!!!!!!!
    do i=1,systemSize
        theBrain(i) = 0
    end do

! Initialize the connections network
    allocate(networkConnected(systemSize,systemSize),degreeK(systemSize))                            ! Prepare N x N matrix for network connections. 

    do i=1,systemSize                                                                               ! Run over columns, i
        do j=i+1,systemSize                                                                         ! Run over rows, j > i
            coin = rnd48(seed) <= probFlip                                                          ! For connection ij, flip a coin
            if (coin.eqv. .True.) then                                                              ! For coin=heads, ij are connected. Else disconnected. 
                networkConnected(i,j) = 1
            end if
        end do
    end do

    do i=1,systemSize                                                                               ! Set lower left triangle to equal upper right triangle
        do j=i+1,systemSize
            networkConnected(j,i) = networkConnected(i,j)
        end do
    end do

    do i=1,systemSize                                                                               ! Print the connection degree distribution
        kSum = 0
        do j=1,systemSize
            kSum = kSum + networkConnected(i,j)
            if (networkConnected(i,j)==1) then 
                write(25,'(I5)',advance="no")  j
                ! 10 FORMAT(' ',I5)
            end if 
        end do
        degreeK(i) = kSum
        write(24,*) kSum 
        write(25,*)  
    end do
    ! write(25,*) degreeK
    close(25)
    ! write(24,*) networkConnected 
    write(6,*) "Finished creating the network"
    allocate(currentNeuron(systemSize))                                                             ! Allocate space for all current neuron's connections 
    allocate(currentNeuronIA(systemSize))                                                           ! Allocate space for all current neuron's status (Inactive/Active) 
    do i=1,systemSize                                                                               ! Run over columns, i
        currentNeuronIA(i) = -1
    end do

! Initial Activation                                                                                !!!!!!!!!!!! Initial Activation !!!!!!!!!!!!
    placeIndex = nint(rnd48(seed)*systemSize)                                                       ! Choose a random neuron site
    theBrain(placeIndex) = 1                                                                        ! Activate the site
    maxK = 1                                                                                        ! Initialize the highest K
    theActives(placeIndex) = maxK                                                                   ! Set avalancheIndex to newly active site
    avalanches(1,maxK) = 1                                                                          ! Initialize the {size,duration} of avalanche
    avalanches(2,maxK) = 1
    avalanches(3,maxK) = 1
    placeIndex = 0                                                                                  ! Reset the placeIndex
    currentTime=0
    sizeCounter=0
    activityTotal=1
    clusterTotal=1

! =========================================================================================================================================================================
! =========================================================================================================================================================================

! Run simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                             !!!!!!!!!!!! Run simulation !!!!!!!!!!!!
    do while (sizeCounter<timer)                                                                    ! Run for at up to 100,000 events, or until end of timer
        currentTime = currentTime+1

!!!!!!!!!!!!!!!!!!!!!!!! move states from refractory !!!!!!!!!!!!!!!!!!!!!!!!
        do j=1,systemSize                                                                           ! loop over every  neuron availible
            if (theBrain(j)<0) then                                                                 ! check if neuron j is refractory
                theBrain(j) = theBrain(j) + 1                                                       ! Change state: refractory -> inactive
            end if
        end do

!!!!!!!!!!!!!!!!!!!!!!!! initiate a new avalanche +k !!!!!!!!!!!!!!!!!!!!!!!!
        if (driving >= rnd48(seed)) then                                                            ! Generate random number [0,1), if driving is greater, start an activation.
            placeIndex = nint(rnd48(seed)*systemSize)
            if (theBrain(placeIndex) == 0) then
                write(6,*) currentTime, activityTotal, clusterTotal, 'attempting to activate'
                maxK = maxK + 1                                                                     ! Initialize a new avalanche index
                k = maxK                                                                            ! Set new avalanche index to maxK (highest existing thus far)
                theBrain(placeIndex) = 1                                                            ! Activate site at 'placeIndex'
                theActives(placeIndex) = maxK                                                       ! Set new active site to avalancheIndex on 'theActives'
                activityTotal = activityTotal + 1
                clusterTotal = clusterTotal + 1
                do l=1,systemSize                                                                   ! Find an empty spot in the avalanches array to save
                    if (avalanches(1,l) == 0) then                                                  ! Is this slot (avalancheIndex tracker) empty?
                        avalanches(1,l) = maxK                                                      ! Initialize the index of the new avalanche
                        avalanches(2,l) = 1                                                         ! Initialize the size of the new avalanche
                        avalanches(3,l) = 1                                                         ! Initialize the active-deactive counter of the new avalanche
                        avalanches(4,l) = currentTime                                               ! Initialize the start time of the new avalanche
                        avalanches(5,l) = currentTime                                               ! Initialize the end time of the new avalanche
                        exit                                                                        ! Exit from the do loop
                    end if
                end do 
            end if 
            placeIndex = 0                                                                          ! Reset place index
            k = 0                                                                                   ! Reset the avalancheIndex
            write(4,*) currentTime, activityTotal, clusterTotal, 'activated'
            write(6,*) currentTime, activityTotal, clusterTotal, 'activated'
        end if

!!!!!!!!!!!!!!!!!!!!!!!! spread an existing avalanche +1 !!!!!!!!!!!!!!!!!!!!!!!!
        if (propogation >= rnd48(seed)) then                                                        ! Generate random number [0,1), if propogation is greater, start an activation.
            placeIndex = 0
            placeIndex = nint(rnd48(seed)*systemSize)                                               ! Generate random site index
            if (theBrain(placeIndex) == 1) then                                                     ! If the randomly chosen site is active, proceed with spreading activity 
                open(unit=25,file=filename6)                                                        ! Reopen the file 25 containing index of connections 
                do j=1,placeIndex-1
                    read(25,*)  k                                                                   ! Read through all lines up to placeIndex
                    write(6,*) "First lines coming up to: ", placeIndex, k
                end do
                write(6,*) placeIndex, degreeK(placeIndex)

                do j=1,degreeK(placeIndex)                                                          ! 
                    read(25,'(I5)',advance='no') currentNeuron(j)                                   ! For each neuronIndex in currentNeuron, read its status (Inactive/Active) to currentNeuronIA  
                    write(6,*) "Connected from-to: ", placeIndex, currentNeuron(j)
                end do
                
                write(6,*) currentTime, activityTotal, clusterTotal, 'attempting to propogate '
                do j=degreeK(placeIndex)+1,systemSize                                               ! 
                    currentNeuron(j) = -1                                                           ! For each neuronIndex in currentNeuron, read its status (Inactive/Active) to currentNeuronIA  
                end do  
                do j=1,degreeK(placeIndex)                                                          ! 
                    currentNeuronIA(j) = theBrain(currentNeuron(j))                                 ! For each neuronIndex in currentNeuron, read its status (Inactive/Active) to currentNeuronIA  
                end do  
                if (any(currentNeuronIA==0)) then
                    do while (placeIndex2 == 0)                                                     ! find the random inactive site to spread to
                        placeIndex2 = nint(rnd48(seed)*degreeK(placeIndex))                         ! generate random site index (from 1 to k, where k-degree of nueron begin spread)
                        if (theBrain(placeIndex2) /= 0) then                                        ! If the randomly chosen site is not inactive, try again
                            placeIndex2 = 0                                                         ! reset the trial index
                        end if 
                    end do 
                    theBrain(placeIndex2) = 1                                                       ! Activate the site
                    theActives(placeIndex2) = theActives(currentNeuron(placeIndex2))                ! Assign the site its avalancheIndex
                    activityTotal = activityTotal + 1
                    do l=1,systemSize                                                               ! Find an empty spot in the avalanches array to save
                        if (avalanches(1,l) == theActives(currentNeuron(placeIndex2))) then         ! Is this slot's avalanche index \equiv to the one being spread?
                            avalanches(2,l) = avalanches(2,l) + 1                                   ! Increase the size of the new avalanche by 1
                            avalanches(3,l) = avalanches(3,l) + 1                                   ! Increase the active-deactive counter of the new avalanche by 1
                            exit                                                                    ! Exit from the do loop
                        end if
                    end do
                    placeIndex = 0                                                                  ! Reset place index
                    placeIndex2 = 0                                                                 ! Reset place index
                    k = 0                                                                           ! Reset the avalancheIndex
                    q = 0                                                                           ! Reset the avalancheExistenceChecker
                    write(4,*) currentTime, activityTotal, clusterTotal, 'propagate'
                    write(6,*) currentTime, activityTotal, clusterTotal, 'propagate'
                end if 
                deallocate(currentNeuron)                                                           ! Deallocate currentNeuron for availibility next propagation
                do j=1,degreeK(placeIndex)                                                          ! 
                    currentNeuronIA(j) = -1                                                         ! Deallocate currentNeuronIA for availibility next propagation
                end do
                close(25)                                                                           ! Close connections (to start reading from top next time there's a propagation)                                                          
            end if
        end if

!!!!!!!!!!!!!!!!!!!!!!!! deactive an active neuron -1 !!!!!!!!!!!!!!!!!!!!!!!!
        if (decay >= rnd48(seed)) then                                                              ! Generate random number [0,1), if decay is greater, start an activation.
            placeIndex = nint(rnd48(seed)*systemSize)                                               ! Generate random site index
            if ( theBrain(placeIndex)==1 ) then                                                     !!! Check if the randomly chosen site is active neuron
                write(6,*) currentTime, activityTotal, clusterTotal, "attempting to decay    "
                do l=1,systemSize                                                                       ! Find an empty spot in the avalanches array to save progress
                    if (avalanches(1,l) == theActives(placeIndex)) then                                 ! Is this slot's avalanche index \equiv to the one being deactivated?
                        avalanches(3,l) = avalanches(3,l) - 1                                           ! Decrease the active-deactive counter of the new avalanche by 1
                        avalanches(5,l) = currentTime                                                   ! Record the end time of the latest deactivation
                        if (avalanches(3,l) == 0) then                                                  ! Check if the active-deactive = 0 (i.e. no more active neurons of this avalanche)
                            write(5,*) avalanches(1,l), avalanches(2,l), avalanches(4,l), currentTime, avalanches(3,l) ! Write to file {avalancheIndex, avalancheSize, and avalancheEndTime, avalancheStartTime, active/deactive counter}
                            avalanches(1,l) = 0
                            avalanches(2,l) = 0
                            avalanches(3,l) = 0
                            avalanches(4,l) = 0
                            avalanches(5,l) = 0
                            sizeCounter = sizeCounter + 1
                            clusterTotal = clusterTotal - 1
                        end if
                        exit                                                                        ! Exit from the do loop
                    end if
                end do
                theBrain(placeIndex) = -1 * refractory                                              ! Deactivate the neuron
                theActives(placeIndex) = -1                                                         ! Remove avalancheIndex from inactive Site
                activityTotal = activityTotal - 1
                placeIndex = 0                                                                      ! Reset place index
                k = 0                                                                               ! Reset avalancheIndex
                write(4,*) currentTime, activityTotal, clusterTotal, "decay    "
                write(6,*) currentTime, activityTotal, clusterTotal, "decay    "
            end if
        end if

!!!!!!!!!!!!!!!!!!!!!!!! record total activity !!!!!!!!!!!!!!!!!!!!!!!!
        if (mod(currentTime,1000)==0) then
            write(4,*) currentTime, activityTotal
        end if
    end do
!!!!!!!!!!!!!!!!!!!!!!!! Finished Simulation !!!!!!!!!!!!!
! =========================================================================================================================================================================
! =========================================================================================================================================================================
!!!!!!!!!!!!!!!!!!!!!!!! Process Data to corresponding files !!!!!!!!!!!!!!!!
    allocate(outputBrain_Array(5,sizeCounter))                                                      ! Initialize an array equivalent to the outputBrain file
    close(5)

    !!! Create Histrogram_Lists of avalanche sizes and durations
    write (filename2, fmt='(a,I0,a)') "brain",seedRecord,"_Sizes.txt"
    write (filename3, fmt='(a,I0,a)') "brain",seedRecord,"_Durations.txt"
    open(unit=11,file=filename2)                                                                    ! FrequenciesSizes_1
    open(unit=12,file=filename3)                                                                    ! FrequenciesDurations_1
    open(unit=5,file=filename1)                                                                     ! outputBrain
    read(5,'(a)') h1                                                                                ! Read through first 3 lines containing simulation information
    read(5,'(a)') h2
    read(5,'(a)') h3
    do l=1,sizeCounter-1                                                                            ! Create outputBrain_Array from file outputBrain
        read(5,*) a1,a2,a3,a4,a5                                                                    ! Read {Index, Size, StartTime, Endtime, Active-deactive = 0?}
            outputBrain_Array(1,l) = a1                                                             ! Index
            outputBrain_Array(2,l) = a2                                                             ! Size
            outputBrain_Array(3,l) = a3                                                             ! Start
            outputBrain_Array(4,l) = a4                                                             ! End
            outputBrain_Array(5,l) = a4-a3+1                                                        ! End-start = duration.
    end do

!!!!!!!!!!!!!!!!!!!!!!!! Create Histrogram_List for Sizes !!!!!!!!!!!!!!!!!!!!!!!!
    write(11,*) "Number of loads: ", timer, "; seed: ", seed, "; sizeCounter: ", sizeCounter
    write(11,*) "SystemSize, driving, propogation, decay, refractory"
    write(11,*) systemSize, driving, propogation, decay, refractory

    allocate(frequencies_Sizes(2,maxval(outputBrain_Array(2,:))))                                   ! Initialize array for freqneucies_Sizes to create Histrogram_List of sizes
    do l=1,maxval(outputBrain_Array(2,:))                                                           ! Loop through all possible avalanche Sizes
        frequencies_Sizes(1,l) = l                                                                  ! Set frequency_Sizes(1,l): x-axis = the size
        frequencies_Sizes(2,l) = 0                                                                  ! Set frequency_Sizes(1,l): y-axis = frequency, starting at 0
    end do

    do l=1,sizeCounter-1                                                                            ! Loop through entire outputBrain_Array
        frequencies_Sizes(2,outputBrain_Array(2,l)) = frequencies_Sizes(2,outputBrain_Array(2,l)) + 1   ! For the size at point 'l' in outputBrain_Array, +1 to that size in the frequencies_Sizes
    end do

    do l=1,maxval(outputBrain_Array(2,:))                                                           ! Loop through all possible avalanche Sizes
        if (frequencies_Sizes(2,l) /= 0) then
            write(11,*) frequencies_Sizes(1,l), frequencies_Sizes(2,l)                              ! Record the output in line l of frequencies_Sizes
        end if
    end do
    ! write(6,*) "Finished sizes..."

!!!!!!!!!!!!!!!!!!!!!!!! Create Histrogram_List of Durations !!!!!!!!!!!!!!!!!!!!!!!!

    write(12,*) "Number of loads: ", timer, "; seed: ", seed, "; sizeCounter: ", sizeCounter
    write(12,*) "SystemSize, driving, propogation, decay, refractory"
    write(12,*) systemSize, driving, propogation, decay, refractory

    allocate(frequencies_Durations(2,1+maxval(outputBrain_Array(5,:))))                             ! Initialize array for freqneucies_Sizes to create Histrogram_List pf sizes
    ! write(6,*) "Allocated durations..."
    do l=1,maxval(outputBrain_Array(5,:))                                                           ! Loop through all possible avalanche Sizes
        frequencies_Durations(1,l) = l                                                              ! Set frequency_Sizes(1,l): x-axis = the size
        frequencies_Durations(2,l) = 0                                                              ! Set frequency_Sizes(1,l): y-axis = frequency, starting at 0
    end do
    ! write(6,*) "durations array ready..."
    ! write(6,*) "Size Counter: ", sizeCounter

    do l=1,sizeCounter-1                                                                            ! Loop through entire outputBrain_Array
        ! write(6,*) l, outputBrain_Array(5,l)
        frequencies_Durations(2,outputBrain_Array(5,l)) = frequencies_Durations(2,outputBrain_Array(5,l)) + 1   ! For the size at point 'l' in outputBrain_Array, +1 to that size in the frequencies_Sizes
    end do
    ! write(6,*) "durations array filled..."

    do l=1,maxval(outputBrain_Array(5,:))                                                           ! Loop through all possible avalanche Sizes
        if (frequencies_Durations(2,l) /= 0) then
            write(12,*) frequencies_Durations(1,l), frequencies_Durations(2,l)                      ! Record the output in line l of frequencies_Sizes
        end if
    end do
    ! write(6,*) "durations array written out..."

    close(11)
    close(12)
    close(5)
    close(4)
    close(24)
    close(25)

end program neuralModel
! ==================================================================================================

!!!!!!!!!!!!!!!!!!!!!!!! Extra Functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION rnd48(seed)

  implicit none
  integer, parameter :: LONG=SELECTED_INT_KIND(18)
  integer, parameter :: REAL8=SELECTED_REAL_KIND(15,300)
  integer(long), parameter :: SEED_M=25214903917_LONG, SEED_A=11_LONG
  real(real8), parameter :: TWONEG48=0.35527136788005009E-14_REAL8
  real rnd48
  integer(LONG) :: seed

  seed=IAND(SEED_M*seed+SEED_A,281474976710655_LONG)
  rnd48=TWONEG48*seed

END FUNCTION rnd48