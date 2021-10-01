program neuronModel                                                             ! Program for running the neaural Model fully connected. 
    integer, parameter :: dp = selected_real_kind(15, 500)
    real(dp) :: driving, propogation, decay, refractory
    integer(kind=8) :: systemSize, maxUnique,timeStep
    integer(kind=8), allocatable :: positions(:), distances(:,:), avalanches(:,:)
    integer(kind=8) :: timer,currentTime, currentStep, clusterTotal, clusterCurrent, eventTotal
    character(len=1024) :: filename1,filename4,filename6     ! f1 - Raw Output Data; f2 - Histrogram_List of Sizes; f3 - Histrogram_List of Durations; f4 - Unfinished Avalanches

    INTEGER, PARAMETER :: LONG=SELECTED_INT_KIND(18)
    INTEGER(LONG) :: seed,seedRecord
    REAL rnd48, r
    integer(kind=8) :: k=1,l
    integer(kind=8) :: activityTotal,activeAttemptCount

    read *,timer, systemSize, driving, propogation, decay, seed
    
    !!! This reads in several parameters:
        ! timer = the number of updates for the lattice
        ! systemSize = the total number of neurons
        ! driving = the value of epsilon driving rate
        ! propagation = the value of the lambda propagation parameter
        ! decay = the value of the mu decay parameter
        ! seed = random number seed

    seedRecord = seed
    
    write (filename1, fmt='(a,I0,a)') "brain",seedRecord,"_rawData.txt"
    write (filename4, fmt='(a,I0,a)') "brain",seedRecord,"_timeData.txt"
    
    open(unit=5,file=filename1)
    open(unit=4,file=filename4)

    ! Write parameters to lines 1-3 of textfile
    write(5,*) "Number of loads: ", timer, "; seed: ", seedRecord
    write(5,*) "SystemSize, driving, propogation, decay"
    write(5,*) systemSize, driving, propogation, decay

    write(4,*) "Number of loads: ", timer, "; seed: ", seedRecord
    write(4,*) "SystemSize, driving, propogation, decay"
    write(4,*) systemSize, driving, propogation, decay
    
    
    maxUnique = systemSize 
    allocate(distances(maxUnique,2), positions(maxUnique),avalanches(5,maxUnique))
    !!! For avalanches, each column is {avalancheIndex, avalancheSize, activated-deactivated/time, end-time, start-time}
    
    do l=1,maxUnique
        positions(l) = 0
        distances(l,1) = 0
        distances(l,2) = 0
    end do 

! Initial Activation                                                            !!!!!!!!!!!! Initial Activation !!!!!!!!!!!!
    positions(1) = 1                                                            ! Activate the first site (position)
    distances(1,1) = 1                                                          ! Activate the first site (distance)
    distances(1,2) = 1                                                          ! Activate the first site (distance)
    currentTime=1
    activityTotal=1
    clusterCurrent=1 
    clusterTotal=1


!!!!!!!!!!!!!!!!!!!!!!! Run simulation !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Run simulation !!!!!!!!!!!!
    do while (currentTime<timer)                                                 ! Run for at up to 100,000 events, or until end of timer
        currentTime = currentTime+1

        write(4,*) currentTime, activityTotal, clusterCurrent               
        !!! Records the current values: (current plate update, the total number of active sites, the total number of unique clusters)


        activeAttemptCount = 0

        do l=1,maxUnique
        !!!!!!!!!!!!!!!!!!!!!!!! Take step right or left for each walker !!!!!!!!!!!!!!!!!!!!!!!!
            if (positions(l) .gt. 0) then                                       ! If walker site is alive
                do k=1,positions(l)
                    if (activityTotal .eq. systemSize) then                     ! If total distance is full, dont allow step forward
                        if (rnd48(seed)<decay) then                             ! If decay probability suceeeds... 
                            positions(l) = positions(l) - 1                     ! Add/subtract step to current position 
                            activityTotal = activityTotal -1                    ! Add to total activity/distance
                        end if 
                    else                                                        ! If total distance is not full, proceed normally 
                        if (rnd48(seed)<decay) then                             ! If decay probability succeeds...
                            positions(l) = positions(l) - 1                     ! Add/subtract step to current position 
                            activityTotal = activityTotal - 1                   ! Add to total activity/distance 
                        else                                                    ! If ddecay probability fails, the propagate
                            positions(l) = positions(l) + 1                     ! Add/subtract step to current position 
                            activityTotal = activityTotal + 1                   ! Add to total activity/distance 
                            distances(l,1) = distances(l,1) + 1.                ! Add step to total distance if +1
                        end if 
                    end if

                    if (positions(l) .eq. 0) then                               ! If walker has returned to origin 
                        if (mod(eventTotal,10000) .eq. 0) then 
                            
                            write(5,*) l, distances(l,1), distances(l,2), currentTime,eventTotal 
                            !!! Write the ending walkers: (distance, starttime, endtime). 
                            !!! The relevant quantities are:
                                !!! distance = size of avalanche
                                !!! endtime-starttime = duration of avalanche 
                            
                        end if 
                        distances(l,1) = 0
                        distances(l,2) = 0
                        clusterCurrent = clusterCurrent - 1
                        eventTotal = eventTotal + 1
                    end if 
                end do 
            
        !!!!!!!!!!!!!!!!!!!!!!!!    birth new walkers   !!!!!!!!!!!!!!!!!!!!!!!!
            else if (positions(l) .eq. 0) then                                  ! If walker site is dead 
                if (driving >= rnd48(seed)) then                                ! Generate random number [0,1), if driving is greater, start an activation.
                    if (activityTotal .ne. systemSize) then                     ! If total distance is full, dont allow step forward 
                        positions(l) = 1                                        ! Initialize walker position 
                        distances(l,1) = 1                                      ! Initialize walker distance  
                        distances(l,2) = currentTime                            ! Initialize walker distance  
                        activityTotal = activityTotal + 1                       ! Add to total activity 
                        clusterCurrent = clusterCurrent + 1                     ! Add to current active walkers 
                        clusterTotal = clusterTotal + 1                         ! Add to total active walkers 
                    end if 
                end if  
            end if 
        end do 

    end do
    close(5)
    close(4)

end program neuronModel

!!!!!!!!!!!!!!!!!!!!!!!! Extra Functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION rnd48(seed)

  IMPLICIT NONE
  INTEGER, PARAMETER :: LONG=SELECTED_INT_KIND(18)
  INTEGER, PARAMETER :: REAL8=SELECTED_REAL_KIND(15,300)
  INTEGER(LONG), PARAMETER :: SEED_M=25214903917_LONG, SEED_A=11_LONG
  REAL(REAL8), PARAMETER :: TWONEG48=0.35527136788005009E-14_REAL8
  REAL rnd48
  INTEGER(LONG) :: seed

  seed=IAND(SEED_M*seed+SEED_A,281474976710655_LONG)
  rnd48=TWONEG48*seed

END FUNCTION rnd48
