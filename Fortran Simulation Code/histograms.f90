program neuronModel
    integer, parameter :: dp = selected_real_kind(15, 500)
    integer :: eventSize = 0, systemSize, placeIndex = 0, maxK,placeIndex2=0
    real(dp) :: driving, propogation, decay, neuronPlaceholder
    integer(kind=16), allocatable :: theBrain(:), avalanches(:,:), theActives(:),outputBrain_Array(:,:)
    integer(kind=16), allocatable :: frequencies_Sizes(:,:),frequencies_Durations(:,:)
    integer(kind=8) :: timer,currentTime,nonZero
    character(len=1024) :: filename0, filename1,filename2,filename3,filename4, h1,h2,h3    ! f1 - Raw Output Data; f2 - Histrogram_List of Sizes; f3 - Histrogram_List of Durations; f4 - Unfinished Avalanchews

    INTEGER, PARAMETER :: LONG=SELECTED_INT_KIND(18)
    INTEGER(LONG) :: seed,seedRecord,reason
    REAL rnd48, r
    integer(kind=8) :: i,k=1,t,tt,p,q,l
    integer(kind=8) :: sizeCounter=0, a1,a2,a3,a4,a5, b1,b2,b3,b4

    write(6,*) "Enter: Seed"
    read*, seedRecord

    write (filename1, fmt='(a,I0,a)') "brain",seedRecord,"_rawData.txt"
    write (filename2, fmt='(a,I0,a)') "brain",seedRecord,"_Sizes.txt"
    write (filename3, fmt='(a,I0,a)') "brain",seedRecord,"_Durations.txt"

    open(unit=5,file=filename1)

    print *,'file1 opened'
    ! count up sizeCounter
    read(5,'(a)') h1                                                            ! Read through first 3 lines containing simulation information
    read(5,'(a)') h2
    read(5,'(a)') h3
    print *,'starting to count'
    sizeCounter = 0
    do
        read(5,*,IOSTAT=reason)
        if (reason<0) then
            Exit
        else
            sizeCounter  = sizeCounter + 1
            ! write(*,*) sizeCounter
        end if
    end do

    allocate(outputBrain_Array(5,sizeCounter))                                  ! Initialize an array equivalent to the outputBrain file
    close(5)
    write (6,*) sizeCounter

    !!! Create Histrogram_Lists of avalanche sizes and durations
    ! write (filename2, fmt='(a,I0,a)') "brain_",seedRecord,"_frequencySizes.txt"
    ! write (filename3, fmt='(a,I0,a)') "brain_",seedRecord,"_frequencyDurations.txt"
    open(unit=11,file=filename2)                                                ! FrequenciesSizes_1
    open(unit=12,file=filename3)                                                ! FrequenciesDurations_1
    open(unit=5,file=filename1)                                                 ! outputBrain
    read(5,'(a)') h1                                                            ! Read through first 3 lines containing simulation information
    read(5,'(a)') h2
    read(5,'(a)') h3
    do l=1,sizeCounter-1                                                        ! Create outputBrain_Array from file outputBrain
        read(5,*) a1,a2,a3,a4,a5                                                ! Read {Index, Size, StartTime, Endtime, Active-deactive = 0?}
            outputBrain_Array(1,l) = a1                                         ! Index
            outputBrain_Array(2,l) = a2                                         ! Size
            outputBrain_Array(3,l) = a3                                         ! Start
            outputBrain_Array(4,l) = a4                                         ! End
            outputBrain_Array(5,l) = a4-a3+1                                    ! End-start = duration.
            ! write(6,*) "Reading: ", l, a1, a2, a3
    end do

!!! Create Histrogram_List for Sizes
    write(11,*) h1
    write(11,*) h2
    write(11,*) h3
    
    allocate(frequencies_Sizes(2,maxval(outputBrain_Array(2,:))))               ! Initialize array for freqneucies_Sizes to create Histrogram_List of sizes
    do l=1,maxval(outputBrain_Array(2,:))                                       ! Loop through all possible avalanche Sizes
        frequencies_Sizes(1,l) = l                                              ! Set frequency_Sizes(1,l): x-axis = the size
        frequencies_Sizes(2,l) = 0                                              ! Set frequency_Sizes(1,l): y-axis = frequency, starting at 0
    end do

    do l=1,sizeCounter-1                                                        ! Loop through entire outputBrain_Array
        ! write(6,*) l
        frequencies_Sizes(2,outputBrain_Array(2,l)) = frequencies_Sizes(2,outputBrain_Array(2,l)) + 1   ! For the size at point 'l' in outputBrain_Array, +1 to that size in the frequencies_Sizes
    end do

    do l=1,maxval(outputBrain_Array(2,:))                                       ! Loop through all possible avalanche Sizes
        if (frequencies_Sizes(2,l) /= 0) then
            ! write(6,*) l
            write(11,*) frequencies_Sizes(1,l), frequencies_Sizes(2,l)              ! Record the output in line l of frequencies_Sizes
        end if
    end do
    write(6,*) "Finished sizes..."
    ! deallocate(frequencies_Sizes)

!!! Create Histrogram_List of Durations
    write(12,*) h1
    write(12,*) h2
    write(12,*) h3

    allocate(frequencies_Durations(2,1+maxval(outputBrain_Array(5,:))))         ! Initialize array for freqneucies_Sizes to create Histrogram_List pf sizes
    ! write(6,*) "Allocated durations..."
    do l=1,maxval(outputBrain_Array(5,:))                                       ! Loop through all possible avalanche Sizes
        frequencies_Durations(1,l) = l                                              ! Set frequency_Sizes(1,l): x-axis = the size
        frequencies_Durations(2,l) = 0                                              ! Set frequency_Sizes(1,l): y-axis = frequency, starting at 0
    end do
    write(6,*) "durations array ready..."
    ! write(6,*) "Size Counter: ", sizeCounter

    do l=1,sizeCounter-1                                                        ! Loop through entire outputBrain_Array
        ! write(6,*) l, outputBrain_Array(5,l)
        frequencies_Durations(2,outputBrain_Array(5,l)) = frequencies_Durations(2,outputBrain_Array(5,l)) + 1   ! For the size at point 'l' in outputBrain_Array, +1 to that size in the frequencies_Sizes
    end do
    write(6,*) "durations array filled..."

    do l=1,maxval(outputBrain_Array(5,:))                                       ! Loop through all possible avalanche Sizes
        if (frequencies_Durations(2,l) /= 0) then
            ! write(6,*) l
            write(12,*) frequencies_Durations(1,l), frequencies_Durations(2,l)              ! Record the output in line l of frequencies_Sizes
        end if
    end do
    write(6,*) "durations array written out..."
    ! deallocate(frequencies_Durations)

    close(11)
    close(12)
    close(5)


end program neuronModel