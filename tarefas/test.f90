program test
    implicit none
    integer:: i

    do i = 1, 5
        print *, i
        if ( i==3 ) exit
    end do
    print *, "saí do loop"
    
end program test