program matrix
    implicit none
    integer,  parameter:: dim = 2
    real, dimension(dim,dim):: A
    integer:: i,j
    real :: ei = -1

    do i = 1, dim
        do j = 1, dim
        
            A(i,j) = ei
            
        end do
    end do

    !print *, A(1,j), A(2,j)
 
    do i = 1, dim
        write(*,*) (A(i,j),j = 1, dim)
    end do

    


end program matrix