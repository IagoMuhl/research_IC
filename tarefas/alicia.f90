program alicia
    implicit none
    real:: Ds, lambda, Y, equis, d
    integer:: i
    
    Ds=475; lambda=6.328e-5

    equis = 2*Ds*lambda

    do i = 1, 30000

    write(*,*) 'Entre com Y1:'
    read(*,*) Y


    d = equis/Y

    print*, d

    end do

end program