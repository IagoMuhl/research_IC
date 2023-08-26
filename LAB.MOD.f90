program lab
    implicit none
    real:: Alfa, Beta, VL, IL, R0, T, RT, R, Rad

    Alfa = 4.82*10.0**(-3)
    Beta = 6.76*10.0**(-7)
    R0 = 0.259


    do
    
    ! T = 0

    ! print*, "Entre com VL, IL:"
    ! read(*,*) VL, IL

    ! RT = VL/IL

    ! write(*,*) 'RT:', RT

    ! R = RT/R0

    ! write(*,*) 'R:', R

    ! T = (1.0/(2.0*Beta))*(sqrt((Alfa)**(2) + (4.0*Beta*( R-1.0))) - Alfa) + 273.15

    ! write(*,*) 'T:', T

    ! RT = log(RT)
    ! T = log(T)
    print*, 'Entre com T'
    read(*,*) T

    T = log10(T)

    write(*,*) '--------LOGs--------'
    write(*,*) 'Radiancia:', T
    ! write(*,*) 'T:', T

    enddo

end program