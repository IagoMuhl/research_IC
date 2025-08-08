program estado_fundamental
    implicit none
    real*8:: J2, J3, H, E_5 , E_4, E_6, E_PO

     J2 = 0.806122449; J3 = -0.111564626;  H = 0.d0
! 0.806122449
     open(15)

    do while (H<=4)

        E_5 = - J2/2 + (3)*J3/4 - 1.d0/4
        
        E_4 = 3*J3 - H/2

        E_6 = -J2/2 + J3/12 - H/3 + 1.d0/12

        E_PO = 3*J2/2 + 3*J3/4 - H -3.d0/4

        WRITE(15,*) H, E_5, E_4, E_6, E_PO

        H = H + 0.001

    enddo

    close(15)
    
end program estado_fundamental