program estado_fundamental
    implicit none
    real*8:: J2, J3, H, E_5 , E_4, E_6, E_PO, J1

    J1 = 1.d0; J2 = 0.806122449; J3 = -0.111564626;  H = 0.d0
! 0.806122449
     open(15)

    do while (H<=10)

        E_4 = (3*J3/4.d0)*2 - H/2
        E_5 = (-J1/4.d0 - J2/2.d0 + 3*J3/4.d0)*2
        E_6 = (J1/12.d0 - J2/2.d0 + J3/12.d0)*2 - H/3.d0
        E_PO = (3*J1/4.d0 + 3*J2/2.d0 + 3*J3/4.d0)*2 - H

        WRITE(15,*) H, E_4, E_5, E_6, E_PO

        H = H + 0.001

    enddo

    close(15)
    
end program estado_fundamental