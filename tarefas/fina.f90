program fina
implicit none
    real:: x, y, final
    integer:: i

    y = 153.9

    do i= 1,200

        print*, 'De o valor de x'
        read(*,*) x



    !print*, 'De o valor de v_max'
    !read(*,*) y

    final = x/y

    print*, 'final',final

    end do


end program