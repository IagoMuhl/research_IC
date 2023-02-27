program list5
    implicit none

    integer:: i0,i,step

    print*, 'Informe o valor inicial:'
    read(*,*) i0

    print*, 'Informe o valor final:'
    read(*,*) i

    print*, 'Informe o passo de variação da sequência:'
    read(*,*) step

    do while (i0<=i)

        print*, i0
        i0 = i0 + step

    end do

end program list5