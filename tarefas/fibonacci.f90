program fibo
    implicit none

    integer:: i,k
    integer, allocatable:: a(:)

    print*, 'Informe quantos termos terá a sequência:'
    read(*,*) i
    
    allocate (a(i))

    a(0) = 0
    a(1) = 1

    do k=2,i,1
        a(k) = a(k-1) + a(k-2)
    end do

    print*, 'Sequência completa:',a
    print*, 'Quinto termo:', a(5)

    print*, 'Digite o termo que gostaria de conferir:'
    read(*,*) k
    print*, a(k)

end program fibo