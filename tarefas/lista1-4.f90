program list
    implicit none

    integer::i
    real::r

    print*, 'Dê valor para i (inteiro):'
    read(*,*) i

    print*, 'Dê valor para r (real):'
    read(*,*) r
    

    print*, 'O valor real de i é:',real(i)
    print*, 'O valor inteiro mais próximo de r é:',nint(r)

    

end program list