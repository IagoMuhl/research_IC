program tarefa
    implicit none
    integer:: n1,n2,n3,n4,n5
    real, dimension(5):: n

print*, 'Dê o valor para n1:'
read(*,*) n1

print*, 'Dê o valor para n2:'
read(*,*) n2

print*, 'Dê o valor para n3:'
read(*,*) n3

print*, 'Dê o valor para n4:'
read(*,*) n4

print*, 'Dê o valor para n5:'
read(*,*) n5

print*, n1,n2,n3,n4,n5

if (n1>n2) then
    if(n1>n3) then
        if(n1>n4) then
            if (n1>n5) then
                print*, 'n1 é o maior valor'
            endif
        endif
    endif
endif
if (n2>n1) then
    if(n2>n3) then
        if(n2>n4) then
            if (n2>n5) then
                print*, 'n2 é o maior valor'
            endif
        endif
    endif
endif
if (n3>n2) then
    if(n3>n1) then
        if(n3>n4) then
            if (n3>n5) then
                print*, 'n3 é o maior valor'
            endif
        endif
    endif
endif
if (n4>n2) then
    if(n4>n3) then
        if(n4>n1) then
            if (n4>n5) then
                print*, 'n4 é o maior valor'
            endif
        endif
    endif
endif
if (n5>n2) then
    if(n5>n3) then
        if(n5>n4) then
            if (n5>n1) then
                print*, 'n5 é o maior valor'
            endif
        endif
    endif
endif



end program tarefa