program lista2
implicit none

integer:: z
real:: y,x,c
double precision, parameter:: pi = 3.14159265359, degpi = pi/180.d0


!-----------------QUESTÃO 1------------------------
 print*, '------QUESTÃO 1------'
! 
 write(*,*) 'Dê o valor da base'
 read(*,*) x

 write(*,*) 'Dê o valor do expoente'
 read(*,*) y
! 
 z = nint(x)
! 
! 
 print *, x**y ,'Valor real'
 print *, nint(z**y) ,'Valor inteiro'

!-----------------------------------------------------------

 write(*,*) 'Dê o valor da base'
 read(*,*) x
! 
 write(*,*) 'Dê a ordem da raíz'
 read(*,*) y



 if (x<0.0) then
! 
 x = (-1)*x
! 
 z = nint(x)
! 
 print *, (-1)*x**(1/y)
 !print *, (-1)*z**(1/y)
! 
 else
! 
 z = nint(x)
! 
 print *, x**(1/y)
 print *, z**(1/y)
! 
 endif

print*, '---FIM DA QUESTÃO 1---'

!------------------------QUESTÃO 2--------------------

print*, '------QUESTÃO 2------'

write(*,*) 'Dê o valor da base'
read(*,*) x

write(*,*) 'Dê o valor do numerador do expoente'
read(*,*) y

write(*,*) 'Dê o valor do denominador do expoente'
read(*,*) c


 if (x<0.0) then
! 
!     x = (-1)*x
! 
!     print *, (-1)*x**(y/z)
! 
     else
! 
     print *, x**(y/c)
! 
 end if

print*, '---FIM DA QUESTÃO 2---'

!-------------------------QUESTÃO 3--------------------

print*, '------QUESTÃO 3------'

write(*,*) 'Dê o valor de em radianos'
read(*,*) x

print *, 'Valor do cosseno de',x,'radianos é:',cos(x)

!------------------------------------------------------

write(*,*) 'Dê o valor de em graus'
read(*,*) x

y = x*(degpi)

print *, 'Valor do cosseno de',x,'graus é:',cos(y)

end program lista2
