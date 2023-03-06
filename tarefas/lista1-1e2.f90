program lista1
   implicit none
   
   integer::TERMOS,i,res
   integer, allocatable:: T(:)

   allocate (T(TERMOS))

   print*, 'Digite quantos termos precisa:'
   read(*,*) TERMOS
   
   do i=1,TERMOS,1
      print*, 'Dê o valor para o termo:',i
      read(*,*) T(i)
      !print*, T(1), T(3), T(5)
   enddo

   call QUESTION_1(TERMOS,T,res)

   call QUESTION_2(res)

   deallocate (T)

end program lista1

subroutine  QUESTION_1(TERMOS,T,res)
   integer, intent(in):: TERMOS
   integer, dimension(TERMOS), intent(in):: T
   integer::g,res

   if(T(1)>=T(2)) then
      res = T(1)
   else
      res = T(2)
   end if

   do g = 1,TERMOS

   if (res>=T(g)) then
      print*, ''
   else
      res = T(g)
   end if

   end do

   print*, 'O maior valor informado foi:',res

end subroutine

subroutine  QUESTION_2(res)
   integer, intent(in)::res
   integer:: i, k

   ult= mod(res,2)
   
   if(ult==1) then
      print*, 'O maior valor informado é ímpar e não divisível por 4'
   else
      print*, 'O maior valor informado é par'
      if(mod(res,4)==0) then
         print*, 'É divisível por 4'
      else
         print*, 'Não é divisível por 4'
      endif
   endif

   k = 0

   do i = 1, res
      if (mod(res,i)==0.d0) then
         k = k+1
      endif
   end do
   
      if(k>2) then
         print*, res,'não é um número primo'
      else
         print*, res,'é um número primo'
      endif
   

end subroutine

