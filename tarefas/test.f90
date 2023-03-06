program test
   implicit none
   character(len=2):: state

   write(*,*) 'Entre com o estado'
   read(*,*) state

   if ( state=='af' ) then
      print*, 'af'
   else
      print*, 'erro'
   end if
   
end program test