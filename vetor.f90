program vector
   implicit none

   real, dimension(4) :: V
   real:: elem
   integer:: i
   V = [1,2,3,4]

   do i = 1, 4
      elem = V(i)

      if ( mod(elem,2.) == 0. ) then
        print *, elem, 'É par'

      end if
   end do
end program vector
