program normal
   implicit none

   integer, parameter:: minConfig=2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites) :: s
   integer :: i,j


   call base(s)

   do i = 1, maxConfig
      write(*,*) (s(i,j),j = 1, num_sites)
   end do


end program normal

subroutine base(s)

   implicit none
   integer, parameter:: minConfig=2
   integer, parameter:: maxConfig= minConfig**4
   integer:: k
   integer:: s1,s2,s3,s4,m1,m2,m3,m4
   integer, dimension(16,4), intent(out):: s(maxConfig,4)


   k=1
   !!!!!!!!!!    CONTRUÇÃO DA BASE     !!!!!!!!!!!!
   do s1 =-1,1,2
      do s2 =-1,1,2
         do s3 =-1,1,2
            do s4 =-1,1,2

               s(k,1)=s1
               s(k,2)=s2
               s(k,3)=s3
               s(k,4)=s4


               m1=s1
               m2=s2
               m3=s3
               m4=s4

               k=k+1

            end do
         end do
      end do
   end do


end subroutine
