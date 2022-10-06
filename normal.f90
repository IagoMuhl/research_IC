program normal
   implicit none

   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites) :: s
   integer :: i
   real, parameter:: J1=1.,J2=0.5,J3=0
   real, dimension(maxConfig):: H_intra
   real, dimension(maxConfig):: H_inter
   call base(s)

   !---------------------HAMILTONIANO INTRA-----------------------------
   do i = 1, maxConfig
      H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,1)*s(i,3) + s(i,4)*s(i,3) + s(i,2)*s(i,4)) &
      & + J2*(s(i,1)*s(i,4) + s(i,3)*s(i,2))

   end do
   !---------------------HAMILTONIANO INTRA-----------------------------

   do i = 1, maxConfig
      H_inter(i) = J1*(((s(i,1)*s(i,2))-(s(i,1)*s(i,2))/2.) &
                   & + ((s(i,1)*s(i,3))-(s(i,1)*s(i,3))/2.)) &
      & + (3*J2)*((s(i,1)*s(i,4))-(s(i,1)*s(i,4))/2.) &
      & + (4*J3)*((s(i,1)**2)-((s(i,1))**2.)/2.)
   end do




   call print_matrix(maxConfig,num_sites,s)

   call print_matrixH(maxConfig,1,H_intra)

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

subroutine print_matrix(row,column,matrix)
   integer, intent(in):: row,column
   integer, dimension(row,column), intent(in):: matrix
   integer:: i,j
   do i = 1, row
      write(*,*) (matrix(i,j),j = 1, column)
   end do

end subroutine

subroutine print_matrixH(rowH,columnH,matrixH)
   integer, intent(in):: rowH,columnH
   real, dimension(rowH,columnH), intent(in):: matrixH
   integer:: i,j
   do i = 1, rowH
      write(*,*) (matrixH(i,j),j = 1, columnH)
   end do

end subroutine
