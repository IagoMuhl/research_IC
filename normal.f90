program normal
   implicit none

   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites) :: s
   real, parameter:: J1=1.,J2=0.,J3=0
   real, dimension(maxConfig):: H_intra, H_inter, H
   real, dimension(num_sites):: m_guess

   call base(s)

   call HAM_INTRA(J2,s,H_intra)

      m_guess = [1.,-1.,-1.,1.]

   call HAM_INTER(J2,J3,s,m_guess,H_inter)

   H = H_intra + H_inter


   !call print_matrix(maxConfig,num_sites,s)

   call print_matrixH(maxConfig,1,H)

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

subroutine HAM_INTRA(J2,s,H_intra)
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites), intent(in):: s
   integer :: i
   real, intent(in):: J2
   real, parameter:: J1=1.
   real, dimension(maxConfig), intent(out):: H_intra

   !---------------------HAMILTONIANO INTRA-----------------------------
   do i = 1, maxConfig
      H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,1)*s(i,3) + s(i,4)*s(i,3) + s(i,2)*s(i,4)) &
      & + J2*(s(i,1)*s(i,4) + s(i,3)*s(i,2))

   end do
   !---------------------HAMILTONIANO INTRA-----------------------------

   ! do i = 1, maxConfig
   !  H_inter(i) = J1*(((s(i,1)*s(i,2))-(s(i,1)*s(i,2))/2.) &
   !  & + ((s(i,1)*s(i,3))-(s(i,1)*s(i,3))/2.)) &
   !  & + (3*J2)*((s(i,1)*s(i,4))-(s(i,1)*s(i,4))/2.) &
   !   & + (4*J3)*((s(i,1)**2)-((s(i,1))**2.)/2.)
   ! end do

end subroutine

subroutine HAM_INTER(J2,J3,s,mag,H_inter)
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites), intent(in):: s
   integer :: i
   real, intent(in):: J2,J3
   real, parameter:: J1=1.
   real, dimension(num_sites), intent(in):: mag
   real, dimension(maxConfig), intent(out):: H_inter
 

    do i = 1, maxConfig
      H_inter(i) = J1*(((s(i,1)*mag(2))-(mag(1)*mag(2))/2.) &
    & +                ((s(i,1)*mag(3))-(mag(1)*mag(3))/2.)) &
    & +         (3*J2)*((s(i,1)*mag(4))-(mag(1)*mag(4))/2.) &
      & +       (4*J3)*((s(i,1)*mag(1))-((mag(1))**2)/2.) +&

      &            J1*(((s(i,2)*mag(4))-(mag(2)*mag(4))/2.) &
    & +                ((s(i,2)*mag(1))-(mag(2)*mag(1))/2.)) &
    & +         (3*J2)*((s(i,2)*mag(3))-(mag(2)*mag(3))/2.) &
      & +       (4*J3)*((s(i,2)*mag(2))-((mag(2))**2)/2.) +&

      &            J1*(((s(i,3)*mag(1))-(mag(3)*mag(1))/2.) &
    & +                ((s(i,3)*mag(4))-(mag(3)*mag(4))/2.)) &
    & +         (3*J2)*((s(i,3)*mag(2))-(mag(3)*mag(2))/2.) &
      & +       (4*J3)*((s(i,3)*mag(3))-((mag(3))**2)/2.) +&

      &            J1*(((s(i,4)*mag(2))-(mag(4)*mag(2))/2.) &
    & +                ((s(i,4)*mag(3))-(mag(4)*mag(3))/2.)) &
    & +         (3*J2)*((s(i,4)*mag(1))-(mag(4)*mag(1))/2.) &
      & +       (4*J3)*((s(i,4)*mag(4))-((mag(4))**2)/2.)

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
