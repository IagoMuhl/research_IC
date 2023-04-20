program nano
   implicit none
   integer:: row, col, linha, coluna


row = 4
col = 4
linha = 2
coluna = 2

!call loop_do(i)

!call charac(row)

!call matrix_charac(linha,coluna)

!call Matrix_num

call pauli_matrices

end 

subroutine loop_do(i)

   do i=1,5
      print*, i
   end do

   write(*,*)(i,i=1,5)
end subroutine

subroutine charac(row)
   implicit none
   integer, intent(in):: row
   character(len=1), dimension(row):: V
   integer:: i

   V = ['A','B','C','D']

   do i=1,4
    print*, V(i)
   end do

   write(*,*)(V(i),i=1,4)

end subroutine

subroutine matrix_charac(linha,coluna)
   implicit none
   integer, intent(in):: linha, coluna
   character(len=1), dimension(linha,coluna):: W
   integer:: i,j



   W = reshape(['A','B','C','D'],[linha,coluna])



   do j=1,2
      do i=1,2
      print*, W(i,j)
      end do
   enddo

   do j=1,2
   write(*,*)(W(i,j),i=1,linha)
   enddo


end subroutine

subroutine Matrix_num
   implicit none
   integer, dimension(2,2)::    M = reshape( (/1,2,3,4/), (/ 2,2/) )
   integer:: i,j


   do j=1,2
      write(*,*)(M(i,j),i=1,2)
   enddo

end subroutine

subroutine pauli_matrices
   implicit none
   integer, dimension(2,2)::    sigma_x, sigma_y, sigma_z
   integer:: i,j

   sigma_x = reshape( (/0,1,1,0/), (/ 2,2/) )

   sigma_y = reshape( (/0,-1,1,0/), (/ 2,2/) )

   sigma_z = reshape( (/1,0,0,-1/), (/ 2,2/) )

   do j=1,2
      write(*,*)(sigma_x(i,j),i=1,2)
   enddo

   print*, '-----------Sigma x------------'
   print*, ''

   do j=1,2
      write(*,*)(sigma_y(i,j),i=1,2)
   enddo

   print*, '-----------Sigma y------------'
   print*, ''
   
   do j=1,2
      write(*,*)(sigma_z(i,j),i=1,2)
   enddo


   print*, '-----------Sigma z------------'

end subroutine