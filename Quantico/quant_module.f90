module QUANTICO
   implicit none
   real, parameter:: J1=1.d0
contains


subroutine print_matrix(A,dim,din)
      implicit none
      integer, intent(in):: dim,din
      real*8, dimension(dim,din), intent(in):: A
      integer:: i, j
      do i=1,dim
         write(*,20)(A(i,j),j=1,din) !para numeros inteiros (4(I3))
      20   format (16(1x, F8.4))
      enddo
end subroutine

subroutine tensorial(A,B,dim,E)
      implicit none
      integer, intent(in):: dim
      real*8, dimension(dim,dim), intent(in):: A
      real*8, dimension(dim,dim), intent(in):: B
      real*8, dimension(dim**2,dim**2), intent(out):: E
      integer:: i, j, k, l, alfa, beta

      do i = 1,dim
         do j = 1,dim
            do k = 1,dim
               do l = 1,dim

                  alfa = dim*(i-1) + k
                  beta = dim*(j-1) + l

                  E(alfa,beta) = A(i,j) * B(k,l)

               end do
            end do
         end do
      end do
end subroutine

subroutine partition(W,T,dim,Z)
      implicit none
      real*8, intent(out):: Z
      integer, intent(in):: dim
      real*8, dimension(dim**4), intent(in):: W
      real*8, intent(in):: T
      real*8:: b
      integer :: i
      b = 1.d0/T
      Z = 0.d0

      do i = 1, dim**4
         Z = Z + (dexp(-b*(W(i))))
      end do

end subroutine

subroutine magnetization(W,Z,T,dim,s,m)
      implicit none
      integer, intent(in):: dim
      real*8, dimension(dim**4,dim**4), intent(in):: s
      real*8, dimension(dim**4):: W
      real*8, intent(in):: T, Z
      real*8, intent(out):: m
      integer :: i
      real*8:: b

      b = 1.d0/T
      m = 0

      do i = 1, dim**4
         m = m + s(i,i)*(dexp(-b*(W(i))))
      end do

      m = m/Z

end subroutine

subroutine mag_vector(state,m1,m)
   implicit none
   character(len=3), intent(in):: state
   real*8, dimension(2),intent(out):: m
   real*8, intent(in):: m1

   select case(state)
   case('AF')
      m = [m1,-m1]

   case ('SAF')
      m = [m1,m1]

   case ('PM')
      m = 0
   
   case default
      print *, 'State inválido'

end select



end subroutine

subroutine Ham_inter_state(state,J2,J3,s1,s2,s3,s4,H,m1,m2,Id_4,H_inter)
      implicit none
      character(len=3), intent(in):: state
      real*8, dimension(16,16), intent(in):: s1,s2,s3,s4,Id_4
      real*8, intent(in):: J2, J3, m1, m2 , H
      real*8, dimension(16,16), intent(out):: H_inter

      H_inter = 0.d0

      ! if (H/=0.0) then

      select case (state)
       case ('AF')

         H_inter = (2*J1)*(m2*(s1+s4)+m1*(s2+s3)-(2*m1*m2*Id_4))+(3*J2 + 4*J3)*(m1*(s1+s4-(m1*Id_4))+m2*(s2+s3-(m2*Id_4)))

       case ('PM')
         H_inter = 0
      
      !  case default
      !    print *, 'State inválido'
!
         !-H*(s1+s2+s3+s4) 
!

! !
!        case ('AF')
!          H_inter = (-2*J1 + 3*J2 + 4*J3)*(s1-s2-s3+s4-2*m1*Id_4)*(m1)

      !  case ('SAF')
      !    H_inter = (-3*J2 + 4*J3)*(s1-s2+s3-s4-2*m1*Id_4)*(m1)

      !  case ('SD')
      !    H_inter = (-2*J1 + J2)*(s1+s2-s3-s4-2*m1*Id_4)*(m1)

      !  case ('PM')
      !    H_inter = 0
       case default
         print *, 'State inválido'
      end select

!

end subroutine



subroutine Free_nrg(T,Z,F)
      implicit none
      real*8, intent(in):: Z,T
      real*8, intent(out):: F

      F = -T*dlog(Z)

end subroutine

subroutine diagonalization(A, V, W)
   IMPLICIT NONE
   integer, parameter :: dim = 2, LWMAX = 1000
   real*8, dimension (dim**4,dim**4), intent(in) :: A
   real*8, dimension(dim**4,dim**4), intent(out):: V
   ! !!! Para diagonalizar !!!!!
   character(len=1):: JOBZ , UPLO
   integer:: N , LDA , INFO , lwork!, i, j
   real*8, dimension (dim**4), intent(out) :: W
   real*8, dimension (LWMAX):: WORK



   ! Consulte o espa ç o de trabalho ideal .
   JOBZ = 'V'; UPLO = 'U'
   N = dim**4 ; LDA = dim**4 ; lwork = -1

   call dsyev ( JOBZ , UPLO , N , A , LDA , W , WORK , LWORK , INFO )

   V = A

   LWORK = MIN ( LWMAX , INT ( WORK (1) ) )

   ! Resolvendo o problema de autovalores .
   call dsyev ( JOBZ , UPLO , N , V , LDA , W , WORK , LWORK , INFO )

   ! Checando a converg ê ncia .
   if ( info > 0 ) then
   write (* ,*) 'O algoritmo falhou em encontrar os autovalores'
   stop
   end if

   ! Print autovalores .
   ! call print_matrix ( 'eigenvalues' , n , 1 , w , 1)

   !read(*,*)



   ! ! Print autovetores .
   ! call print_matrix ( 'Eigenvectors' , n , n , A , LDA )

   !call print_matrix(A,dim**4,dim**4)
   !call print_matrix(D,dim**4,dim**4)

      

end subroutine

subroutine decompSpectral ( dim , A , P , matrixDiagonal )
    implicit none
    ! Inputs
    integer ,intent (in):: dim
    real*8 ,dimension (dim,dim), intent (in) :: P , A
    real*8 ,dimension (dim,dim), intent (out) :: matrixDiagonal
    ! Variables
    real*8 ,dimension (dim,dim) :: AP
   
    AP = matmul (A , P)

    matrixDiagonal = matmul ( transpose ( P ) , AP )
end subroutine

subroutine magnetization_diag(W,Z,T,dim,s,V,m)
   implicit none
   integer, intent(in):: dim
   real*8, dimension(dim**4,dim**4), intent(in):: s,V
   real*8, dimension(dim**4,dim**4):: m_prime
   real*8, dimension(dim**4), intent(in):: W
   real*8, intent(in):: T, Z
   real*8, intent(out):: m
   integer :: i
   real*8:: b

   b = 1.d0/T
   m = 0
   m_prime = 0

   m_prime = matmul(transpose(V),matmul(s,V))
   !call print_matrix(m_prime,dim**4,dim**4)
   !read(*,*)

   ! call print_matrix(m_prime,16,16)
   ! read(*,*)

   do i = 1, dim**4

      m = m + m_prime(i,i)*(dexp(-b*(W(i))))

   end do

   m = m/Z

end subroutine

end module QUANTICO
