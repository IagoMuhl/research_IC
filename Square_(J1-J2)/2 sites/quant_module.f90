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
      real*8, dimension(dim*dim,dim*dim), intent(out):: E
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
      real*8, dimension(dim), intent(in):: W
      real*8, intent(in):: T
      real*8:: b
      integer :: i
      b = 1.d0/T
      Z = 0.d0

      do i = 1, dim
         Z = Z + (dexp(-b*(W(i))))
      end do

end subroutine

! subroutine magnetization(W,Z,T,dim,s,m)
!       implicit none
!       integer, intent(in):: dim
!       real*8, dimension(4,4), intent(in):: s
!       real*8, dimension(4):: W
!       real*8, intent(in):: T, Z
!       real*8, intent(out):: m
!       integer :: i
!       real*8:: b

!       b = 1.d0/T
!       m = 0

!       do i = 1, dim*2
!          m = m + s(i,i)*(dexp(-b*(W(i))))
!       end do

!       m = m/Z

! end subroutine

subroutine mag_vetor(state,m1,m2,m,m_order)
   implicit none
   character(len=*), intent(in):: state
   real*8, intent(in):: m1,m2
   real*8, dimension(2), intent(out):: m
   real*8, intent(out):: m_order


   select case (state)

    case ('AF')
      m = [m1, m2]


      m_order = abs(m(1)-m(2))/2.d0

      case ('2AF')
        m = [m1, m2]

  
        m_order = abs(m(1)-m(2))/2.d0

   !  case ('SAF')
   !    m = [m1, m2, m1, m2]


   !    m_order = abs(m(sigma1) - m(sigma2))/2.d0

    case ('PM')
      m = [m1,m2]


      m_order = abs(m(1)+m(2))/2.d0

   !  case ('SD')
   !    m = [m1, m1, m2, m2]


   !    m_order = abs(m(sigma1) - m(sigma2))/2.d0

   !  case ('ZIG')
   !    m = [m1, m2, m3, m4]

   !    m_order = abs(m(1) - m(2) - m(3) - m(4))/4.d0

    case default
      write(*,*) 'Inaccurate State'
   end select
end subroutine


subroutine Ham_inter_state(J2,s1,s2,m_guess,Id,H_inter)
   implicit none
   real*8, dimension(4,4), intent(in):: s1,s2,Id
   real*8, dimension(2):: m_guess
   real*8, intent(in):: J2
   real*8, dimension(4,4), intent(out):: H_inter


   H_inter = 0.d0


   H_inter =  3*J1*((s1*m_guess(2)-(m_guess(1)*m_guess(2)*Id)/2.d0)+(s2*m_guess(1)-(m_guess(1)*m_guess(2)*Id)/2.d0))


   H_inter = H_inter + 4*J2*((s1*m_guess(1)-(m_guess(1)*m_guess(1)*Id)/2.d0)+(s2*m_guess(2)-(m_guess(2)*m_guess(2)*Id)/2.d0))
   

   ! &  + 3*J2*(s1*m_guess(4)+s4*m_guess(1)-m_guess(1)*m_guess(4)*Id_4 &
   ! &  +   s2*m_guess(3)+s3*m_guess(2)-m_guess(2)*m_guess(3)*Id_4)

      !  H_inicial = 1.5d0
      !  H_final = 8.d0
end subroutine

! subroutine Ham_inter_state(state,J2,J3,s1,s2,s3,s4,H,m1,m2,Id_4,H_inter)
!       implicit none
!       character(len=3), intent(in):: state
!       real*8, dimension(16,16), intent(in):: s1,s2,s3,s4,Id_4
!       real*8, intent(in):: J2, J3, m1, m2 , H
!       real*8, dimension(16,16), intent(out):: H_inter

!       H_inter = 0.d0

!       ! if (H/=0.0) then

!       select case (state)
!        case ('AF')

!          H_inter = (2*J1)*(m2*(s1+s4)+m1*(s2+s3)-(2*m1*m2*Id_4))+(3*J2 + 4*J3)*(m1*(s1+s4-(m1*Id_4))+m2*(s2+s3-(m2*Id_4)))

!        case ('PM')
!          H_inter = 0
      
!       !  case default
!       !    print *, 'State inválido'
! !
!          !-H*(s1+s2+s3+s4) 
! !

! ! !
! !        case ('AF')
! !          H_inter = (-2*J1 + 3*J2 + 4*J3)*(s1-s2-s3+s4-2*m1*Id_4)*(m1)

!       !  case ('SAF')
!       !    H_inter = (-3*J2 + 4*J3)*(s1-s2+s3-s4-2*m1*Id_4)*(m1)

!       !  case ('SD')
!       !    H_inter = (-2*J1 + J2)*(s1+s2-s3-s4-2*m1*Id_4)*(m1)

!       !  case ('PM')
!       !    H_inter = 0
!        case default
!          print *, 'State inválido'
!       end select

! !

! end subroutine

subroutine chose(state,m_fe,m_af,m_guess)
real*8, intent(in):: m_fe, m_af
real*8, dimension(2),intent(out):: m_guess
character(len=3), intent(in):: state

select case(state)

case('AF')
m_guess = [m_fe,m_af]

case('2AF')
m_guess = [m_fe,m_af]

case('PM')
m_guess = [m_fe,m_fe]

! case('ZIG')
! m_guess = [m_fe,m_af,m_af,m_af]

case default
   print*, 'Estado não encontrado'
end select

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
   real*8, dimension (4,4), intent(in) :: A
   real*8, dimension(4,4), intent(out):: V
   ! !!! Para diagonalizar !!!!!
   character(len=1):: JOBZ , UPLO
   integer:: N , LDA , INFO , lwork!, i, j
   real*8, dimension (4), intent(out) :: W
   real*8, dimension (LWMAX):: WORK



   ! Consulte o espa ç o de trabalho ideal .
   JOBZ = 'V'; UPLO = 'U'
   N = 4 ; LDA = 4 ; lwork = -1

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

! subroutine decompSpectral ( dim , A , P , matrixDiagonal )
!     implicit none
!     ! Inputs
!     integer ,intent (in):: dim
!     real*8 ,dimension (dim,dim), intent (in) :: P , A
!     real*8 ,dimension (dim,dim), intent (out) :: matrixDiagonal
!     ! Variables
!     real*8 ,dimension (dim,dim) :: AP
   
!     AP = matmul (A , P)

!     matrixDiagonal = matmul ( transpose ( P ) , AP )
! end subroutine

subroutine magnetization_diag(W,Z,T,dim,s,V,m)
   implicit none
   integer, intent(in):: dim
   real*8, dimension(dim*2,dim*2), intent(in):: s,V
   real*8, dimension(dim*2,dim*2):: m_prime
   real*8, dimension(dim*2), intent(in):: W
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

   do i = 1, dim*2

      m = m + m_prime(i,i)*(dexp(-b*(W(i))))

   end do

   m = m/Z

end subroutine

subroutine int_energy(dim,W,T,Z,U)
integer, intent(in):: dim
real*8, dimension (dim), intent(in) :: W
real*8, intent(in):: T, Z
real*8, intent(out):: U
integer:: i
U = 0.d0

   do i = 1,dim

   U = U + W(i)*(dexp(-W(i)/T))

   enddo

   U = U/Z

end subroutine


subroutine especific_heat(step,T,U,Uo,Cv)
real*8, intent(in):: step, T, U, Uo
real*8, intent(out):: Cv


Cv = (U-Uo)/(T -(step/2.d0))


end subroutine

subroutine entropy(U,F,T,S)
real*8, intent(in):: U,F,T
real*8, intent(out):: S

S = (U - F)/T

end subroutine

end module QUANTICO
