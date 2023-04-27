module QUANTICO
   implicit none
   integer, parameter:: J1=1
contains

subroutine print_matrix(A,dim)
      implicit none
      integer, intent(in):: dim
      real*8, dimension(dim,dim), intent(in):: A
      integer:: i, j
      do i=1,dim
         write(*,20)(A(i,j),j=1,dim) !para numeros inteiros (4(I3))
      20   format (16(F3.0))
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

subroutine partition(H,T,dim,Z)
      implicit none
      real*8, intent(out):: Z
      integer, intent(in):: dim
      real*8, dimension(dim**4,dim**4), intent(in):: H
      real*8, intent(in):: T
      real*8:: b
      integer :: i
      b = 1.d0/T
      Z = 0.d0

      do i = 1, dim**4
         Z = Z + (dexp(-b*(H(i,i))))
      end do

end subroutine

subroutine magnetization(H,Z,T,dim,s1,m)
      implicit none
      integer, intent(in):: dim
      real*8, dimension(dim**4,dim**4), intent(in):: H,s1
      real*8, intent(in):: T, Z
      real*8, intent(out):: m
      integer :: i
      real*8:: b

      b = 1.d0/T
      m = 0

      do i = 1, dim**4
         m = m + s1(i,i)*(dexp(-b*(H(i,i))))
      end do

      m = m/Z

end subroutine

subroutine Ham_inter_state(state,J2,J3,s1,s2,s3,s4,m_guess,Id_4,H_inter)
      implicit none
      character(len=3), intent(in):: state
      real*8, dimension(16,16), intent(in):: s1,s2,s3,s4, Id_4
      real*8, intent(in):: J2, J3, m_guess
      real*8, dimension(16,16), intent(out):: H_inter

      H_inter = 0.d0

      select case (state)
       case ('AF')
         H_inter = (-2*J1 + 3*J2 + 4*J3)*(s1-s2-s3+s4-2*m_guess*Id_4)*(m_guess)

       case ('SAF')
         H_inter = (-3*J2 + 4*J3)*(s1-s2+s3-s4-2*m_guess*Id_4)*(m_guess)

       case ('SD')
         H_inter = (-2*J1 + J2)*(s1+s2-s3-s4-2*m_guess*Id_4)*(m_guess)

       case ('PM')
         H_inter = 0
       case default
         print *, 'State inválido'
      end select

end subroutine

subroutine Free_nrg(T,Z,F)
      implicit none
      real*8, intent(in):: Z,T
      real*8, intent(out):: F

      F = -T*dlog(Z)

end subroutine

SUBROUTINE diagonalization(A, V, D, N)
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: N
   REAL*8, INTENT(INOUT) :: A(N,N)
   REAL*8, INTENT(INOUT) :: V(N,N)
   REAL*8, INTENT(INOUT) :: D(N,N)
   INTEGER :: i, j, k, l
   REAL*8 :: c, s, t, tau, term1, term2

   ! Inicialização das matrizes D e V
   D = A
   V = 0.0
   DO i = 1, N
       V(i,i) = 1.0
   END DO

   ! Algoritmo QR com shift de Wilkinson
   DO l = 1, N-1
       DO k = 1, 1000*N**2 ! Número máximo de iterações permitidas
           ! Teste de convergência
           term1 = ABS(D(l,l)) + ABS(D(l+1,l+1))
           term2 = 0.0
           DO i = l+1, N
               term2 = term2 + ABS(D(i,l))
           END DO
           IF (term2 <= 1e-12*term1) EXIT
           ! Shift de Wilkinson
           tau = (D(l,l) - D(l+1,l+1))/(2.0*D(l+1,l))
           IF (tau >= 0.0) THEN
               t = 1.0/(tau + SQRT(1.0 + tau**2))
           ELSE
               t = -1.0/(-tau + SQRT(1.0 + tau**2))
           END IF
           c = 1.0/SQRT(1.0 + t**2)
           s = t*c
           ! Atualização das matrizes D e V
           DO i = l+1, N
               term1 = D(l,i)
               term2 = D(l+1,i)
               D(l,i) = c*term1 - s*term2
               D(l+1,i) = s*term1 + c*term2
           END DO
           DO i = 1, N
               term1 = V(i,l)
               term2 = V(i,l+1)
               V(i,l) = c*term1 - s*term2
               V(i,l+1) = s*term1 + c*term2
           END DO
           DO j = l+2, N
               term1 = D(j,l)
               term2 = D(j,l+1)
               D(j,l) = c*term1 - s*term2
               D(j,l+1) = s*term1 + c*term2
           END DO
       END DO
       IF (k == 1000*N**2) STOP "QR algorithm did not converge"
   END DO


   call print_matrix(A,N)
   read(*,*)

   call print_matrix(V,N)
   read(*,*)

   call print_matrix(D,N)
   read(*,*)

END SUBROUTINE


end module QUANTICO
