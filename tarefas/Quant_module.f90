module QUANTICO
   implicit none
   integer, parameter:: J1=1
   contains

   subroutine print_matrix(A,m)
      implicit none
      integer, intent(in):: m
      real*8, dimension(m,m), intent(in):: A
      integer:: i, j
      do i=1,m
         write(*,20)(A(i,j),j=1,m) !para numeros inteiros (4(I3))
         20 format (16(F3.0))
      enddo
   end subroutine
   
   subroutine tensorial(A,B,m,E)
      implicit none
      integer, intent(in):: m
      real*8, dimension(m,m), intent(in):: A
      real*8, dimension(m,m), intent(in):: B
      real*8, dimension(m**2,m**2), intent(out):: E
      integer:: i, j, k, l, alfa, beta
   
      do i = 1,m
         do j = 1,m
            do k = 1,m
               do l = 1,m
   
                  alfa = m*(i-1) + k
                  beta = m*(j-1) + l
   
                  E(alfa,beta) = A(i,j) * B(k,l)
   
               end do
            end do
         end do
      end do
   end subroutine
   
   subroutine partition(H,T,m,Z)
      real*8, intent(out):: Z
      integer, intent(in):: m
      real*8, dimension(m**4,m**4), intent(in):: H
      real*8, intent(in):: T
      real*8:: b
      integer :: i
      b = 1.d0/T
      Z = 0.d0
   
   
      do i = 1, m**4
         Z = Z + (dexp(-b*(H(i,i))))
      end do
      
   end subroutine
   
   subroutine magnetization(H,Z,T,m,s1,mag)
      real*8:: b
      integer, intent(in):: m
      real*8, dimension(m**4,m**4), intent(in):: H,s1
      real*8, intent(out):: mag
      real*8, intent(in):: T, Z
      integer :: i
     
      
   
      b = 1.d0/T
      mag = 0

      do i = 1, m**4
         PRINT*, mag
         mag = mag + s1(i,i)*(dexp(-b*(H(i,i))))
         !READ(*,*)
         !print*, 'HA MULQKE',i
      end do
   
      mag = mag/Z
 
   end subroutine

   subroutine mag_vetor(state,sigma,mag)
      implicit none
      character(len=*), intent(in):: state
      real*8, intent(in):: sigma
      real*8, dimension(4), intent(out):: mag


      if ( state=='AF' ) then
         mag = [sigma, -sigma, -sigma, sigma]
      else if (state=='SAF') then
         mag = [sigma, -sigma, sigma, -sigma]
      else if (state=='PM') then
         mag = 0.d0
      else if (state=='SD') then
         mag = [sigma, sigma, -sigma, -sigma]
      end if

   end subroutine

   subroutine Ham_inter_state(state,J2,s1,s2,s3,s4,m_guess,Id_4,H_inter)
      implicit none
            real*8, dimension(16,16), intent(in):: s1,s2,s3,s4, Id_4
            integer :: i
            real*8, intent(in):: J2
            real*8, dimension(4), intent(in):: m_guess
            character(len=3), intent(in):: state
            real*8, dimension(16,16), intent(out):: H_inter

      select case (state)
      case ('AF')
         do i = 1, 16
            H_inter(i,i) = (-2*J1 + 3*J2)*(s1(i,i)-s2(i,i)-s3(i,i)+s4(i,i)-2*m_guess(1)*Id_4(i,i))*(m_guess(1)*Id_4(i,i))
         enddo
      case ('SAF')
         do i = 1, 16
            H_inter(i,i) = (-3*J2)*(s1(i,i)-s2(i,i)+s3(i,i)-s4(i,i)-2*m_guess(1)*Id_4(i,i))*(m_guess(1)*Id_4(i,i))
         enddo
      case ('SD')
         do i = 1, 16
            H_inter(i,i) = (-2*J1+J2)*(s1(i,i)+s2(i,i)-s3(i,i)-s4(i,i)-2*m_guess(1)*Id_4(i,i))*(m_guess(1)*Id_4(i,i))
         enddo
      case ('PM')
         H_inter = 0
      case default
         print *, 'State inválido'
      end select
      
         end subroutine

end module QUANTICO
