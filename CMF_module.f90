module CMF
   implicit none
   integer, parameter:: db=8
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   real(kind=db), parameter:: J1=1.d0
contains

   subroutine base(s)

      implicit none
      integer:: k
      integer:: s1,s2,s3,s4
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

                  k=k+1

               end do
            end do
         end do
      end do


   end subroutine

   subroutine HAM_INTRA(J2,s,H_intra)
      integer, dimension(maxConfig,num_sites), intent(in):: s
      integer :: i
      real(kind=db), intent(in):: J2
      real(kind=db), dimension(maxConfig), intent(out):: H_intra

      !---------------------HAMILTONIANO INTRA-----------------------------
      do i = 1, maxConfig
         H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,1)*s(i,3) + s(i,4)*s(i,3) + s(i,2)*s(i,4)) &
         & + J2*(s(i,1)*s(i,4) + s(i,3)*s(i,2))

      end do
      !---------------------HAMILTONIANO INTRA-----------------------------


   end subroutine

   subroutine HAM_INTER(J2,J3,s,m_guess,H_inter)
      integer, dimension(maxConfig,num_sites), intent(in):: s
      integer :: i
      real(kind=db), intent(in):: J2,J3
      real(kind=db), dimension(num_sites), intent(in):: m_guess
      real(kind=db), dimension(maxConfig), intent(out):: H_inter



      do i = 1, maxConfig
         H_inter(i) = 4*(J3*sumJ3()) + 3*(J2*sumJ2()) + J1*sumJ1()

      end do
   contains
      real(kind=db) function sumJ1()
         implicit none
         sumJ1 = s(i,1)*m_guess(2)-m_guess(1)*m_guess(2)/2.+ &
         & s(i,1)*m_guess(3)-m_guess(1)*m_guess(3)/2.+ &
         & s(i,2)*m_guess(4)-m_guess(2)*m_guess(4)/2.+ &
         & s(i,2)*m_guess(1)-m_guess(2)*m_guess(1)/2.+ &
         & s(i,3)*m_guess(1)-m_guess(3)*m_guess(1)/2.+ &
         & s(i,3)*m_guess(4)-m_guess(3)*m_guess(4)/2.+ &
         & s(i,4)*m_guess(2)-m_guess(4)*m_guess(2)/2.+ &
         & s(i,4)*m_guess(3)-m_guess(4)*m_guess(3)/2.
      end function
      real(kind=db) function sumJ2()
         implicit none
         sumJ2 = s(i,1)*m_guess(4)-m_guess(1)*m_guess(4)/2.+ &
         & s(i,2)*m_guess(3)-m_guess(2)*m_guess(3)/2.+ &
         & s(i,3)*m_guess(2)-m_guess(3)*m_guess(2)/2.+ &
         & s(i,4)*m_guess(1)-m_guess(4)*m_guess(1)/2.
      end function
      real(kind=db) function sumJ3()
         implicit none
         sumJ3 = s(i,1)*m_guess(1)-(m_guess(1)*m_guess(1))/2.+ &
         & s(i,2)*m_guess(2)-(m_guess(2)*m_guess(2))/2.+ &
         & s(i,3)*m_guess(3)-(m_guess(3)*m_guess(3))/2.+ &
         & s(i,4)*m_guess(4)-(m_guess(4)*m_guess(4))/2.
      end function

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
      real(kind=db), dimension(rowH,columnH), intent(in):: matrixH
      integer:: i,j
      do i = 1, rowH
         write(*,*) (matrixH(i,j),j = 1, columnH)
      end do

   end subroutine

   subroutine partition(H,T,Z)
      real(kind=db), intent(out):: Z
      real(kind=db), dimension(maxConfig), intent(in):: H
      real(kind=db), intent(in):: T
      real(kind=db):: b
      integer :: i
      b = 1.d0/T
      Z = 0.d0

      !T = 0.1d0
      ! do while (T<=0.999d0)


      do i = 1, maxConfig
         Z = Z + (dexp(-b*(H(i))))
      end do
      !print*, Z,T

      !T = T + 0.1d0

      !end do
   end subroutine

   subroutine magnetization(H,Z,s,T,m)
      integer, dimension(maxConfig,num_sites), intent(in) :: s
      real(kind=db):: b
      real(kind=db), dimension(maxConfig), intent(in):: H
      real(kind=db), intent(out):: m
      real(kind=db), intent(in):: T, Z
      integer :: i
      b = 1.d0/T
      m = 0.d0

      !do while (T<=0.999d0)


      do i = 1, maxConfig
         m = m + s(i,1)*(dexp(-b*(H(i))))

      end do

      m = m/Z


      !  print*, Z, m, T
      ! T = T + 0.1d0

      !end do
   end subroutine

   subroutine mag_vetor(state,sigma,m)
      implicit none
      character(len=*), intent(in):: state
      real(kind=db), intent(in):: sigma
      real(kind=db), dimension(num_sites), intent(out):: m

      if ( state=='AF' ) then
         m = [sigma, -sigma, -sigma, sigma]
      else if (state=='SAF') then
         m = [sigma, -sigma, sigma, -sigma]
      end if

   end subroutine

   subroutine F_helm(Z,F)
      implicit none
      real(kind=db), intent(in):: Z
      real(kind=db), intent(out):: F

      F = dlog(Z)

   end subroutine
end module CMF