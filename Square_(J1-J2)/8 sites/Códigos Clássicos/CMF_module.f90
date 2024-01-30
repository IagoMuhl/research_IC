module CMF
   implicit none
   integer, parameter:: db=8
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 8
   integer, parameter:: maxConfig= minConfig**num_sites
   real(kind=db), parameter:: J1=1.d0
contains

   subroutine base(s)

      implicit none
      integer:: k
      integer:: s1,s2,s3,s4,s5,s6,s7,s8
      integer, intent(out):: s(maxConfig,num_sites)


      k=1
      !!!!!!!!!!    CONTRUÇÃO DA BASE     !!!!!!!!!!!!
      do s1 =-1,1,2
         do s2 =-1,1,2
            do s3 =-1,1,2
               do s4 =-1,1,2
                  do s5 =-1,1,2
                     do s6 =-1,1,2
                        do s7 =-1,1,2
                           do s8 =-1,1,2

                              s(k,1)=s1
                              s(k,2)=s2
                              s(k,3)=s3
                              s(k,4)=s4
                              s(k,5)=s5
                              s(k,6)=s6
                              s(k,7)=s7
                              s(k,8)=s8

                              k=k+1

                           end do
                        end do
                     end do
                  end do
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
         H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,1)*s(i,5 ) + s(i,2)*s(i,6) + s(i,5)*s(i,6) + s(i,2)*s(i,3) + &
         &      s(i,6)*s(i,7) + s(i,3)*s(i,7)+ s(i,3)*s(i,4) + s(i,4)*s(i,8) + s(i,7)*s(i,8)) + &
         &  J2*(s(i,1)*s(i,6) + s(i,2)*s(i,5) + s(i,2)*s(i,7) + s(i,6)*s(i,3) + s(i,3)*s(i,8) + s(i,4)*s(i,7))

      end do
      !---------------------HAMILTONIANO INTRA-----------------------------


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


      do i = 1, maxConfig
         Z = Z + (dexp(-b*(H(i))))
      end do

   end subroutine



   subroutine magnetization(H,Z,s,inf,T,m)
      integer, dimension(maxConfig,num_sites), intent(in) :: s
      real(kind=db):: b
      real(kind=db), dimension(maxConfig), intent(in):: H
      real(kind=db), intent(out):: m
      real(kind=db), intent(in):: T, Z
      integer, intent(in):: inf
      integer :: i
      b = 1.d0/T
      m = 0.d0



      do i = 1, maxConfig

         m = m + s(i,inf)*dexp(-b*(H(i)))

      end do

      m = m/Z


      !  print*, Z, m, T
      ! T = T + 0.1d0

      !end do
   end subroutine



   subroutine mag_vetor(state,m1,m2,m3,m4,m,m_order)
      implicit none
      character(len=*), intent(in):: state
      real(kind=db), intent(in):: m1,m2,m3,m4
      real(kind=db), dimension(4), intent(out):: m
      real(kind=db), intent(out):: m_order

      select case (state)

       case ('AF')
         m = [m1, m2, m3, m4]

         m_order = abs(m(1) - m(2) + m(3) - m(4))/4.d0

       case ('2AF')
         m = [m1, m2, m3, m4]

         m_order = abs(m(1) - m(2) + m(3) - m(4))/4.d0

       case ('PM')
         m = [m1,m3,m3,m1]

         m_order = abs(m(1) + m(2) + m(3) + m(4))/4.d0

       case default
         write(*,*) 'Inaccurate State'
      end select

   end subroutine


   subroutine F_helm(T,Z,F)
      implicit none
      real(kind=db), intent(in):: Z,T
      real(kind=db), intent(out):: F

      F = -T*dlog(Z)

   end subroutine



   subroutine print_m(state,J2,m,tol,T,n)
      real(kind=db), intent(in):: J2,m,tol,T
      character(len=*):: state
      integer:: n

      if (m<tol) then
         if (n==1) then
            print *, J2,T,state
            !write(20,*) J2, T
            n =  n + 1
         endif
      endif

   end subroutine



   subroutine Ham_inter_state(J2,s,m_guess,H_inter)
      implicit none
      integer, dimension(maxConfig,num_sites), intent(in):: s
      integer :: i
      real(kind=db), intent(in):: J2
      real(kind=db), dimension(4), intent(in):: m_guess
      real(kind=db), dimension(maxConfig), intent(out):: H_inter
      !character(len=3), intent(in) :: state

      H_inter = 0.d0

      do i = 1, maxConfig

         H_inter(i) =  J1*(s(i,1)*m_guess(4)+s(i,4)*m_guess(1)-m_guess(1)*m_guess(4)  &
         &   +   s(i,1)*m_guess(4)+s(i,5)*m_guess(1)-m_guess(1)*m_guess(4)  &
         &   +   s(i,2)*m_guess(3)+s(i,6)*m_guess(2)-m_guess(1)*m_guess(3)  &
         &   +   s(i,3)*m_guess(2)+s(i,7)*m_guess(3)-m_guess(2)*m_guess(3)  &
         &   +   s(i,4)*m_guess(1)+s(i,8)*m_guess(4)-m_guess(1)*m_guess(4)  &
         &   +   s(i,5)*m_guess(1)+s(i,8)*m_guess(4)-m_guess(1)*m_guess(4)) &


         & + J2*(s(i,1)*m_guess(3)+s(i,6)*m_guess(1)-m_guess(1)*m_guess(3)  &
         &   +   s(i,2)*m_guess(4)+s(i,5)*m_guess(2)-m_guess(2)*m_guess(4) &
         &   +   s(i,2)*m_guess(2)+s(i,7)*m_guess(2)-m_guess(2)*m_guess(2)  &
         &   +   s(i,3)*m_guess(3)+s(i,6)*m_guess(3)-m_guess(3)*m_guess(3)   &
         &   +   s(i,3)*m_guess(1)+s(i,8)*m_guess(3)-m_guess(1)*m_guess(3)   &
         &   +   s(i,4)*m_guess(2)+s(i,7)*m_guess(4)-m_guess(2)*m_guess(4)  &
         &  + 2*(s(i,1)*m_guess(1)+s(i,8)*m_guess(1)-m_guess(1)*m_guess(1)   &
         &   +   s(i,4)*m_guess(2)+s(i,5)*m_guess(2)-m_guess(2)*m_guess(2)))

      enddo


   end subroutine

   ! subroutine order_parameter(state,m,m_order)
   !    real(kind=db), dimension(minConfig), intent(in):: m
   !    character(len=3), intent(in):: state
   !    real(kind=db), intent(out):: m_order

   !    select case(state)

   !     case('AF')

   !       m_order = (m(1)-m(2))/2.d0

   !     case('PM')
   !       m_order   = (m(1)+m(2))/2.d0

   !     case default
   !       print*, 'ERRO'

   !    end select

   ! end subroutine


end module CMF
