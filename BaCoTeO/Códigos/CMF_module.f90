module CMF
   implicit none
   integer, parameter:: db=8
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 6
   integer, parameter:: maxConfig= minConfig**num_sites
   real(kind=db), parameter:: J1=1.d0

contains

   subroutine base(s)

      implicit none
      integer:: k
      integer:: s1,s2,s3,s4,s5,s6
      integer, intent(out):: s(maxConfig,num_sites)

      k=1
      !!!!!!!!!!    CONTRUÇÃO DA BASE     !!!!!!!!!!!!
      do s1 =-1,1,2
         do s2 =-1,1,2
            do s3 =-1,1,2
               do s4 =-1,1,2
                  do s5 =-1,1,2
                     do s6 =-1,1,2

                     s(k,1)=s1
                     s(k,2)=s2
                     s(k,3)=s3
                     s(k,4)=s4
                     s(k,5)=s5
                     s(k,6)=s6
                     
                     k=k+1

                     end do
                  end do
               end do
            end do
         end do
      end do                                    

   end subroutine


   subroutine HAM_INTRA(J2,J3,s,H_intra)
      integer, dimension(maxConfig,num_sites), intent(in):: s
      integer :: i
      real(kind=db), intent(in):: J2, J3
      real(kind=db), dimension(maxConfig), intent(out):: H_intra


      ! N = 0

      !---------------------HAMILTONIANO INTRA-----------------------------
      do i = 1, maxConfig

         H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,2)*s(i,3) + s(i,3)*s(i,4) + s(i,4)*s(i,5) + s(i,5)*s(i,6) + s(i,6)*s(i,1)) &
         & + J2*(s(i,1)*s(i,3) + s(i,3)*s(i,5) + s(i,5)*s(i,1) + s(i,2)*s(i,4) + s(i,4)*s(i,6) + s(i,6)*s(i,2)) &
         & + J3*(s(i,1)*s(i,4) + s(i,2)*s(i,5) + s(i,3)*s(i,6))

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

         m = m + (s(i,inf)*dexp(-b*(H(i))))

      end do

      ! print*, m, Z
      ! read*,

      m = m/Z

   end subroutine



   subroutine mag_vetor(state,m)
      implicit none
      character(len=*), intent(in):: state
      real(kind=db), dimension(6), intent(out):: m
      real*8:: up, down


      up = 1.d0
      down = -1.d0

      m = 0.d0

      select case (state)

      case('4')

         m = [down, up, up, down, up, up]

      case ('5')

         m = [up, down, down, up, down, down]

      case ('6')

         m = [down, down, up, up, up, up]

      case ('pm')

         m = [up, up, up, up, up, up]

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

   subroutine Inner_energy(H,Z,T,U)
      implicit none
      
      real(kind=db), intent(in):: Z,T,H(maxConfig)
      real(kind=db), intent(out):: U
      real(kind=db):: b
      integer:: i
      b = 1.d0/T
      U = 0.d0


      do i = 1, maxConfig

         U = U + H(i)*dexp(-b*(H(i)))

      end do
      
      U = U/Z

   end subroutine



   subroutine print_m(state,J2,m,tol,T,N)
      real(kind=db), intent(in):: J2,m,tol,T
      character(len=*):: state
      integer:: N

      if (m<tol) then
         if (N==1) then
            print *, J2,T,state
            !write(20,*) J2, T
            N =  N + 1
         endif
      endif

   end subroutine


   subroutine Ham_inter(state,J2,J3,s,m,m_prime,H_inter)
      implicit none
      integer, intent(in):: s(maxConfig,num_sites)
      character(len=*), intent(in):: state
      real(kind=db), intent(in):: J2,J3
      real(kind=db), dimension(6), intent(in):: m,m_prime
      real(kind=db), dimension(maxConfig), intent(out):: H_inter
      integer :: i

      H_inter = 0.d0

      select case (state)

      case('4')

      do i = 1, maxConfig

         H_inter(i) = J1*((s(i,1)-m(1)/2)*(m_prime(4)) + (s(i,2)-m(2)/2)*(m(6)) + (s(i,3)-m(3)/2)*(m(5)) &
         & + (s(i,4)-m(4)/2)*(m_prime(1)) + (s(i,5)-m(5)/2)*(m(3)) + (s(i,6)-m(6)/2)*(m(2))) &
         & + J2*((s(i,1)-m(1)/2)*(m(2)+m_prime(5)+m_prime(3)+m(6)) + (s(i,2)-m(2)/2)*(m_prime(4)+m(1)+2*m(5)) &
         & + (s(i,3)-m(3)/2)*(m(4)+m_prime(1)+2*m(6)) + (s(i,4)-m(4)/2)*(m(5)+m_prime(2)+m_prime(6)+m(3)) &
         & + (s(i,5)-m(5)/2)*(m(4)+m_prime(1)+2*m(2)) + (s(i,6)-m(6)/2)*(m(1)+m_prime(4)+2*m(3))) &
         & + J3*((s(i,1)-m(1)/2)*(2*m(1)) + (s(i,2)-m(2)/2)*(m(6)+m_prime(3)) + (s(i,3)-m(3)/2)*(m(5)+m_prime(2)) &
         & + (s(i,4)-m(4)/2)*(2*m(4)) + (s(i,5)-m(5)/2)*(m(3)+m_prime(6)) + (s(i,6)-m(6)/2)*(m(2)+m_prime(5)))


      enddo

      case('5')

      do i = 1, maxConfig

         H_inter(i) = J1*((s(i,1)*m(4)+s(i,4)*m(1)-m(1)*m(4))+(s(i,2)*m_prime(5)-(m(2)*m_prime(5)/2)) &
         & + (s(i,3)*m_prime(6)-(m(3)*m_prime(6))/2)+(s(i,5)*m_prime(2)-(m(5)*m_prime(2))/2) &
         & + (s(i,6)*m_prime(3)-(m(6)*m_prime(3))/2)) + J2*((s(i,1)-m(1)/2)*(m(3)+m_prime(3)+m(5)+m_prime(5)) &
         & + (s(i,2)-m(2)/2)*(m(4)+m_prime(4)+2*m_prime(6)) + (s(i,3)-m(3)/2)*(m(1)+m_prime(1)+2*m_prime(5)) &
         & + (s(i,5)-m(5)/2)*(m(1)+m_prime(1)+2*m_prime(3)) + (s(i,6)-m(6)/2)*(m(4)+m_prime(4)+2*m_prime(2)) &
         & + (s(i,4)-m(4)/2)*(m(2)+m_prime(2)+m(6)+m_prime(6))) + J3*((s(i,1)-m(1)/2)*(m_prime(2)+m_prime(6)) & 
         & + (s(i,2)-m(2)/2)*(m(3)+m_prime(1)) + (s(i,3)-m(3)/2)*(m(2)+m_prime(4)) &
         & + (s(i,4)-m(4)/2)*(m_prime(3)+m_prime(5)) + (s(i,5)-m(5)/2)*(m(6)+m_prime(4)) & 
         & + (s(i,6)-m(6)/2)*(m(5)+m_prime(1)))

      enddo

      case('6')

      do i = 1, maxConfig

         H_inter(i) = J1*((s(i,1)*m(4)+s(i,4)*m(1)-m(1)*m(4))+(s(i,2)*m_prime(5)-(m(2)*m_prime(5)/2)) &
         & + (s(i,3)*m_prime(6)-(m(3)*m_prime(6))/2)+(s(i,5)*m_prime(2)-(m(5)*m_prime(2))/2) &
         & + (s(i,6)*m_prime(3)-(m(6)*m_prime(3))/2)) + J2*((s(i,1)-m(1)/2)*(m(3)+m_prime(3)+m(5)+m_prime(5)) &
         & + (s(i,2)-m(2)/2)*(m(4)+m_prime(4)+2*m_prime(6)) + (s(i,3)-m(3)/2)*(m(1)+m_prime(1)+2*m_prime(5)) &
         & + (s(i,5)-m(5)/2)*(m(1)+m_prime(1)+2*m_prime(3)) + (s(i,6)-m(6)/2)*(m(4)+m_prime(4)+2*m_prime(2)) &
         & + (s(i,4)-m(4)/2)*(m(2)+m_prime(2)+m(6)+m_prime(6))) + J3*((s(i,1)-m(1)/2)*(m_prime(2)+m_prime(6)) & 
         & + (s(i,2)-m(2)/2)*(m(3)+m_prime(1)) + (s(i,3)-m(3)/2)*(m(2)+m_prime(4)) &
         & + (s(i,4)-m(4)/2)*(m_prime(3)+m_prime(5)) + (s(i,5)-m(5)/2)*(m(6)+m_prime(4)) & 
         & + (s(i,6)-m(6)/2)*(m(5)+m_prime(1)))

      enddo

      case('pm')

      do i = 1, maxConfig

         H_inter(i) = (s(i,1)+s(i,2)+s(i,3)+s(i,4)+s(i,5)+s(i,6)-3*m(1))*(J1*m(1) + 4*J2*m(1) + 2*J3*m(1))

      enddo

      end select

   end subroutine

      subroutine Ham_inter_fase_4(J2,J3,s,m,m_prime,H_inter,H_inter_prime)
      implicit none
      integer, intent(in):: s(maxConfig,num_sites)
      real(kind=db), intent(in):: J2,J3
      real(kind=db), dimension(6), intent(in):: m,m_prime
      real(kind=db), dimension(maxConfig), intent(out):: H_inter, H_inter_prime
      integer :: i

      H_inter = 0.d0; H_inter_prime = 0.d0

      do i = 1, maxConfig

      H_inter(i) = J1*((s(i,1)-m(1)/2)*(m_prime(4)) + (s(i,2)-m(2)/2)*(m(6)) + (s(i,3)-m(3)/2)*(m(5)) &
      & + (s(i,4)-m(4)/2)*(m_prime(1)) + (s(i,5)-m(5)/2)*(m(3)) + (s(i,6)-m(6)/2)*(m(2))) &
      & + J2*((s(i,1)-m(1)/2)*(m(2)+m_prime(5)+m_prime(3)+m(6)) + (s(i,2)-m(2)/2)*(m_prime(4)+m(1)+2*m(5)) &
      & + (s(i,3)-m(3)/2)*(m(4)+m_prime(1)+2*m(6)) + (s(i,4)-m(4)/2)*(m(5)+m_prime(2)+m_prime(6)+m(3)) &
      & + (s(i,5)-m(5)/2)*(m(4)+m_prime(1)+2*m(2)) + (s(i,6)-m(6)/2)*(m(1)+m_prime(4)+2*m(3))) &
      & + J3*((s(i,1)-m(1)/2)*(2*m(1)) + (s(i,2)-m(2)/2)*(m(6)+m_prime(3)) + (s(i,3)-m(3)/2)*(m(5)+m_prime(2)) &
      & + (s(i,4)-m(4)/2)*(2*m(4)) + (s(i,5)-m(5)/2)*(m(3)+m_prime(6)) + (s(i,6)-m(6)/2)*(m(2)+m_prime(5)))

      ! H_inter(i) = J1*(m_prime(1)*(s(i,1)+s(i,4)-m(1)) + m(2)*(s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m(2))) &
      ! & +  J2*(m_prime(1)*(2*s(i,1)+2*s(i,4)-2*m(1)+s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m(2)) &
      ! & + m(1)*(s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m(2)) + 2*m(2)*(s(i,1)+s(i,4)-m(1)+s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m(2))) &
      ! & + J3*(2*m(1)*(s(i,1)+s(i,4)-m(1)) + (m_prime(1) + m(2))*(s(i,2)+s(i,3)+s(i,5)+s(i,6)-2*m(2)))

      H_inter_prime(i) = J1*((s(i,1)-m_prime(1)/2)*(m(4)) + (s(i,2)-m_prime(2)/2)*(m(4)) + (s(i,3)-m_prime(3)/2)*(m(1)) &
      & + (s(i,4)-m_prime(4)/2)*(m(1)) + (s(i,5)-m_prime(5)/2)*(m(1)) + (s(i,6)-m_prime(6)/2)*(m(4))) &
      & + J2*((s(i,1)-m_prime(1)/2)*(m(4)+m(5)+m(3)+m(4)) + (s(i,2)-m_prime(2)/2)*(m(4)+m(5)+m(3)+m(1)) &
      & + (s(i,3)-m_prime(3)/2)*(m(4)+m(2)+m(6)+m(1)) + (s(i,4)-m_prime(4)/2)*(m(1)+m(2)+m(6)+m(1)) &
      & + (s(i,5)-m_prime(5)/2)*(m(1)+m(2)+m(6)+m(4)) + (s(i,6)-m_prime(6)/2)*(m(1)+m(5)+m(3)+m(4))) &
      & + J3*((s(i,1)-m_prime(1)/2)*(m(5)+m(3)) + (s(i,2)-m_prime(2)/2)*(m(2)+m(3)) + (s(i,3)-m_prime(3)/2)*(m(2)+m(3)) &
      & + (s(i,4)-m_prime(4)/2)*(m(6)+m(2)) + (s(i,5)-m_prime(5)/2)*(m(5)+m(6)) + (s(i,6)-m_prime(6)/2)*(m(5)+m(6)))

      ! H_inter_prime(i) = (s(i,1)+s(i,2)+s(i,3)+s(i,4)+s(i,5)+s(i,6)-3*m_prime(1))*(J1*m(1) + J2*2*(m(1)+m(2)) + J3*2*m(2))

      enddo

   end subroutine



   subroutine order_parameter(state,m,m_order)
      real(kind=db), dimension(6), intent(in):: m
      character(len=3), intent(in):: state
      real(kind=db), intent(out):: m_order
      real*8:: m_a,m_b

      m_order = 0.d0; m_a = 0.d0; m_a = 0.d0

      select case(state)

      case('5')

         m_a = m(1) + m(4)

         m_b = m(2) + m(3) + m(5) + m(6)

         m_order = abs(m_a - m_b)/num_sites
      
      case('6')

         m_a = m(3) + m(4) + m(5) + m(6)

         m_b = m(1) + m(2)

         m_order = abs(m_a - m_b)/num_sites

      case('4')

         m_a = m(2) + m(3) + m(5) + m(6)

         m_b = m(1) + m(4)

         m_order = abs(m_a - m_b)/num_sites

       case('pm')

         m_order = 0.d0

         ! m_order = (m(1)+m(2)+m(3)+m(4)+m(5)+m(6)+m(7)+m(8)+2*m(1)+2*(10))


       case default
         print*, 'ERRO'

      end select

   end subroutine


end module CMF
