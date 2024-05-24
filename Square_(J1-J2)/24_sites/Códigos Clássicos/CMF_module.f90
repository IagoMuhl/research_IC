module CMF
   implicit none
   integer, parameter:: db=8
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 24
   integer, parameter:: maxConfig= minConfig**num_sites
   real(kind=db), parameter:: J1=1.d0

contains

   subroutine base(s,s_sub)

      implicit none
      integer:: k
      integer:: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16
      integer:: s17,s18,s19,s20,s21,s22,s23,s24
      integer, intent(out):: s(maxConfig,num_sites)
      integer, intent(out):: s_sub(maxConfig,8)

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
                              do s9 =-1,1,2
                                 do s10 =-1,1,2
                                    do s11 =-1,1,2
                                       do s12 =-1,1,2
                                          do s13 =-1,1,2
                                             do s14 =-1,1,2
                                                do s15 =-1,1,2
                                                   do s16 =-1,1,2
                                                      do s17 =-1,1,2
                                                         do s18 =-1,1,2
                                                            do s19 =-1,1,2
                                                               do s20 =-1,1,2
                                                                  do s21 =-1,1,2
                                                                     do s22 =-1,1,2
                                                                        do s23 =-1,1,2
                                                                           do s24 =-1,1,2

                                                                              s(k,1)=s1; s_sub(k,1)=s1
                                                                              s(k,2)=s2; s_sub(k,2)=s2
                                                                              s(k,3)=s3; s_sub(k,3)=s3
                                                                              s(k,4)=s4; s_sub(k,4)=s4
                                                                              s(k,5)=s5; s_sub(k,5)=s5
                                                                              s(k,6)=s6; s_sub(k,6)=s6
                                                                              s(k,7)=s7
                                                                              s(k,8)=s8
                                                                              s(k,9)=s9
                                                                              s(k,10)=s10
                                                                              s(k,11)=s11
                                                                              s(k,12)=s12
                                                                              s(k,13)=s13
                                                                              s(k,14)=s14
                                                                              s(k,15)=s15
                                                                              s(k,16)=s16
                                                                              s(k,17)=s17; s_sub(k,7)=s17
                                                                              s(k,18)=s18; s_sub(k,8)=s18
                                                                              s(k,19)=s19
                                                                              s(k,20)=s20
                                                                              s(k,21)=s21
                                                                              s(k,22)=s22
                                                                              s(k,23)=s23
                                                                              s(k,24)=s24

                                                                              k=k+1

                                                                           end do
                                                                        end do
                                                                     end do
                                                                  end do
                                                               end do
                                                            end do
                                                         end do
                                                      end do
                                                   end do
                                                end do
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do


   end subroutine


   subroutine HAM_INTRA(J2,s,H_intra,N)
      integer, dimension(maxConfig,num_sites), intent(in):: s
      integer :: i
      real(kind=db), intent(in):: J2
      real(kind=db), dimension(maxConfig), intent(out):: H_intra
      real(kind=db), dimension(maxConfig,6), intent(out):: N

      ! N = 0

      !---------------------HAMILTONIANO INTRA-----------------------------
      do i = 1, maxConfig

         H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,2)*s(i,3) + s(i,3)*s(i,4) + s(i,4)*s(i,5) + s(i,5)*s(i,6) + &
         &  s(i,16)*s(i,17) + s(i,17)*s(i,18) + s(i,18)*s(i,19) + s(i,19)*s(i,20) + s(i,20)*s(i,7) + s(i,15)*s(i,24) + &
         &  s(i,24)*s(i,23) + s(i,23)*s(i,22) + s(i,22)*s(i,21) + s(i,21)*s(i,8) + s(i,14)*s(i,13) + s(i,13)*s(i,12) + &
         &  s(i,12)*s(i,11) + s(i,11)*s(i,10) + s(i,10)*s(i,9) + s(i,1)*s(i,16) + s(i,2)*s(i,17) + s(i,3)*s(i,18) + &
         &  s(i,4)*s(i,19) + s(i,5)*s(i,20) + s(i,6)*s(i,7) + s(i,16)*s(i,15) + s(i,17)*s(i,24) + s(i,18)*s(i,23) + &
         &  s(i,19)*s(i,22) + s(i,20)*s(i,21) + s(i,7)*s(i,8) + s(i,15)*s(i,14) + s(i,24)*s(i,13) + s(i,23)*s(i,12) + &
         &  s(i,22)*s(i,11) + s(i,21)*s(i,10) + s(i,8)*s(i,9)) + &

         &  J2*(s(i,1)*s(i,17) + s(i,2)*s(i,18) + s(i,3)*s(i,19) + s(i,4)*s(i,20) + s(i,5)*s(i,7) + s(i,16)*s(i,24) + &
         &  s(i,17)*s(i,23) + s(i,18)*s(i,22) + s(i,19)*s(i,21) + s(i,20)*s(i,8) + s(i,15)*s(i,13) + s(i,24)*s(i,12) + &
         &  s(i,23)*s(i,11) + s(i,22)*s(i,10) + s(i,21)*s(i,9) + s(i,2)*s(i,16) + s(i,3)*s(i,17) + s(i,4)*s(i,18) + &
         &  s(i,5)*s(i,19) + s(i,6)*s(i,20) + s(i,17)*s(i,15) + s(i,18)*s(i,24) + s(i,19)*s(i,23) + s(i,20)*s(i,22) + &
         &  s(i,7)*s(i,21) + s(i,24)*s(i,14) + s(i,23)*s(i,13) + s(i,22)*s(i,12) + s(i,21)*s(i,11) + s(i,8)*s(i,10))

         !---------------------HAMILTONIANO INTRA-----------------------------

         N(i,1) = 2*J1*(s(i,6)+s(i,14)) + J2*(s(i,1)+s(i,5)+s(i,7)+s(i,9)+s(i,13)+s(i,15))
         N(i,2) = J1*(s(i,13)+s(i,5)+s(i,7)+s(i,15)) + J2*(2*s(i,6)+s(i,4)+s(i,8)+s(i,12)+2*s(i,14)+s(i,16))
         N(i,3) = J1*(s(i,4)+s(i,12)) + J2*(s(i,3)+s(i,5)+s(i,11)+s(i,13))
         N(i,4) = J1*(s(i,3)+s(i,11)) + J2*(s(i,2)+s(i,4)+s(i,10)+s(i,12))
         N(i,5) = J1*(s(i,2)+s(i,8)+s(i,10)+s(i,16)) + J2*(2*s(i,1)+s(i,3)+s(i,7)+s(i,11)+2*s(i,9)+s(i,15))
         N(i,6) = 2*J1*(s(i,1)+s(i,9)) + J2*(s(i,2)+s(i,6)+s(i,8)+s(i,10)+s(i,14)+s(i,16))



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

         ! read*,

      end do


      m = m/Z

      ! print*, m, Z, inf


      !  print*, Z, m, T
      ! T = T + 0.1d0

      !end do
   end subroutine



   subroutine mag_vetor(state,m)
      implicit none
      character(len=*), intent(in):: state
      real(kind=db), dimension(8), intent(out):: m
      real*8:: m_fe, m_af
      integer:: i

      m_fe = 1.d0;
      m_af = -1.d0

      m = 0.d0

      select case (state)

       case ('AF')

         do i = 1,7,2
            m(i) = m_fe
         enddo

         do i = 2,8,2
            m(i) = m_af
         enddo

         ! do i = 1,11,2
         !    m(i) = m_fe
         ! enddo

         ! do i = 2,12,2
         !    m(i) = m_af
         ! enddo

       case ('2AF')
         m_fe = 0.84719110987579493
         m_af = 0.65197042353076307

         do i = 1,7,2
            m(i) = m_fe
         enddo

         do i = 2,8,2
            m(i) = m_af
         enddo

         ! do i = 1,11,2
         !    m(i) = m_fe
         ! enddo

         ! do i = 2,12,2
         !    m(i) = m_af
         ! enddo


       case ('PM')

         do i = 1,8
            m(i) = m_fe
         enddo

         ! do i = 1,12
         !    m(i) = m_fe
         ! enddo


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



   subroutine Ham_inter_state(J2,N,m,H_inter)
   !subroutine Ham_inter_state(J2,s,m,H_inter)
      implicit none
      !integer, dimension(maxConfig,num_sites), intent(in):: s
      real(kind=db), dimension(maxConfig,6), intent(in):: N
      integer :: i
      real(kind=db), intent(in):: J2
      real(kind=db), dimension(8), intent(in):: m
      ! integer, dimension(maxConfig), intent(out):: N11,N12,N21,N22,N31,N32,N41,N42,N51,N52,N61,N62
      real(kind=db), dimension(maxConfig), intent(out):: H_inter
      !character(len=3), intent(in) :: state

      H_inter = 0.d0

      do i = 1, maxConfig


         ! H_inter(i) = J1*(2*m(1)*(s(i,6)+s(i,14)) &
         ! &        +        (m(2)*(s(i,13)+s(i,5)+s(i,7)+s(i,15))) &
         ! &        +        (m(3)*(s(i,4)+s(i,12))) &
         ! &        +        (m(4)*(s(i,3)+s(i,11))) &
         ! &        +        (m(5)*(s(i,2)+s(i,10)+s(i,8)+s(i,16))) &
         ! &        +        (2*m(6)*(s(i,1)+s(i,9))) &
         ! &        -      (4*m(1)*m(6)+4*m(2)*m(5)+2*m(3)*m(4))) + &


         ! & J2*(m(1)*(s(i,1)+s(i,5)+s(i,7)+s(i,9)+s(i,13)+s(i,15)-m(1)) + &
         ! &     m(2)*(2*s(i,6)+s(i,4)+s(i,8)+2*s(i,14)+s(i,12)+s(i,16)-m(2)) + &
         ! &     m(3)*(s(i,3)+s(i,5)+s(i,11)+s(i,13)-m(3)) + &
         ! &     m(4)*(s(i,2)+s(i,4)+s(i,10)+s(i,12)-m(4)) + &
         ! &     m(5)*(2*s(i,1)+s(i,3)+s(i,7)+2*s(i,9)+s(i,11)+s(i,15)-m(5)) + &
         ! &     m(6)*(s(i,2)+s(i,6)+s(i,8)+s(i,10)+s(i,14)+s(i,16)-m(6)) - &
         ! &     2*(2*m(1)*m(5)+2*m(2)*m(6)+m(2)*m(4)+m(3)*m(5)))

         H_inter(i) = m(1)*(N(i,1) - J2*m(1)) + m(2)*(N(i,2) - J2*m(2)) &
         &  + m(3)*(N(i,3) - J2*m(3)) + m(4)*(N(i,4) - J2*m(4)) &
         &  + m(5)*(N(i,5) - J2*m(5)) + m(6)*(N(i,6) - J2*m(6)) &
         &  - 2*(J1*(2*m(1)*m(6) + 2*m(2)*m(5) + m(3)*m(4)) &
         &  + J2*(2*m(1)*m(5) + 2*m(2)*m(6) + m(2)*m(4) + m(3)*m(5)))

      enddo


   end subroutine

   subroutine order_parameter(state,m,m_order)
      real(kind=db), dimension(8), intent(in):: m
      character(len=3), intent(in):: state
      real(kind=db), intent(out):: m_order
      integer:: i

      m_order = 0.d0

      select case(state)

       case('AF')

         do i = 1,7,2
            m_order = m_order + (m(i) - m(i+1))
         enddo

         ! m_order = (m(1)-2*m(2)+m(3)-m(4)+2*m(5)-m(6)+2*m(7)-2*m(8))

         ! do i = 7,11,2
         !    m_order = m_order + (m(i+1) - m(i))
         ! enddo

         !m_order = abs(3*m_order/num_sites)

         m_order = abs(m_order/8)

       case('2AF')

         do i = 1,7,2
            m_order = m_order + (m(i) - m(i+1))
         enddo

         ! m_order = (m(1)-2*m(2)+m(3)-m(4)+2*m(5)-m(6)+2*m(7)-2*m(8))

         !m_order = abs(2*m_order/num_sites)

         m_order = abs(m_order/8)

       case('PM')

         do i = 1,8
            m_order = m_order + abs(m(i))
         enddo

         ! m_order = (m(1)+2*m(2)+m(3)+m(4)+2*m(5)+m(6)+2*m(7)+2*m(8))

         m_order = abs(m_order/8)

       case default
         print*, 'ERRO'

      end select

   end subroutine

   ! subroutine simplify(J2,s,N,S_sub)
   !    integer, dimension(maxConfig,num_sites), intent(in):: s
   !    real*8:: J2
   !    real(kind=db), dimension(maxConfig,6), intent(out):: N
   !    integer, dimension(maxConfig,8), intent(out):: S_sub
   !    integer:: i,j

   !    S_sub = 0

   !    do j = 1, 8
         
   !    do i = 1, maxConfig

   !    N(i,1) = 2*J1*(s(i,6)+s(i,14)) + J2*(s(i,1)+s(i,5)+s(i,7)+s(i,9)+s(i,13)+s(i,15))
   !    N(i,2) = J1*(s(i,13)+s(i,5)+s(i,7)+s(i,15)) + J2*(2*s(i,6)+s(i,4)+s(i,8)+s(i,12)+2*s(i,14)+s(i,16))
   !    N(i,3) = J1*(s(i,4)+s(i,12)) + J2*(s(i,3)+s(i,5)+s(i,11)+s(i,13))
   !    N(i,4) = J1*(s(i,3)+s(i,11)) + J2*(s(i,2)+s(i,4)+s(i,10)+s(i,12))
   !    N(i,5) = J1*(s(i,2)+s(i,8)+s(i,10)+s(i,16)) + J2*(2*s(i,1)+s(i,3)+s(i,7)+s(i,11)+2*s(i,9)+s(i,15))
   !    N(i,6) = 2*J1*(s(i,1)+s(i,9)) + J2*(s(i,2)+s(i,6)+s(i,8)+s(i,10)+s(i,14)+s(i,16))

   !    S_sub(i,j) = S_sub(i,j) + s(i,j)

   !    enddo
   !    enddo





   ! end subroutine


end module CMF
