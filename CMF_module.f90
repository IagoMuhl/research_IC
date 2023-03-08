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
      else if (state=='PM') then
         m = 0.d0
      else if (state=='SD') then
         m = [sigma, sigma, -sigma, -sigma]
      end if


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

   function condition(stateOne,stateTwo,J2,J2Final,J2Inicial)
      implicit none
      logical:: condition
      character(len=*):: stateOne,stateTwo
      real(kind=db), intent(in):: J2,J2Final,J2Inicial

      select case (stateOne)
       case ('AF')

         if ( J2<=J2Final) then
            condition = .true.
         else
            condition = .false.
         end if

      case ('SD')
         if ( J2<=J2Final ) then
            condition = .true.
         else
            condition = .false.
         end if

      case ('PM')
         if ( J2<=J2Final ) then
            condition = .true.
         else
            condition = .false.
         end if

      case ('SAF')
         if ( J2<=J2Final ) then
            condition = .true.
         else
            condition = .false.
         end if

      case default
         print *, 'State inválido'

      end select

      select case (stateTwo)
      case ('AF')

        if ( J2>=J2Inicial) then
           condition = .true.
        else
           condition = .false.
        end if

     case ('SD')
        if ( J2>=J2Inicial ) then
           condition = .true.
        else
           condition = .false.
        end if

     case ('PM')
        if ( J2>=J2Inicial ) then
           condition = .true.
        else
           condition = .false.
        end if

     case ('SAF')
        if ( J2>=J2Inicial ) then
           condition = .true.
        else
           condition = .false.
        end if

     case default
        print *, 'State inválido'

     end select

   end function

subroutine Ham_inter_state(state,J2,J3,s,m_guess,H_inter)
implicit none
      integer, dimension(maxConfig,num_sites), intent(in):: s
      integer :: i
      real(kind=db), intent(in):: J2,J3
      real(kind=db), dimension(num_sites), intent(in):: m_guess
      real(kind=db), dimension(maxConfig), intent(out):: H_inter
      character(len=3), intent(in) :: state


select case (state)
case ('AF')
   do i = 1, maxConfig
      H_inter(i) = (-2*J1 + 3*J2 + 4*J3)*(s(i,1)-s(i,2)-s(i,3)+s(i,4)-2*m_guess(1))*(m_guess(1))
   enddo
case ('SAF')
   do i = 1, maxConfig
      H_inter(i) = (-3*J2+4*J3)*(s(i,1)-s(i,2)+s(i,3)-s(i,4)-2*m_guess(1))*(m_guess(1))
   enddo
case ('SD')
   do i = 1, maxConfig
      H_inter(i) = (-2*J1+J2)*(s(i,1)+s(i,2)-s(i,3)-s(i,4)-2*m_guess(1))*(m_guess(1))
   enddo
case ('PM')
   H_inter = 0
case default
   print *, 'State inválido'
end select

   end subroutine

end module CMF
