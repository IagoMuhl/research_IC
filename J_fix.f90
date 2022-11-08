program normal
   implicit none
   integer, parameter:: db=8
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites) :: s
   integer:: mJ
   real(kind=db), parameter:: J1=1.d0,J2=0.5d0,J3=0.0d0
   real(kind=db), dimension(maxConfig):: H_intra, H_inter, H
   real(kind=db), dimension(num_sites):: m_guess
   real(kind=db):: Z,T,m,step,tol,error,tolJ
   character(len=5) :: nameFileJ2, nameFileJ3
   character(len=3) :: state

   state = 'SAF'
   mJ = 1
   m = 1.d0
   T = 1.d0
   step = (10.d0)**(-3)
   tol = (10.d0)**(-5)
   tolJ = (10.d0)**(-3)


   call base(s)

   !Automatização criação de arquivos
   !Conversão real -> string
   WRITE (nameFileJ2, '(F5.2)') j2
   WRITE (nameFileJ3, '(F5.2)') j3


   open(unit=20, file=trim(state) // "_T-m-mJ_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")

   do while (T<=4.0d0)

      m_guess = [m,-m,m,-m]
      call mag_vetor(num_sites,state,m,m_guess)

      call HAM_INTRA(J2,s,H_intra)

      call HAM_INTER(J2,J3,s,m_guess,H_inter)

      H = H_intra + H_inter

      call partition(H,T,Z)

      call magnetization(H,Z,s,T,m)
      m = abs(m)
      error = abs(m - m_guess(1))

      do while (error >= tol)
         !ATUALIZANDO O CHUTE
        
         call mag_vetor(num_sites,state,m,m_guess)

         call HAM_INTER(J2,J3,s,m_guess,H_inter)

         H = H_intra + H_inter

         call partition(H,T,Z)

         call magnetization(H,Z,s,T,m)

         error = abs(m - m_guess(1))

         if ( error <= tol ) exit

      enddo

      if ( m<=tolJ ) then
         mJ = m 

      end if


      write(20,*) T,m,mJ

      !print*, m_guess
      T = T + step

   enddo

   close(20)

   ! call print_matrix(maxConfig,num_sites,s)
   !print *, T,Z
   !call print_matrixH(maxConfig,1,H)

end program normal

subroutine base(s)

   implicit none
   integer, parameter:: minConfig=2
   integer, parameter:: maxConfig= minConfig**4
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
   integer, parameter:: db=8
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites), intent(in):: s
   integer :: i
   real(kind=db), intent(in):: J2
   real(kind=db), parameter:: J1=1.d0
   real(kind=db), dimension(maxConfig), intent(out):: H_intra

   !---------------------HAMILTONIANO INTRA-----------------------------
   do i = 1, maxConfig
      H_intra(i) = J1*(s(i,1)*s(i,2) + s(i,1)*s(i,3) + s(i,4)*s(i,3) + s(i,2)*s(i,4)) &
      & + J2*(s(i,1)*s(i,4) + s(i,3)*s(i,2))

   end do
   !---------------------HAMILTONIANO INTRA-----------------------------


end subroutine

subroutine HAM_INTER(J2,J3,s,m_guess,H_inter)
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: db=8
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites), intent(in):: s
   integer :: i
   real(kind=db), intent(in):: J2,J3
   real(kind=db), parameter:: J1=1.d0
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
   integer, parameter:: db=8
   integer, intent(in):: rowH,columnH
   real(kind=db), dimension(rowH,columnH), intent(in):: matrixH
   integer:: i,j
   do i = 1, rowH
      write(*,*) (matrixH(i,j),j = 1, columnH)
   end do

end subroutine

subroutine partition(H,T,Z)
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: db=8
   integer, parameter:: maxConfig= minConfig**num_sites
   real(kind=db), intent(out):: Z
   real(kind=db), dimension(maxConfig), intent(in):: H
   real(kind=db), intent(in):: T
   real(kind=db):: b
   b = 1.d0/T
   Z = 0.d0
   m = 0.d0
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
   integer, parameter:: minConfig= 2
   integer, parameter:: num_sites= 4
   integer, parameter:: db=8
   integer, parameter:: maxConfig= minConfig**num_sites
   integer, dimension(maxConfig,num_sites), intent(in) :: s
   real(kind=db):: b
   real(kind=db), dimension(maxConfig), intent(in):: H
   real(kind=db), intent(out):: m
   real(kind=db), intent(in):: T, Z

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
subroutine mag_vetor(num_sites,state,sigma,m)
   implicit none
   integer, parameter:: db=8
   character(len=*), intent(in):: state
   integer, intent(in):: num_sites
   real(kind=db), intent(in):: sigma
   real(kind=db), dimension(num_sites), intent(out):: m

   if ( state=='AF' ) then
      m = [sigma, -sigma, -sigma, sigma]
   else if (state=='SAF') then
      m = [sigma, -sigma, sigma, -sigma]
   end if

end subroutine