program normal
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites) :: s
   integer:: mJ
   real(kind=db), parameter:: J2=0.61d0,J3=0.0d0
   real(kind=db), dimension(maxConfig):: H_intra, H_inter, H
   real(kind=db), dimension(num_sites):: m_guess
   real(kind=db):: Z,T,m,step,tol,error,tolJ,F
   character(len=5) :: nameFileJ2, nameFileJ3
   character(len=3) :: state

   state = 'PM'
   mJ = 1
   m = 0.d0
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
      call mag_vetor(state,m,m_guess)

      call HAM_INTRA(J2,s,H_intra)

      call HAM_INTER(J2,J3,s,m_guess,H_inter)

      H = H_intra + H_inter

      call partition(H,T,Z)

      call magnetization(H,Z,s,T,m)
      m = abs(m)
      error = abs(m - m_guess(1))

      do while (error >= tol)
         !ATUALIZANDO O CHUTE
        
         call mag_vetor(state,m,m_guess)

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

      call F_helm(Z,F)

      write(20,*) T,F

      !print*, m_guess
      T = T + step

   enddo

   close(20)

   ! call print_matrix(maxConfig,num_sites,s)
   !print *, T,Z
   !call print_matrixH(maxConfig,1,H)

end program normal
