program normal
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites) :: s
   integer:: mJ
   real(kind=db), parameter:: J3=0.0d0,T=1.7d0
   real(kind=db):: J2
   real(kind=db), dimension(maxConfig):: H_intra, H_inter, H
   real(kind=db), dimension(num_sites):: m_guess
   real(kind=db):: Z,m,step,tol,error,F
   character(len=5) :: nameFileT, nameFileJ3
   character(len=3) :: state

   state = 'SAF'
   mJ = 1
   m = 1.d0
   step = (10.d0)**(-5)
   tol = (10.d0)**(-5)
   !tolJ = (10.d0)**(-3)


   if ( state=='AF' ) then
      J2 = 0.4d0
   else
      J2 = 0.6d0
   end if
   call base(s)

   !Automatização criação de arquivos
   !Conversão real -> string
   WRITE (nameFileT, '(F5.2)') T
   WRITE (nameFileJ3, '(F5.2)') j3


   open(unit=20, file=trim(state) // "_J2-F_T(" // trim(adjustl(nameFileT)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")

   do while (J2>0.4d0)

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



      call F_helm(Z,F)

      write(20,*) J2,F

      if ( state=='AF' ) then
         J2 = J2 + step
      else 
         J2 = J2 - step
      end if

  

   enddo

   close(20)

end program normal



