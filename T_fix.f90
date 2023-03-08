program normal
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites) :: s
   real(kind=db), parameter:: T=1.d0, J3=0.0d0
   real(kind=db):: J2, J2Final,J2Inicial
   real(kind=db), dimension(maxConfig):: H_intra, H_inter, H
   real(kind=db), dimension(num_sites):: m_guess
   real(kind=db):: Z,m,step,tol,error,F
   character(len=5) :: nameFileT, nameFileJ3
   character(len=3) :: state, stateOne, stateTwo
   integer:: i

   !  não faça a cagada de mudar o STEP. Ass.: Matheus


   step = (10.d0)**(-5)
   tol = (10.d0)**(-8)

   print*, 'Entre com a fase crescente(AF,SAF,SD)'
      read(*,*)   stateOne

      print*, 'Entre com J2 inicial'
      read(*,*)   J2Inicial

      print*, 'Entre com J2 final'
      read(*,*)   J2Final

   print*, 'Entre com a fase decrescente(AF,SAF,SD)'
      read(*,*)   stateTwo


   call base(s)

   do i = 1, 2
      if ( i==1 ) then
         state = stateOne
         J2 = J2Inicial
      else
         state = stateTwo
         J2 = J2Final
         step = step*(-1)
      end if


      ! devemos 'zerar' a mag a cada loop das fases
      m = 1.d0




      !Automatização criação de arquivos
      !Conversão real -> string
      WRITE (nameFileT, '(F5.2)') T
      WRITE (nameFileJ3, '(F5.2)') j3

      ! AF --> (J2<0.6d0)
      ! SAF --> (J2>0.4d0)


      !open(unit=20, file=trim(state) // "-PM_J2-F_T(" // trim(adjustl(nameFileT)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
      open(unit=20, file=trim(state) // "_J2-F-m.dat")




      do while (condition(stateOne,stateTwo,J2,J2Final,J2Inicial))


         m_guess = [m,-m,m,-m]
         call mag_vetor(state,m,m_guess)

         call HAM_INTRA(J2,s,H_intra)

         call  Ham_inter_state(state,J2,J3,s,m_guess,H_inter)

         H = H_intra + H_inter

         call partition(H,T,Z)

         call magnetization(H,Z,s,T,m)
         m = abs(m)
         error = abs(m - m_guess(1))

         do while (error >= tol)
            !ATUALIZANDO O CHUTE

            call mag_vetor(state,m,m_guess)

            call  Ham_inter_state(state,J2,J3,s,m_guess,H_inter)

            H = H_intra + H_inter

            call partition(H,T,Z)

            call magnetization(H,Z,s,T,m)

            error = abs(m - m_guess(1))

         enddo


         call F_helm(T,Z,F)

         write(20,*) J2,F,m
         !print*, F,m
         !read(*,*)
         !write(27,*) J2,m



         J2 = J2 + step


      enddo

      close(20)
   end do
end program normal



