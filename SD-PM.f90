program normal
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites) :: s
   real(kind=db), parameter:: J3 = 0.2d0, tol = (10.d0)**(-8)
   real(kind=db), dimension(maxConfig):: H_intra, H_inter, H
   real(kind=db), dimension(num_sites):: m_guess
   real(kind=db):: Z,T,m,step,error,F
   real(kind=db):: J2
   character(len=5) :: nameFileJ2, nameFileJ3
   character(len=3) :: state, stateOne, stateTwo
   integer:: i

   step = (10.d0)**(-5)

   stateOne = 'SD'
   stateTwo = 'PM'

   do

      write(*,*) 'Entre com J2'
      read(*,*) J2
      if ( J2==-1 ) stop 'Fim da rotina'

      do i = 1, 2
         if ( i==1 ) then
            state = stateOne
         else
            state = stateTwo
         end if

         call base(s)

         call HAM_INTRA(J2,s,H_intra)

         !Automatização criação de arquivos
         !Conversão real -> string
         WRITE (nameFileJ2, '(F5.2)') j2
         WRITE (nameFileJ3, '(F5.2)') j3

         !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
         open(unit=20, file=trim(state) // "_T-F-m.dat")

         T = 0.1d0
         m = 1.d0    ! chute magnetização!

         do while (T<=4.0d0)

            call mag_vetor(state,m,m_guess)

            if ( state=='SD' ) then
               call HAM_INTER_SD(J2,J3,s,m_guess,H_inter)
            else
               call HAM_INTER(J2,J3,s,m_guess,H_inter)
            end if

            H = H_intra + H_inter

            call partition(H,T,Z)

            call magnetization(H,Z,s,T,m)

            error = abs(m - m_guess(1))

            !AUTO-CONSISTÊNCIA

            do while (error >= tol)
               !ATUALIZANDO O CHUTE

               call mag_vetor(state,m,m_guess)

               if ( state=='SD' ) then
                  call HAM_INTER_SD(J2,J3,s,m_guess,H_inter)
               else
                  call HAM_INTER(J2,J3,s,m_guess,H_inter)
               end if

               H = H_intra + H_inter

               call partition(H,T,Z)

               call magnetization(H,Z,s,T,m)

               error = abs(m - m_guess(1))

            enddo

            call F_helm(T,Z,F)

            write(20,*) T,F,m

            T = T + step

         enddo

         close(20)

      end do
   enddo
end program normal
