program normal
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites) :: s
   real(kind=db), parameter:: J3 = -0.2d0, tol = (10.d0)**(-8)
   real(kind=db), dimension(maxConfig):: H_intra, H_inter, H
   real(kind=db), dimension(num_sites):: m_guess
   real(kind=db):: Z,T,m,step,error,F,Tfinal
   real(kind=db):: J2 
   character(len=5) :: nameFileJ2, nameFileJ3
   character(len=3) :: state
   integer::i
   !REAL:: time_begin, time_end

   step = (10.d0)**(-5)

   
   
   do

      J2 = 0.0

      print*, 'Entre com J2'
      read(*,*) J2
      if ( J2==-1 ) stop 'Fim da rotina'
         
      print*, 'Entre com a fase (AF,SAF,SD,PM)'
      read(*,*)   state

      print*, 'Entre com T final'
      read(*,*)   Tfinal

      i = 0      

  !CALL CPU_TIME ( time_begin )
      !-----COMPUTA O TEMPO DE COMPILAÇÃO--------- 


      call base(s)

      call HAM_INTRA(J2,s,H_intra)
      !call HAM_INTRA(J2second,s,H_intra2)

      !Automatização criação de arquivos
      !Conversão real -> string
      WRITE (nameFileJ2, '(F5.2)') j2
      WRITE (nameFileJ3, '(F5.2)') j3

      !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
      open(unit=20, file=trim(state) // "_T-F-m.dat")
      !open(unit=21, file=trim(state2) // "_T-m.dat")

      T = 0.1d0
      m = 1.d0    ! chute magnetização!
      
      do while (T<=Tfinal)
         !print *, T

         call mag_vetor(state,m,m_guess)


         call Ham_inter_state(state,J2,J3,s,m_guess,H_inter)

         H = H_intra + H_inter

         call partition(H,T,Z)

         call magnetization(H,Z,s,T,m)

         error = abs(m - m_guess(1))
         !AUTO-CONSISTÊNCIA

         do while (error >= tol)
            !ATUALIZANDO O CHUTE

            call mag_vetor(state,m,m_guess)

            call Ham_inter_state(state,J2,J3,s,m_guess,H_inter)


            H = H_intra + H_inter

            call partition(H,T,Z)

            call magnetization(H,Z,s,T,m)

            error = abs(m - m_guess(1))

         enddo

         call F_helm(T,Z,F)

         !if (m<=tol) then
            !i = i + 1
            !if (i==1) then
           ! print *, T,m
            !endif
         !endif

         

         write(20,*) T,F,m

         T = T + step

      enddo

      write(*,*) '---------------'
      print*, 'J2:',J2
      print*, 'T', T
      print*, 'Fase', state
      write(*,*) '------FIM------'
      !-----COMPUTA O TEMPO DE COMPILAÇÃO--------- 
      !CALL CPU_TIME ( time_end )
      !WRITE (*,*) 'Time of operation was ', time_end - time_begin, ' seconds'
   
      close(20)

   end do



   ! call print_matrix(maxConfig,num_sites,s)
   !print *, T,Z
   !call print_matrixH(maxConfig,1,H)

end program normal
