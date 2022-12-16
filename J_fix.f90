program normal
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites) :: s
   real(kind=db), parameter:: J3 = -0.1d0, tol = (10.d0)**(-8)
   !real(kind=db), parameter:: J2 = 0.51d0, J2second = 0.51d0
   real(kind=db), dimension(maxConfig):: H_intra, H_inter, H
   ! H_intra2, H_inter2, H2
   real(kind=db), dimension(num_sites):: m_guess
   !,m_guess2
   real(kind=db):: Z,T,m,step,error,F
   !m2,error2,Z2,F2,Dif,stepJ2
   real(kind=db):: J2 !J2second
   character(len=5) :: nameFileJ2, nameFileJ3
   character(len=3) :: state !state2
   integer:: first_m0 ! first_m02

   step = (10.d0)**(-3)
   !state2 = 'SAF'
   state = 'AF'
   J2 = 0.53d0
   !J2second = 1.d0

   call base(s)

   call HAM_INTRA(J2,s,H_intra)
   !call HAM_INTRA(J2second,s,H_intra2)

   !Automatização criação de arquivos
   !Conversão real -> string
   WRITE (nameFileJ2, '(F5.2)') j2
   WRITE (nameFileJ3, '(F5.2)') j3

   !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
   open(unit=20, file=trim(state) // "_T-m.dat")
   !open(unit=21, file=trim(state2) // "_T-m.dat")

   !do while (J2<1.d0)

      !stepJ2 = (10.d0)**(-1)  !STEP DE J2
      T = 0.1d0


      m = 1.d0    ! chute magnetização!
     ! m2 = 1.d0


      first_m0 = 1    !PRINT PARA M<=TOL
     !first_m02 = 1



   do while (T<=4.0d0)
      !print*, m2, T
      !read(*,*)
      !call print_m(state,J2,m,tol,T,first_m0) !SAF
      !call print_m(state2,J2second,m2,tol,T,first_m02) !AF


      call mag_vetor(state,m,m_guess)
      !call mag_vetor(state2,m2,m_guess2)


      call HAM_INTER(J2,J3,s,m_guess,H_inter)
      !call HAM_INTER(J2second,J3,s,m_guess2,H_inter2)

      H = H_intra + H_inter
     ! H2 = H_intra2 + H_inter2

      call partition(H,T,Z)
      !call partition(H2,T,Z2)

      call magnetization(H,Z,s,T,m)
      !call magnetization(H2,Z2,s,T,m2)

      error = abs(m - m_guess(1))
      !error2 = abs(m2 - m_guess2(1))
         !AUTO-CONSISTÊNCIA

      do while (error >= tol)
         !ATUALIZANDO O CHUTE
        
         call mag_vetor(state,m,m_guess)

         call HAM_INTER(J2,J3,s,m_guess,H_inter)

         H = H_intra + H_inter

         call partition(H,T,Z)

         call magnetization(H,Z,s,T,m)

         error = abs(m - m_guess(1))

      enddo
      !do while (error2 >= tol)
         !ATUALIZANDO O CHUTE
        
        ! call mag_vetor(state2,m2,m_guess2)

        ! call HAM_INTER(J2second,J3,s,m_guess2,H_inter2)

         !H2 = H_intra2 + H_inter2

        ! call partition(H2,T,Z2)

        ! call magnetization(H2,Z,s,T,m2)

        ! error2 = abs(m2 - m_guess2(1))

      !enddo


      
      




      call F_helm(T,Z,F)
      !call F_helm(T,Z2,F2)

      !---------------------------  PONTO CRÍTICO
     ! Dif = abs(F) - abs(F2)  

      !if(Dif<=tol) then
       !  print*, T,J2,Dif
        ! stop
      !end if
      !---------------------------
      write(20,*) T,m
     ! write(21,*) T,m2
      !print*, m_guess
      T = T + step

   enddo

   !read(*,*)


   !J2 = J2 + stepJ2
   !J2second = J2second - stepJ2

   !enddo


   close(20)

   

   ! call print_matrix(maxConfig,num_sites,s)
   !print *, T,Z
   !call print_matrixH(maxConfig,1,H)

end program normal
