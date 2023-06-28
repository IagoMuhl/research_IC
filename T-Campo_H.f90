program normal_1
    use CMF
    implicit none
    integer, dimension(maxConfig,num_sites) :: s
    real(kind=db), parameter:: J3 = 0.0d0, tol = (10.d0)**(-8)
    real(kind=db), dimension(maxConfig):: H_intra, H_inter, Ham
    real(kind=db), dimension(num_sites):: m_guess
    real(kind=db):: Z,T,m,step,error_tot,error_1,error_2,F,H_final
    real(kind=db):: J2 , H, m1, m2
    character(len=5) :: nameFileJ2, nameFileJ3
    character(len=3) :: state
    integer::k,i
    !REAL:: time_begin, time_end
 
    step = (10.d0)**(-3);       J2 = -1.0
 
    do
 
       T = 0.0
       k = 0
       print*, 'Entre com T'
       read(*,*) T
       if ( T==-1 ) stop 'Fim da rotina'
          
       print*, 'Entre com a fase (AF,SAF,SD,PM)'
       read(*,*)   state
 
    !    print*, 'Entre com T final'
    !    read(*,*)   Tfinal
 
       i = 0      
 
   !CALL CPU_TIME ( time_begin )
       !-----COMPUTA O TEMPO DE COMPILAÇÃO--------- 
 
 
       call base(s)
 

       !call HAM_INTRA(J2second,s,H_intra2)
 
       !Automatização criação de arquivos
       !Conversão real -> string
       WRITE (nameFileJ2, '(F5.2)') j2
       WRITE (nameFileJ3, '(F5.2)') j3
 
       !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
       open(unit=20, file=trim(state) // "_H_T-F-m.dat")
       !open(unit=21, file=trim(state2) // "_T-m.dat")
 
       H = 0.d0
       H_final = 10.d0
       m1 = 1.d0    ! chute magnetização!
       m2 = -1.d0
 
       do while (H<=H_final)
 
         call HAM_INTRA(J2,H,s,H_intra)
 
          call mag_vetor(state,m1,m2,m_guess)
 
 
 
          call Ham_inter_state(state,J2,J3,s,m_guess,H_inter)
 
          Ham = H_intra + H_inter 
 
 
 
          call partition(Ham,T,Z)
 
          call magnetization(Ham,Z,s,1,T,m1)
 
          call magnetization(Ham,Z,s,2,T,m2)
 
 
 
          error_1 = abs(m1 - m_guess(1))
          error_2 = abs(m2 - m_guess(2))
 
 
 
          error_tot = max(error_1,error_2)
 
          !AUTO-CONSISTÊNCIA
 
          do while (error_tot >= tol)
             !ATUALIZANDO O CHUTE
            call HAM_INTRA(J2,H,s,H_intra)
 
             call mag_vetor(state,m1,m2,m_guess)
 
             call Ham_inter_state(state,J2,J3,s,m_guess,H_inter)
 
 
             Ham = H_intra + H_inter
 
             call partition(Ham,T,Z)
 
             call magnetization(Ham,Z,s,1,T,m1)
 
             call magnetization(Ham,Z,s,2,T,m2)
 
             m = (abs(m1)+abs(m2))/2.d0
 
             ! print*, m
             ! read(*,*)
 
             error_1 = abs(m1 - m_guess(1))
             error_2 = abs(m2 - m_guess(2))
    
             error_tot = max(error_1,error_2)
 
             ! print *, T, m1, m_guess(1), error_tot
 
             m_guess(1) = m1
             m_guess(2) = m2
 
          enddo
 
          call F_helm(T,Z,F)
 
          !if (m<=tol) then
             !i = i + 1
             !if (i==1) then
            ! print *, T,m
             !endif
          !endif
 
          
 
          write(20,*) H,F,m,m1,m2

          print*, H, T, m
 
!           if (k==0) then
!              if (abs(m)<=10.d0**(-4)) then
!                 print*, '------------'
!                 write(*,18) T, H
!  18             format ((F8.5))
!                 k = 1
!              end if
!           end if
 
          H = H + step
 
       enddo
 
       ! write(*,*) '---------------'
       ! print*, 'J2:',J2
       ! print*, 'T', T
       ! print*, 'Fase', state
       ! write(*,*) '------FIM------'
       !-----COMPUTA O TEMPO DE COMPILAÇÃO--------- 
       !CALL CPU_TIME ( time_end )
       !WRITE (*,*) 'Time of operation was ', time_end - time_begin, ' seconds'
    
       close(20)
 
    end do
 
 
 
    ! call print_matrix(maxConfig,num_sites,s)
    !print *, T,Z
    !call print_matrixH(maxConfig,1,Ham)
 
 end program normal_1
 