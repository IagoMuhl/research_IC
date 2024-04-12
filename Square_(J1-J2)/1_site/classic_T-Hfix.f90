program normal
    use CMF
    implicit none
    integer, dimension(maxConfig,num_sites) :: s
    real(kind=db):: J2
    real(kind=db), dimension(maxConfig):: H_intra, H_inter, H_long, H_inter_2
    real(kind=db), dimension(num_sites*2):: m_guess
    real(kind=db):: Z,m_fe,m_af,step,tol,erro,F,T_inicial,T_final,cd,m1,m2,erro1,erro2,Z2,m_order,H
    character(len=5) :: nameFileT, nameFileJ3
    character(len=3) :: state
    integer:: i
 
    tol = (10.d0)**(-8); J2 = 0.d0

    call base(s)
 
       ! devemos 'zerar' a mag a cada loop das fases
       m_fe = 1.d0
       m_af = -1.d0
 


      do
         i = 0

         H_intra = 0.d0

         print*, 'Entre com H e Step(-5,-3):'
         read(*,*)  H, cd

         step = (10.d0)**(cd)

         ! print*, 'Entre com T_inicial'
         ! read(*,*) T_inicial
   
         ! print*, 'Entre com T_final'
         ! read(*,*) T_final

         T_inicial = 0.5d0
         T_final = 8.d0

         print*, 'Entre com a fase (AF,SAF,SD,PM)'
         read(*,*)   state
 
       !open(unit=20, file=trim(state) // "-PM_J2-F_T(" // trim(adjustl(nameFileT)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
       open(unit=20, file=trim(state) // "_Hfix_T-F-m.dat")
 
 
          call mag_vetor(state,m_fe,m_af,m_guess)
 
         do while (T_inicial/=T_final) 

         erro1 = 1.d0; erro2 = 1.d0; erro = 1.d0

         do while (erro >= tol)

            
 
          call  Ham_inter_state(J2,s,H,m_guess,H_inter)
          call  Ham_inter_state(J2,s,H,[m_guess(2),m_guess(1)],H_inter_2)

          !H = H_intra + H_inter
 
          call partition(H_inter,T_inicial,Z)
          call partition(H_inter_2,T_inicial,Z2)

         !  print*, m_guess, erro
         !    read(*,*)
 
         !    print*, H_inter_2
         !    read(*,*)

          call magnetization(H_inter,Z,s,1,T_inicial,m1)
          call magnetization(H_inter_2,Z2,s,1,T_inicial,m2)

          erro1 = abs(m_guess(1) - m1)
          erro2 = abs(m_guess(2) - m2)

          erro = max(abs(erro1),abs(erro2))
 
          call mag_vetor(state,m1,m2,m_guess)

          m_order = abs(m1 - m2)/2.d0

          enddo
 
          call F_helm(T_inicial,Z,F)
 
          write(20,*) T_inicial, F, m_order, m1, m2
          !print*, T_inicial, m_order, m1, m2
          !read(*,*)
          !write(27,*) J2,m
 
          if (i==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '\/------------\/'
               write(*,18) H, T_inicial
               print*, '/\------------/\'
18             format ((F8.5))
               i = 1
            end if
         end if
 
          T_inicial = T_inicial + step
 
          if ((max(T_inicial,T_final))-(min(T_inicial,T_final))<=abs(step)) then
            T_inicial = T_final
         endif

         enddo

         close(20)

      enddo
 
       

 end program normal
 
 
 
 