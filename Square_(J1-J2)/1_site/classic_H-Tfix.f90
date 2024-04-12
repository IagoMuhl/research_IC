program normal
    use CMF
    implicit none
    integer, dimension(maxConfig,num_sites) :: s
    real(kind=db):: J2
    real(kind=db), dimension(maxConfig):: H_intra, H_inter, H_long, H_inter_2
    real(kind=db), dimension(num_sites*2):: m_guess
    real(kind=db):: Z,m_fe,m_af,step,tol,erro,F,H_inicial,H_final,cd,m1,m2,erro1,erro2,Z2,m_order,T
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

         print*, 'Entre com T e Step(-5,-3):'
         read(*,*)  T, cd

         step = (10.d0)**(cd)

         ! print*, 'Entre com H_inicial'
         ! read(*,*) H_inicial
   
         ! print*, 'Entre com H_final'
         ! read(*,*) H_final

         H_inicial = 0.5d0
         H_final = 8.d0

         print*, 'Entre com a fase (AF,SAF,SD,PM)'
         read(*,*)   state
 
       !open(unit=20, file=trim(state) // "-PM_J2-F_T(" // trim(adjustl(nameFileT)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
       open(unit=20, file=trim(state) // "_Tfix_H-F-m.dat")
 
 
          call mag_vetor(state,m_fe,m_af,m_guess)
 
         do while (H_inicial/=H_final) 

         erro1 = 1.d0; erro2 = 1.d0; erro = 1.d0

         do while (erro >= tol)

          call  Ham_inter_state(J2,s,H_inicial,m_guess,H_inter)
          call  Ham_inter_state(J2,s,H_inicial,[m_guess(2),m_guess(1)],H_inter_2)

          !H = H_intra + H_inter
 
          call partition(H_inter,T,Z)
          call partition(H_inter_2,T,Z2)

         !  print*, m_guess, erro
         !    read(*,*)
 
         !    print*, H_inter_2
         !    read(*,*)

          call magnetization(H_inter,Z,s,1,T,m1)
          call magnetization(H_inter_2,Z2,s,1,T,m2)

          erro1 = abs(m_guess(1) - m1)
          erro2 = abs(m_guess(2) - m2)

          erro = max(abs(erro1),abs(erro2))
 
          call mag_vetor(state,m1,m2,m_guess)

          m_order = abs(m1 - m2)/2.d0

          enddo
 
          call F_helm(T,Z,F)
 
          write(20,*) H_inicial, F, m_order, m1, m2
          !print*, T_inicial, m_order, m1, m2
          !read(*,*)
          !write(27,*) J2,m
 
          if (i==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '\/------------\/'
               write(*,18) H_inicial, T
               print*, '/\------------/\'
18             format ((F8.5))
               i = 1
            end if
         end if
 
          H_inicial = H_inicial + step
 
          if ((max(H_inicial,H_final))-(min(H_inicial,H_final))<=abs(step)) then
            H_inicial = H_final
         endif

         enddo

         close(20)

      enddo
 
       

 end program normal
 
 
 
 