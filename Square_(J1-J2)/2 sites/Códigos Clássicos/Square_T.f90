program square_T
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_x
   real*8, dimension(num_sites):: m
   real*8:: J2, erro, m1, m2
   real*8:: T_inicial,T_final,H,m_fe,m_af,step,Z,m_order,tol,F,erro1,erro2
   character(len=3):: state
   integer:: j, up, down, cd,i

   up = 1; down = 2
   tol = 10.d0**(-8); J2 = -0.465d0

   call base(s)
!--------------------------------------------------------------
   do i = 1, maxConfig
      s_x(i) = s(i,1) + s(i,2)
   enddo

   do

      j = 0


      write(*,*) 'Entre com H, step:'
      read*, H, cd

      write(*,*) 'Entre com T_inicial:'
      read*, T_inicial

      write(*,*) 'Entre com T_final:'
      read*, T_final




      write(*,*) 'Entre com a fase:'
      read*, state


      ! -
      if ( T_inicial>T_final ) then
         step = -10.d0**(cd)
      else
         step = 10.d0**(cd)
      end if
!-

      if(state/='2AF') then
         m_fe = 1.0d0;
         m_af = -1.0d0;
      else
         m_fe = 0.84719110987579493
         m_af = 0.65197042353076307
         !   m_fe = 0.99099828895786968!1.0d0;
         !   m_af =  0.44969820018918100 !-1.0d0;
      endif

      ! - - - - - - - - - - - - - - - - - - - - - - -

      ! call HAM_INTRA(J2,H,s,H_intra)

      ! call mag_vetor(state,m_fe,m_af,m,m_order)

      ! call Ham_inter_state(J2,s,m,H_inter)

      ! H_total = H_intra + H_inter

      ! call partition(H_total,T_inicial,Z)

      ! call magnetization(H_total,Z,s,up,T_inicial,m_fe)

      ! call magnetization(H_total,Z,s,down,T_inicial,m_af)

      ! call mag_vetor(state,m_fe,m_af,m,m_order)

      ! - - - - - - - - - - - - - - - - - - - - - - -


      call mag_vetor(state,m_fe,m_af,m,m_order)

      call HAM_INTRA(J2,s,H_intra)
      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_H_T-F-m.dat")


      do while (T_inicial/=T_final)

         erro1 = 1.d0; erro2 = 1.d0; erro = 1.d0

         do while(erro >= tol)

            call Ham_inter_state(J2,s,m,H_inter)

            H_total = H_intra + H_inter - H*s_x

            call partition(H_total,T_inicial,Z)

            call magnetization(H_total,Z,s,up,T_inicial,m1)

            call magnetization(H_total,Z,s,down,T_inicial,m2)


            erro1 = abs(m(1) - m1)
            erro2 = abs(m(2) - m2)

            erro = max(abs(erro1),abs(erro2))

            !  print*, m_fe, m_af, m_order
            !  read(*,*)

            call mag_vetor(state,m1,m2,m,m_order)


         end do

         call F_helm(T_inicial,Z,F)

         !print*, T_inicial, m_order


         write(20,*) T_inicial, F, m_order, m1, m2
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '------------'
               write(*,18) H, T_inicial
18             format ((F8.5))
               j = 1
            end if
         end if



         T_inicial = T_inicial + step

         if ((max(T_inicial,T_final))-(min(T_inicial,T_final))<=abs(step)) then
            T_inicial = T_final
         endif

         ! if (max(T_final,T_inicial)>=12.d0) then
         !    exit
         ! endif

      end do

      close(20)


      !----------------- NOVAS ENTRADAS  -----------------


      print*, '------',State,'------'
      write(*,19) T_inicial, H
19    format ((F8.5))



   enddo

end program


