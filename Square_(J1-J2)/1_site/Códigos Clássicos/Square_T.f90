program square_T
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_inter_1, H_inter_2, H_total_1, H_total_2 ,s_z
   real*8, dimension(minConfig):: m
   real*8:: J2, erro, m1, m2, Z1, Z2, F2, U1, U2, U
   real*8:: T_inicial,T_final,H,m_fe,m_af,step,m_order,tol,F,erro1,erro2
   real*8:: St, So, Cv
   character(len=3):: state
   integer:: j, up, down, cd,i, conts

   up = 1; down = 2
   tol = 10.d0**(-8); J2 = -0.4d0; conts = 1

   call base(s)
!--------------------------------------------------------------
   do i = 1, maxConfig
      s_z(i) = s(i,1)
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


      ! - - - - - - - - - - - - - - - - - - - - - - -


      call mag_vetor(state,m_fe,m_af,m,m_order)

      !call HAM_INTRA(J2,s,H_intra)
      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_H_T-F-m.dat")


      do while (T_inicial/=T_final)

         erro1 = 1.d0; erro2 = 1.d0; erro = 1.d0

         do while(erro >= tol)

            call Ham_inter_state(J2,s,[m(1),m(2)],H_inter_1)

            call Ham_inter_state(J2,s,[m(2),m(1)],H_inter_2)

            
            H_total_1 =  H_inter_1 - H*s_z

            H_total_2 =  H_inter_2 - H*s_z


            call partition(H_total_1,T_inicial,Z1)

            call partition(H_total_2,T_inicial,Z2)


            call magnetization(H_total_1,Z1,s,T_inicial,m1)

            call magnetization(H_total_2,Z2,s,T_inicial,m2)

            !print*, T_inicial, m1, m2

            erro1 = abs(m(1) - m1)
            erro2 = abs(m(2) - m2)

            erro = max(abs(erro1),abs(erro2))

            !  print*, m_fe, m_af, m_order
            !  read(*,*)

            call mag_vetor(state,m1,m2,m,m_order)


         end do

         call F_helm(T_inicial,Z1,F)
         call F_helm(T_inicial,Z2,F2)

         call Inner_energy(H_total_1,Z1,T_inicial,U1)
         call Inner_energy(H_total_2,Z2,T_inicial,U2)

         F = (F + F2)/2.d0
         U = (U1 + U2)/2.d0

         St = (U - F)/(num_sites*T_inicial)

         !Cv = T_inicial*((St - So)/(step))

      
         If(conts == 1) then
            So = St 
            Cv = 0
         else
            Cv = T_inicial*((St - So)/(step))
            So = St
         end if

         conts = 0
         !print*, T_inicial, m_order


         write(20,*) T_inicial, F, m_order, St, Cv
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '\/---------\/'
               write(*,18) H, T_inicial
               print*, '/\---------/\'
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


