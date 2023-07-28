program square_T
   use CMF
   implicit none
   real*8, parameter:: J2 = -1.d0
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total
   real*8, dimension(num_sites):: m
   real*8:: T_inicial,T_final,H,m_fe,m_af,step,Z,m_order,tol,F,erro1,erro2
   character(len=3):: state
   integer:: j, up, down, cd

   tol = 10.d0**(-8)

   call base(s)
!--------------------------------------------------------------

   cd = -3
   step = 10.d0**(cd)


   write(*,*) 'Entre com T_inicial:'
   read*, T_inicial

   write(*,*) 'Entre com T_final:'
   read*, T_final

   ! -
   if ( T_inicial>T_final ) then
      step = -10.d0**(cd)
   else
      step = 10.d0**(cd)
   end if
!-


   do

      j = 0


      write(*,*) 'Entre com H:'
      read*, H

      write(*,*) 'Entre com a fase:'
      read*, state

      m_fe = 1.d0; m_af = -1.d0

      ! - - - - - - - - - - - - - - - - - - - - - - -

      call HAM_INTRA(J2,H,s,H_intra)

      call mag_vetor(state,m_fe,m_af,up,down,m,m_order)

      call Ham_inter_state(J2,s,m,H_inter)

      H_total = H_intra + H_inter

      call partition(H_total,T_inicial,Z)

      call magnetization(H_total,Z,s,up,T_inicial,m_fe)

      call magnetization(H_total,Z,s,down,T_inicial,m_af)

      call mag_vetor(state,m_fe,m_af,up,down,m,m_order)

      ! - - - - - - - - - - - - - - - - - - - - - - -



      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_H_T-F-m.dat")


      do while (T_inicial/=T_final)


         erro1 = 1.d0; erro2 = 1.d0

         do while(max(erro1,erro2) >= tol)


            call Ham_inter_state(J2,s,m,H_inter)

            H_total = H_intra + H_inter

            call partition(H_total,T_inicial,Z)

            call magnetization(H_total,Z,s,up,T_inicial,m_fe)

            call magnetization(H_total,Z,s,down,T_inicial,m_af)


            erro1 = abs(m_fe - m(up))
            erro2 = abs(m_af - m(down))

            if(state=='PM') then
               erro2 = abs(m_fe - m(down))
            endif

            !  print*, m_fe, m_af, m_order
            !  read(*,*)

            call mag_vetor(state,m_fe,m_af,up,down,m,m_order)


         end do

         call F_helm(T_inicial,Z,F)

         !print*, T_inicial, m_order


         write(20,*) T_inicial, F, m_order
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


         write(*,*) 'Entre com T_inicial:'
         read*, T_inicial

         write(*,*) 'Entre com T_final:'
         read*, T_final

         ! -
         if ( T_inicial>T_final ) then
            step = -10.d0**(cd)
         else
            step = 10.d0**(cd)
         end if
         !-------------------------------



   enddo

end program


