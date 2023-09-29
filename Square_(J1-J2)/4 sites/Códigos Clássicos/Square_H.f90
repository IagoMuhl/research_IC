program square_H
   use CMF
   implicit none
   real*8, parameter:: J2 = -1.d0
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total
   real*8, dimension(num_sites):: m
   real*8:: H_inicial,H_final,T,m_fe,m_af,step,Z,m_order,tol,F,erro1,erro2
   character(len=3):: state
   integer:: j, up, down, cd

   tol = 10.d0**(-8)

   call base(s)
!--------------------------------------------------------------

   cd = -5
   step = 10.d0**(cd)


   write(*,*) 'Entre com H_inicial:'
   read*, H_inicial

   write(*,*) 'Entre com H_final:'
   read*, H_final

   ! -
   if ( H_inicial>H_final ) then
      step = -10.d0**(cd)
   else
      step = 10.d0**(cd)
   end if
!-


   do

      j = 0


      write(*,*) 'Entre com T:'
      read*, T

      write(*,*) 'Entre com a fase:'
      read*, state

      m_fe = 1.d0; m_af = -1.d0

      ! - - - - - - - - - - - - - - - - - - - - - - -

      call HAM_INTRA(J2,H_inicial,s,H_intra)

      call mag_vetor(state,m_fe,m_af,up,down,m,m_order)

      call Ham_inter_state(J2,s,m,H_inter)

      H_total = H_intra + H_inter

      call partition(H_total,T,Z)

      call magnetization(H_total,Z,s,up,T,m_fe)

      call magnetization(H_total,Z,s,down,T,m_af)

      call mag_vetor(state,m_fe,m_af,up,down,m,m_order)

      ! - - - - - - - - - - - - - - - - - - - - - - -



      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_T_H-F-m.dat")


      do while (H_inicial/=H_final)

         call HAM_INTRA(J2,H_inicial,s,H_intra)

         ! do j=1, maxConfig
         !    write(*,*) H_inter(j)
         ! end do

         erro1 = 1.d0; erro2 = 1.d0

         do while(max(erro1,erro2) >= tol)


            call Ham_inter_state(J2,s,m,H_inter)

            H_total = H_intra + H_inter

            call partition(H_total,T,Z)

            call magnetization(H_total,Z,s,up,T,m_fe)

            call magnetization(H_total,Z,s,down,T,m_af)


            erro1 = abs(m_fe - m(up))
            erro2 = abs(m_af - m(down))

            if(state=='PM') then
               erro2 = abs(m_fe - m(down))
            endif

            !  print*, m_fe, m_af, m_order
            !  read(*,*)

            call mag_vetor(state,m_fe,m_af,up,down,m,m_order)


         end do

         call F_helm(T,Z,F)

         !print*, H_inicial, m_order


         write(20,*) H_inicial, F, m_order
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '------------'
               write(*,18) H_inicial, T
18             format ((F8.5))
               j = 1
            end if
         end if



         H_inicial = H_inicial + step

         if ((max(H_inicial,H_final))-(min(H_inicial,H_final))<=abs(step)) then
            H_inicial = H_final
         endif

         ! if (max(H_final,H_inicial)>=12.d0) then
         !    exit
         ! endif

      end do

      close(20)

      print*, '------',State,'------'
      write(*,19) H_inicial, T
19    format ((F8.5))

      !----------------- NOVAS ENTRADAS  -----------------


      write(*,*) '-----Entre com H_inicial-----:'
      read*, H_inicial

      write(*,*) 'Entre com H_final:'
      read*, H_final

      ! -
      if ( H_inicial>H_final ) then
         step = -10.d0**(cd)
      else
         step = 10.d0**(cd)
      end if
      !-------------------------------



   enddo

end program


