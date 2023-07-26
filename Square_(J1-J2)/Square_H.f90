program square_H
   use CMF
   implicit none
   real*8, parameter:: J2 = -1.d0
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total
   real*8, dimension(num_sites):: m
   real*8:: H,T,m_fe,m_af,step,Z,m_order,tol,F,erro1,erro2
   character(len=3):: state
   integer:: j

   step = 10.d0**(-3); tol = 10.d0**(-8)

   call base(s)
!--------------------------------------------------------------
   do

      H = 2.0d0
      j = 0

      write(*,*) 'Entre com T:'
      read*, T

      write(*,*) 'Entre com a fase:'
      read*, state

      m_fe = 1.d0; m_af = -1.d0

      ! - - - - - - - - - - - - - - - - - - - - - - -

      call HAM_INTRA(J2,H,s,H_intra)

      call mag_vetor(state,m_fe,m_af,m)

      call Ham_inter_state(state,J2,0.d0,s,m,H_inter)

      H_total = H_intra + H_inter

      call partition(H_total,T,Z)

      call magnetization(H_total,Z,s,1,4,T,m_fe)

      call magnetization(H_total,Z,s,2,3,T,m_af)

      call mag_vetor(state,m_fe,m_af,m)

      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_T_H-F-m.dat")

      do while (H<=20.d0)

         call HAM_INTRA(J2,H,s,H_intra)

         erro1 = 1.d0; erro2 = 1.d0

         do while(max(erro1,erro2) >= tol)

            call Ham_inter_state(state,J2,0.d0,s,m,H_inter)

            H_total = H_intra + H_inter

            call partition(H_total,T,Z)

            call magnetization(H_total,Z,s,1,4,T,m_fe)

            call magnetization(H_total,Z,s,2,3,T,m_af)


            erro1 = abs(m_fe - m(1))
            erro2 = abs(m_af - m(2))

            !  print*, m_fe, m_af, m_order
            !  read(*,*)


            call mag_vetor(state,m_fe,m_af,m)

            m_order = abs(m_fe-m_af)/2.d0


         end do

         call F_helm(T,Z,F)

         print*, H

         write(20,*) H, F, m_order
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '------------'
               write(*,18) H, T
18             format ((F8.5))
               j = 1
            end if
         end if



         H = H + step

      end do

      close(20)

   enddo



end program


