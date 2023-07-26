program square_T
   use CMF
   implicit none
   real*8, parameter:: J2 = -1.d0
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total
   real*8, dimension(num_sites):: m
   real*8:: H,T,m1,m2,step,Z,m_order,tol,F,erro1,erro2
   character(len=3):: state
   integer:: j

   step = 10.d0**(-3); tol = 10.d0**(-8)

   call base(s)
!--------------------------------------------------------------
   do

      T = 4.d0
      j = 0

      write(*,*) 'Entre com o campo:'
      read*, H

      write(*,*) 'Entre com a fase:'
      read*, state

      m1 = 1.d0; m2 = -1.d0

      ! - - - - - - - - - - - - - - - - - - - - - - -

      call HAM_INTRA(J2,H,s,H_intra)

      call mag_vetor(state,m1,m2,m)

      call Ham_inter_state(state,J2,0.d0,s,m,H_inter)

      H_total = H_intra + H_inter

      call partition(H_total,T,Z)

      call magnetization(H_total,Z,s,1,4,T,m1)

      call magnetization(H_total,Z,s,2,3,T,m2)

      call mag_vetor(state,m1,m2,m)

      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_H_T-F-m.dat")

      do while (T<=14.d0)

         erro1 = 1.d0; erro2 = 1.d0

         do while(max(erro1,erro2) >= tol)

            call Ham_inter_state(state,J2,0.d0,s,m,H_inter)

            H_total = H_intra + H_inter

            call partition(H_total,T,Z)

            call magnetization(H_total,Z,s,1,4,T,m1)

            call magnetization(H_total,Z,s,2,3,T,m2)


            erro1 = abs(m1 - m(1))
            erro2 = abs(m2 - m(2))

            !  print*, m1, m2, m_order
            !  read(*,*)


            call mag_vetor(state,m1,m2,m)

            m_order = abs(m1-m2)/2.d0



         end do

         call F_helm(T,Z,F)

         write(20,*) T, F, m_order
         ! print*, T, m_order, m1, m2

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '------------'
               write(*,18) H, T
18             format ((F8.5))
               j = 1
            end if
         end if



         T = T + step

      end do

      close(20)

   enddo



end program


