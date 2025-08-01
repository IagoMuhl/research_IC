program square_T
   use CMF
   implicit none

   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total ,s_z, H_inter_prime, H_total_prime
   real*8:: m(6), error(6), mag_prev(6), m_prime(6)
   real*8:: J2, J3, erro, Alfa, Beta, mag
   real*8:: T,H,step,Z,Z_prime,m_order,m_second,tol,F,T_max,F_prime
   character(len=3):: state
   character(len=5):: temp
   integer:: j, cd, i, p



   tol = 10.d0**(-8); J2 = 1.0d0; J3 = -0.1d0; s_z = 0; H = 0.d0
!----------------------------BASE-------------------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
   call HAM_INTRA(J2,J3,s,H_intra)

!--------------------------------------------------------------
   ! write(*,*) 'Entre com H e H_max:'
   ! read*, H,H_max

   write(*,*) 'Entre com a fase:'
   read*, state

   ! ! open(unit=20, file= 'SO_T_' // trim(state) // "_T-H.dat")

   ! write(*,*) 'Entre com o step(-3,-5):'
   ! read*, cd

   T = 2; T_max = 5; cd = -3; step = 10.d0**(cd)
   

   !do

         j = 0; Alfa = 0.d0; Beta = 0.d0

         ! - - - - - - - - - - - - - - - - - - - - - - -

         !open(unit=20, file=trim(state) // "_H_T-F-m.dat")

         WRITE (temp, '(F5.3)') H

         open(unit=20, file= 'Campo_' // trim(temp) // '_Var-Fase_' // trim(state) // "_T-H.dat")
         open(21)

         do while (T/=T_max)

            error = 1.d0; erro = 1.d0

            call mag_vetor(state,m)

            do i = 1,6
            m_prime(i) = -m(i)
            enddo

            do while(erro >= tol)

               call Ham_inter(J2,J3,s,m,m_prime,H_inter)

               call Ham_inter(J2,J3,s,m_prime,m,H_inter_prime)

               H_total = H_intra + H_inter - H*s_z
               H_total_prime = H_intra + H_inter_prime - H*s_z

               !---------------------- SHIFT DA HAMILTONIANA ----------------
               if (T<=10.d0**(-1)) then

                  Alfa = minval(H_total)
                  Beta = minval(H_total_prime)

                  H_total = H_total - Alfa
                  H_total_prime = H_total_prime - Beta

               endif
               !---------------------- SHIFT DA HAMILTONIANA ----------------

               ! print*, m
               ! read*,
               ! print*, m_prime
               ! read*,

               call partition(H_total,T,Z)
               call partition(H_total,T,Z_prime)

               do i = 1, 6

                  mag_prev(i) = m(i)

                  call magnetization(H_total,Z,s,i,T,m(i))

                  error(i) = abs(mag_prev(i) - m(i))

                  ! m1 = m1*0.5 + 0.5*m(1)

                  erro = maxval(error)

               end do

               do i = 1, 6

                  mag_prev(i) = m_prime(i)

                  call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

                  error(i) = abs(mag_prev(i) - m_prime(i))

                  ! m1 = m1*0.5 + 0.5*m(1)

                  erro = maxval(error)

               end do


            end do

            call order_parameter(state,m,m_order)
            call order_parameter(state,m_prime,m_second)


            call F_helm(T,Z,F)
            call F_helm(T,Z_prime,F_prime)

            F = F + Alfa
            F_prime = F_prime + Beta

            F = (F + F_prime)/2

            mag = sum(m)

            write(20,*) T, F, m_order, m_second, mag
            write(21,*) m


!             if (j==0) then
!                if (m_order<=10.d0**(-4)) then
!                   print*, '\/---------\/'
!                   write(*,18) H, T_inicial
!                   print*, '/\---------/\'
! 18                format ((F8.5))
!                   j = 1
!                   !write(20,*) T_inicial, H
!                end if
!             end if



            T = T + step

            if ((max(T,T_max))-(min(T,T_max))<=abs(step)) then
               T = T_max
            endif

            ! if (max(T_final,T_inicial)>=12.d0) then
            !    exit
            ! endif

         end do




         !----------------- NOVAS ENTRADAS  -----------------


         print*, '------------'
         Print*, 'State =','', state
         print*, 'H =','', H
         print*, 'T','',T
         print*, 'J2','',J2
         print*, 'J3','',J3
         print*, '----END-----'


      !enddo


      print*, '=== FIM ==='

      close(20)
      close(21)
      stop

   !enddo

end program


