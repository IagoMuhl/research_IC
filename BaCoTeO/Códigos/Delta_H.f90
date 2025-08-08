program hexa_H
   use CMF
   implicit none

   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total ,s_z, H_inter_prime, H_total_prime
   real*8:: m(6), error(12), mag_prev(12), m_prime(6)
   real*8:: J2, J3, erro, Alfa, Beta, mag, k
   real*8:: T,H,step,Z,Z_prime,m_order,m_second,tol,F,H_final,F_prime
   character(len=3):: state
   character(len=5):: temp
   integer:: j, cd, i, p



   tol = 10.d0**(-8); J2 = 0.806122449; J3 = -0.111564626; s_z = 0; T = 0.01d0
!----------------------------BASE-------------------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
   call HAM_INTRA(J2,J3,s,H_intra)

!--------------------------------------------------------------
   ! write(*,*) 'Entre com T e T:'
   ! read*, T,T

   ! ! open(unit=20, file= 'SO_T_' // trim(state) // "_T-H.dat")

   ! write(*,*) 'Entre com o step(-3,-5):'
   ! read*, cd

   WRITE (temp, '(F5.3)') T

   do

   H = 2; H_final = 5; cd = -3; k = 0

      if (H < H_final) then
         step = 10.d0**(cd)
      else
         step = -10.d0**(cd)
      endif

      write(*,*) 'Entre com a fase:'
         read*, state

      select case(state) 

         case('pm')
!-------------------------- FASE PM ------------------------------------

      open(unit=30, file= 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
      open(unit=21, file = 'Mag_'// trim(state) // '.dat')


         Alfa = 0.d0

         call mag_vetor(state,m)

         ! - - - - - - - - - - - - - - - - - - - - - - -

         do while (H/=H_final)

            error = 1.d0; erro = 1.d0; 

            do while(erro >= tol)

               call Ham_inter(J2,J3,s,m,m,H_inter)

               H_total = H_intra + H_inter - H*s_z

               !---------------------- SHIFT DA HAMILTONIANA ----------------
               if (T<=10.d0**(-1)) then

                  Alfa = minval(H_total)

                  H_total = H_total - Alfa

               endif
               !---------------------- SHIFT DA HAMILTONIANA ----------------

               call partition(H_total,T,Z)

               do i = 1, 6

                  mag_prev(i) = m(i)

                  call magnetization(H_total,Z,s,i,T,m(i))

                  error(i) = abs(mag_prev(i) - m(i))
                  error(i+6) = 0

                  m(i) = 0.5*m(i) + 0.5*mag_prev(i)

               end do

               erro = maxval(error)

            ! print*, m
            ! print*, erro
            ! read*,
            
            k = k + abs(step)

            if (k>=10**(-(cd+1))) then
            print*, 'Rodando... H = ', H
            k = 0
            endif

            end do

            call order_parameter(state,m,m_order)

            call F_helm(T,Z,F)

            F = F + Alfa

            mag = sum(m)

            write(30,*) H, F, m_order, mag
            write(21,*) H, m


            if (k>=10**(-(cd+1))) then
               print*, 'Rodando... H = ', H
               k = 0
            endif

             H = H + step; k = k + 1

            if ((max(H,H_final))-(min(H,H_final))<=abs(step)) then
               H = H_final
            endif

         end do

         !----------------- NOVAS ENTRADAS  -----------------

         print*, '------------'
         Print*, 'State = ','', state
         print*, 'H =','', H
         print*, 'T','',T
         print*, 'J2','',J2
         print*, 'J3','',J3
         print*, '----END-----'


      !enddo

      print*, '=== FIM ==='

      close(30)
      close(21)
!----------------------------------------------------------------------------

!          case('6')
! !-------------------------- FASE 6 ------------------------------------

!       open(unit=25, file= 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
!       open(unit=21, file = 'Magnetização_cluster_1_'// trim(state) // '.dat')
!       open(unit=22, file = 'Magnetização_cluster_2_'// trim(state) // '.dat')

!          Alfa = 0.d0

!          call mag_vetor(state,m)

!          ! - - - - - - - - - - - - - - - - - - - - - - -

!          ! do while (T/=T)
!          do while (T/=T)
!             error = 1.d0; erro = 1.d0; 

!             do while(erro >= tol)

!                ! print*, m
!                ! read*,

!                call Ham_inter(J2,J3,s,m,m,H_inter)

!                H_total = H_intra + H_inter - H*s_z

!             !---------------------- SHIFT DA HAMILTONIANA ----------------
!                if (T<=10.d0**(-1)) then

!                   Alfa = minval(H_total)

!                   H_total = H_total - Alfa

!                endif
!             !---------------------- SHIFT DA HAMILTONIANA ----------------

!                !call partition(H_total,T,Z)
!                call partition(H_total,T,Z)

!                do i = 1, 6

!                   mag_prev(i) = m(i)

!                   ! call magnetization(H_total,Z,s,i,T,m(i))
!                   call magnetization(H_total,Z,s,i,T,m(i))

!                   error(i) = abs(mag_prev(i) - m(i))
!                   error(i+6) = 0

!                end do

!                erro = maxval(error)

!             end do

!             call order_parameter(state,m,m_order)

!             ! call F_helm(T,Z,F)
!             call F_helm(T,Z,F)

!             F = F + Alfa

!             mag = sum(m)

!             ! write(25,*) T, F, m_order, mag
!             ! write(21,*) T, m
!             ! write(22,*) T, m_prime

!             write(25,*) T, F, m_order, mag
!             write(21,*) T, m
!             write(22,*) T, m_prime


!             if (k==10**(-(cd+1))) then
!                print*, 'Rodando... T = ', T
!                k = 0
!             endif

!             ! T = T + step; k = k + 1
!             T = T - step; k = k + 1

!             if ((max(T,T))-(min(T,T))<=abs(step)) then
!                T = T
!             endif

!          end do

!          !----------------- NOVAS ENTRADAS  -----------------

!          print*, '------------'
!          Print*, 'State = ','', state
!          print*, 'H =','', H
!          print*, 'T','',T
!          print*, 'J2','',J2
!          print*, 'J3','',J3
!          print*, '----END-----'

!       !enddo

!       print*, '=== FIM ==='

!       close(25)
!       close(21)
!       close(22)
! !----------------------------------------------------------------------------

         case('5')
!-------------------------- FASE 5 -----------------------------------
      open(unit=20, file = 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
      open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')
      open(unit=22, file = 'Mag_2_'// trim(state) // '.dat')

         Alfa = 0.d0; Beta = 0.d0

         call mag_vetor(state,m)

         m_prime = -m

         j = 0

         ! - - - - - - - - - - - - - - - - - - - - - - -

         do while (H/=H_final)

            error = 1.d0; erro = 1.d0; 

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


               call partition(H_total,T,Z)
               call partition(H_total_prime,T,Z_prime)

               do i = 1, 6

                  mag_prev(i) = m(i)
                  mag_prev(i+6) = m_prime(i)

                  call magnetization(H_total,Z,s,i,T,m(i))
                  call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

                  error(i) = abs(mag_prev(i) - m(i))

                  error(i+6) = abs(mag_prev(i+6) - m_prime(i))

                  m(i) = 0.5*m(i) + 0.5*mag_prev(i)
                  m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev(i+6)

               end do

            erro = maxval(error)

            k = k + abs(step)

            if (k>=10**(-(cd+1))) then
               print*, 'Rodando... H = ', H
               k = 0
            endif

            end do

            call order_parameter(state,m,m_order)
            call order_parameter(state,m_prime,m_second)


            call F_helm(T,Z,F)
            call F_helm(T,Z_prime,F_prime)

            F = F + Alfa
            F_prime = F_prime + Beta

            F = (F + F_prime)/2

            mag = sum(m)

            write(21,*) H, m
            write(22,*) H, m_prime
            write(20,*) H, F, m_order, mag


            if (j==0) then
               if (m_order<=10.d0**(-4)) then
                  print*, '============='
                  print*, 'Salto parâmetro de ordem (H-T)'
                  print*, '============='
                  write(*,11) (H-step/2), T
   11               format ((F8.5))
                  j = 1
                  !write(20,*) T_inicial, H
               end if
            end if

            if (k>=10**(-(cd+1))) then
               print*, 'Rodando... H = ', H
               k = 0
            endif

            H = H + step; k = k + 1.d0


            if ((max(H,H_final))-(min(H,H_final))<=abs(step)) then
               H = H_final
            endif

         end do


         !----------------- NOVAS ENTRADAS  -----------------


         print*, '============='
         Print*, 'State = ','', state
         print*, 'H =','', H
         print*, 'T =','',T
         print*, 'J2 =','',J2
         print*, 'J3 =','',J3
         print*, '==== END ===='


      !enddo


      close(20)
      close(21)
      close(22)
!-------------------------------------------------------------------------      

         case('4')
!-------------------------- FASE 4 -----------------------------------
      open(unit=28, file = 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
      open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')
      open(unit=22, file = 'Mag_2_'// trim(state) // '.dat')

         Alfa = 0.d0; Beta = 0.d0

         call mag_vetor(state,m)

            do i = 1,6
               m_prime(i) = 1.d0
            enddo

         j = 0
         ! - - - - - - - - - - - - - - - - - - - - - - -

         do while (H/=H_final)
         !   do while (T/=T)

            error = 1.d0; erro = 1.d0; 

            do while(erro >= tol)

               ! print*, m
               ! read*,

               call Ham_inter_fase_4(J2,J3,s,m,m_prime,H_inter,H_inter_prime)

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


               call partition(H_total,T,Z)
               call partition(H_total_prime,T,Z_prime)


               do i = 1, 6

                  mag_prev(i) = m(i)
                  mag_prev(i+6) = m_prime(i)

                  call magnetization(H_total,Z,s,i,T,m(i))
                  call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

                  error(i) = abs(mag_prev(i) - m(i))

                  error(i+6) = abs(mag_prev(i+6) - m_prime(i))

                  m(i) = 0.5*m(i) + 0.5*mag_prev(i)
                  m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev(i+6)

               end do

               ! mag_prev(1) = m(1)
               ! mag_prev(2) = m(2)
               ! mag_prev(7) = m_prime(1)

               ! call magnetization(H_total,Z,s,1,T,m(1))
               ! call magnetization(H_total,Z,s,2,T,m(2))
               ! call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(1))

               ! error(1) = abs(mag_prev(1) - m(1))
               ! error(2) = abs(mag_prev(2) - m(2))
               ! error(7) = abs(mag_prev(7) - m_prime(1))

               ! m(4) = m(1)
               ! m(3) = m(2); m(5) = m(2); m(6) = m(2)
               
               ! do i = 2,6
               !    m_prime(i) = m_prime(1)
               ! enddo

            erro = maxval(error)

            ! print*, m
            ! print*, m_prime
            ! print*, erro
            ! print*, k
            ! read*,

            k = k + abs(step)

            if (k>=10**(-(cd+1))) then
               print*, 'Rodando... H = ', H, erro
               k = 0
            endif

            end do

            call order_parameter(state,m,m_order)
            call order_parameter(state,m_prime,m_second)


            call F_helm(T,Z,F)
            call F_helm(T,Z_prime,F_prime)

            F = F + Alfa
            F_prime = F_prime + Beta

            F = (3*F + 1*F_prime)/4

            mag = sum(m)

            write(28,*) H, F, m_order, mag
            write(21,*) H, m
            write(22,*) H, m_prime
            

            if (j==0) then
               if (m_order<=10.d0**(-4)) then
                  print*, '============='
                  print*, 'Salto parâmetro de ordem (H-T)'
                  print*, '============='
                  write(*,10) (H-step/2), T
   10                format ((F8.5))
                  j = 1
                  !write(20,*) T_inicial, H
               end if
            end if

            if (k>=10**(-(cd+1))) then
               print*, 'Rodando... H = ', H
               k = 0
            endif

            H = H + step; k = k + 1.d0

            if ((max(H,H_final))-(min(H,H_final))<=abs(step)) then
               H = H_final
            endif

         end do


         !----------------- NOVAS ENTRADAS  -----------------


         print*, '============='
         Print*, 'State = ','', state
         print*, 'H =','', H
         print*, 'T =','',T
         print*, 'J2 =','',J2
         print*, 'J3 =','',J3
         print*, '==== END ===='


      !enddo


      close(28)
      close(21)
      close(22)
!-------------------------------------------------------------------------      


         case('fim')

            stop

         case default
            write(*,*) 'Inaccurate State'
      end select

   enddo

end program


