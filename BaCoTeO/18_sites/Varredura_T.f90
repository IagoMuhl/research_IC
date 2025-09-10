program hexa_T_varre
   use CMF
   implicit none

   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total ,s_z, H_inter_prime, H_total_prime
   real*8:: m(num_sites), error(2*num_sites), mag_prev, m_prime(num_sites), mag_prev_prime
   real*8:: J2, J3, erro, Alfa, Beta, mag, h_final
   real*8:: T,H,step,Z,Z_prime,m_order,tol,F,T_max,F_prime, passo
   character(len=3):: state
   character(len=5):: temp
   integer:: j, cd, i, p, k, n



   tol = 10.d0**(-8); J2 = 0.806122449; J3 = -0.111564626; s_z = 0; 

   H = 0.d0; H_final = 0.05d0
!----------------------------BASE-------------------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
   call HAM_INTRA(J2,J3,s,H_intra)

!--------------------------------------------------------------
   ! write(*,*) 'Entre com T e T_max:'
   ! read*, T,T_max

   ! ! open(unit=20, file= 'SO_T_' // trim(state) // "_T-H.dat")

   ! write(*,*) 'Entre com o step(-3,-5):'
   ! read*, cd

   do while (H/=H_final)

   ! write(*,*) 'Entre com T e T_max:'
   ! read*, T,T_max

   passo = 10.d0**(-1);

   cd = -5; step = 10.d0**(cd); k = 0

      WRITE (temp, '(F5.3)') H

      ! if (T < T_max) then
      !    step = 10.d0**(cd)
      ! else
      !    step = -10.d0**(cd)
      ! endif

      ! write(*,*) 'Entre com a fase:'
      ! read*, state

      

!       state = 'pm'
! !-------------------------- FASE PM ------------------------------------

!       open(unit=30, file= 'Fase_' // trim(state) // '_Campo_' // trim(temp) // "_T-H.dat")
!       open(unit=21, file = 'Mag_'// trim(state) // '.dat')


!          Alfa = 0.d0

!          call mag_vetor(state,m)

!          ! - - - - - - - - - - - - - - - - - - - - - - -

!          do while (T/=T_max)

!             error = 0.d0; erro = 1.d0; 

!             do while(erro >= tol)

!                call Ham_inter(state,J2,J3,s,m,m,H_inter)

!                H_total = H_intra + H_inter - H*s_z

!                !---------------------- SHIFT DA HAMILTONIANA ----------------
!                if (T<=10.d0**(-1)) then

!                   Alfa = minval(H_total)

!                   H_total = H_total - Alfa

!                endif
!                !---------------------- SHIFT DA HAMILTONIANA ----------------

!                call partition(H_total,T,Z)

!                 do i = 1, num_sites

!                   mag_prev = m(i)

!                   call magnetization(H_total,Z,s,i,T,m(i))

!                   error(i) = abs(mag_prev - m(i))
!                   error(i+num_sites) = 0
                  
!                   m(i) = 0.5*m(i) + mag_prev*0.5

!                end do

!                erro = maxval(error)

!             ! print*, m
!             ! print*, erro
!             ! read*,

!             end do

!             ! call order_parameter(state,m,m_order)

!             call F_helm(T,Z,F)

!             F = F + Alfa

!             mag = sum(m)/num_sites

!             print*, T, mag
!             write(30,*) T, F, mag
!             write(21,*) T, m

!             T = T + step

!             if ((max(T,T_max))-(min(T,T_max))<=abs(step)) then
!                T = T_max
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

!       close(30)
!       close(21)
! !----------------------------------------------------------------------------

!       state = '6'
! !-------------------------- FASE 6 ------------------------------------

!       open(unit=25, file= 'Fase_' // trim(state) // '_Campo_' // trim(temp) // "_T-H.dat")
!       open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')
!       open(unit=22, file = 'Mag_2_'// trim(state) // '.dat')

!          Alfa = 0.d0

!          call mag_vetor(state,m)

!          ! - - - - - - - - - - - - - - - - - - - - - - -

!          do while (T/=T_max)
!             error = 0.d0; erro = 1.d0; 

!             do while(erro >= tol)


!                call Ham_inter(state,J2,J3,s,m,m,H_inter)

!                H_total = H_intra + H_inter - H*s_z

!             !---------------------- SHIFT DA HAMILTONIANA ----------------
!                if (T_max<=10.d0**(-1)) then

!                   Alfa = minval(H_total)

!                   H_total = H_total - Alfa

!                endif
!             !---------------------- SHIFT DA HAMILTONIANA ----------------

!                call partition(H_total,T,Z)

!                ! do i = 1, num_sites

!                !    mag_prev = m(i)

!                !    call magnetization(H_total,Z,s,i,T,m(i))

!                !    error(i) = abs(mag_prev - m(i))
!                !    error(i+num_sites) = 0

!                ! end do

!                mag_prev = m(1)

!                   call magnetization(H_total,Z,s,1,T,m(1))

!                   error(1) = abs(mag_prev - m(1))

!                   m(1) = 0.5*m(1) + mag_prev*0.5

!                   m(2) = m(1)
               
!                n = 15

!                do i = 3, 10

!                   mag_prev = m(i)

!                   call magnetization(H_total,Z,s,i,T,m(i))

!                   error(i) = abs(mag_prev - m(i))

!                   m(i) = 0.5*m(i) + mag_prev*0.5

!                   m(i+n) = m(i)
!                   n = n - 2

!                end do

!                erro = maxval(error)

!             end do


!             call F_helm(T,Z,F)

!             F = F + Alfa

!             mag = sum(m)/num_sites

!             print*, T, mag
!             write(25,*) T, F, mag
!             write(21,*) T, m
!             write(22,*) T, m_prime



!             T = T + step


!             if ((max(T,T_max))-(min(T,T_max))<=abs(step)) then
!                T = T_max
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
      T = 1.7; T_max = 2.2

      if (T < T_max) then
         step = 10.d0**(cd)
      else
         step = -10.d0**(cd)
      endif

      state = '5'
!-------------------------- FASE 5 -----------------------------------
      open(unit=20, file = 'Fase_' // trim(state) // '_Campo_' // trim(temp) // "_T-H.dat")
      ! open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')
      ! open(unit=22, file = 'Mag_2_'// trim(state) // '.dat')

         Alfa = 0.d0; Beta = 0.d0

         call mag_vetor(state,m)

         m_prime = - m

         j = 0

         ! - - - - - - - - - - - - - - - - - - - - - - -

         do while (T/=T_max)

            error = 0.d0; erro = 1.d0; 

            do while(erro >= tol)

               call Ham_inter(state,J2,J3,s,m,m_prime,H_inter)

               call Ham_inter(state,J2,J3,s,m_prime,m,H_inter_prime)

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

               ! do i = 1, num_sites

               !    mag_prev = m(i)
               !    mag_prev_prime = m_prime(i)

               !    call magnetization(H_total,Z,s,i,T,m(i))
               !    call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

               !    error(i) = abs(mag_prev - m(i))

               !    error(i+num_sites) = abs(mag_prev_prime - m_prime(i))

               ! end do
               n = 7

                  error = 0.d0

               do i = 1, 4

                  mag_prev = m(i)
                  mag_prev_prime = m_prime(i)

                  call magnetization(H_total,Z,s,i,T,m(i))
                  call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

                  error(i) = abs(mag_prev - m(i))
                  error(i+num_sites) = abs(mag_prev_prime - m_prime(i))

                  m(i) = 0.5*m(i) + 0.5*mag_prev
                  m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev_prime

                  m(i+n) = m(i)
                  m_prime(i+n) = m_prime(i)

                  n = n - 2
               end do

                  n = 9

               do i = 9, 13

                  mag_prev = m(i)
                  mag_prev_prime = m_prime(i)

                  call magnetization(H_total,Z,s,i,T,m(i))
                  call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

                  error(i) = abs(mag_prev - m(i))
                  error(i+num_sites) = abs(mag_prev_prime - m_prime(i))

                  m(i) = 0.5*m(i) + 0.5*mag_prev
                  m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev_prime

                  m(i+n) = m(i)
                  m_prime(i+n) = m_prime(i)

                  n = n - 2
               end do

            erro = maxval(error)

            end do


            call F_helm(T,Z,F)
            call F_helm(T,Z_prime,F_prime)

            F = F + Alfa
            F_prime = F_prime + Beta

            F = (F + F_prime)/2

            m_order = ((m(1)-m(2)-m(3)+m(4)-m(9)-m(10)+m(11)-m(12)+m(13))&
            & + (-m_prime(1)+m_prime(2)+m_prime(3)-m_prime(4)+m_prime(9)+m_prime(10)&
            & -m_prime(11)+m_prime(12)-m_prime(13)))/num_sites

            mag = (sum(m) + sum(m_prime))/(2*num_sites)

            Print*, T, mag, m_order
            ! write(21,*) T, m
            ! write(22,*) T, m_prime
            write(20,*) T, F, mag, m_order


            if (j==0) then
               if (m_order<=10.d0**(-4)) then
                  print*, '============='
                  print*, 'Salto parâmetro de ordem (H-T)'
                  print*, '============='
                  write(*,11) H, (T-step/2)
   11               format ((F8.5))
                  j = 1
                  !write(20,*) T_inicial, H
               end if
            end if


            T = T + step


            if ((max(T,T_max))-(min(T,T_max))<=abs(step)) then
               T = T_max
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
      ! close(21)
      ! close(22)
!-------------------------------------------------------------------------      

!       state = '4'
! !-------------------------- FASE 4 -----------------------------------
!       open(unit=28, file = 'Fase_' // trim(state) // '_Campo_' // trim(temp) // "_T-H.dat")
!       open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')
!       open(unit=22, file = 'Mag_2_'// trim(state) // '.dat')

!          Alfa = 0.d0; Beta = 0.d0

!          call mag_vetor(state,m)

!          m_prime = [1.d0, -1.d0, 1.d0, 1.d0, 1.d0, 1.d0, -1.d0, 1.d0, 1.d0, &
!                   & -1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, -1.d0, -1.d0]

!          j = 0
!          ! - - - - - - - - - - - - - - - - - - - - - - -

!          do while (T/=T_max)
!          !   do while (T_max/=T)

!             error = 0.d0; erro = 1.d0; 

!             do while(erro >= tol)

!                call Ham_inter(state,J2,J3,s,m,m_prime,H_inter)

!                H_inter_prime = 0.d0

!                do i = 1, maxConfig

!                   H_inter_prime(i) = J1*((s(i,1)-m_prime(1)/2)*(m_prime(18)) + &
!                   & + (s(i,2)-m_prime(2)/2)*(m_prime(15)) + (s(i,3)-m_prime(3)/2)*(m(14)) &
!                   & + (s(i,6)-m_prime(6)/2)*(m(13)) + (s(i,7)-m_prime(7)/2)*(m_prime(12)) &
!                   & + (s(i,8)-m_prime(8)/2)*(m_prime(9)) &
!                   & + (s(i,9)-m_prime(9)/2)*(m_prime(8)) + (s(i,12)-m_prime(12)/2)*(m_prime(7)) &
!                   & + (s(i,13)-m_prime(13)/2)*(m(6)) &
!                   & + (s(i,14)-m_prime(14)/2)*(m(3)) + (s(i,15)-m_prime(15)/2)*(m_prime(2)) &
!                   & + (s(i,18)-m_prime(18)/2)*(m_prime(1))) &
!                   & + J2*((s(i,1)-m_prime(1)/2)*(m_prime(1)+m_prime(1)+m_prime(17)+m_prime(15)) &
!                   & + (s(i,2)-m_prime(2)/2)*(m_prime(18)+m_prime(16)+m_prime(14)+m(14)) &
!                   & + (s(i,3)-m_prime(3)/2)*(m(15)+m_prime(15)+m(13)) + (s(i,4)-m_prime(4)/2)*(m(14)) &
!                   & + (s(i,5)-m_prime(5)/2)*(m(13)) + (s(i,6)-m_prime(6)/2)*(m(14)+m_prime(12)+m(12)) &
!                   & + (s(i,7)-m_prime(7)/2)*(m(13)+m_prime(13)+m_prime(11)+m_prime(9)) &
!                   & + (s(i,8)-m_prime(8)/2)*(m_prime(12)+m_prime(10)+m_prime(8)+m_prime(8)) &
!                   & + (s(i,9)-m_prime(9)/2)*(m_prime(9)+m_prime(9)+m_prime(7)) &
!                   & + (s(i,10)-m_prime(10)/2)*(m_prime(8)) + (s(i,11)-m_prime(11)/2)*(m_prime(7)) &
!                   & + (s(i,12)-m_prime(12)/2)*(m_prime(8)+m(6)+m_prime(6)) &
!                   & + (s(i,13)-m_prime(13)/2)*(m(7)+m_prime(7)+m(5)+m(3)) &
!                   & + (s(i,14)-m_prime(14)/2)*(m(6)+m(4)+m_prime(2)+m(2)) &
!                   & + (s(i,15)-m_prime(15)/2)*(m_prime(3)+m(3)+m_prime(1)) &
!                   & + (s(i,16)-m_prime(16)/2)*(m_prime(2)) + (s(i,17)-m_prime(17)/2)*(m_prime(1)) &
!                   & + (s(i,18)-m_prime(18)/2)*(m_prime(2)+m_prime(18)+m_prime(18))) &
!                   & + J3*((s(i,1)-m_prime(1)/2)*(m_prime(18)+m_prime(16)) &
!                   & + (s(i,2)-m_prime(2)/2)*(m_prime(17)+m(15)) &
!                   & + (s(i,3)-m_prime(3)/2)*(m_prime(14)) + (s(i,4)-m_prime(4)/2)*(m(13)) &
!                   & + (s(i,5)-m_prime(5)/2)*(m(14)) + (s(i,6)-m_prime(6)/2)*(m_prime(13)) &
!                   & + (s(i,7)-m_prime(7)/2)*(m(12)+m_prime(10)) &
!                   & + (s(i,8)-m_prime(8)/2)*(m_prime(11)+m_prime(9)) &
!                   & + (s(i,9)-m_prime(9)/2)*(m_prime(8)) &
!                   & + (s(i,10)-m_prime(10)/2)*(m_prime(7)) + (s(i,11)-m_prime(11)/2)*(m_prime(8)) &
!                   & + (s(i,12)-m_prime(12)/2)*(m(7)) &
!                   & + (s(i,13)-m_prime(13)/2)*(m_prime(6)+m(4)) &
!                   & + (s(i,14)-m_prime(14)/2)*(m(5)+m_prime(3)) &
!                   & + (s(i,15)-m_prime(15)/2)*(m(2)) &
!                   & + (s(i,16)-m_prime(16)/2)*(m_prime(1)) + (s(i,17)-m_prime(17)/2)*(m_prime(2)) &
!                   & + (s(i,18)-m_prime(18)/2)*(m_prime(1)))

!                enddo

!                H_total = H_intra + H_inter - H*s_z
!                H_total_prime = H_intra + H_inter_prime - H*s_z

!                !---------------------- SHIFT DA HAMILTONIANA ----------------
!                if (T<=10.d0**(-1)) then

!                   Alfa = minval(H_total)
!                   Beta = minval(H_total_prime)

!                   H_total = H_total - Alfa
!                   H_total_prime = H_total_prime - Beta

!                endif
!                !---------------------- SHIFT DA HAMILTONIANA ----------------


!                call partition(H_total,T,Z)
!                call partition(H_total_prime,T,Z_prime)


!                ! do i = 1, num_sites

!                !    mag_prev = m(i)
!                !    mag_prev_prime = m_prime(i)

!                !    call magnetization(H_total,Z,s,i,T,m(i))
!                !    call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))


!                !    error(i) = abs(mag_prev - m(i))
!                !    error(i+num_sites) = abs(mag_prev_prime - m_prime(i))

!                !    m(i) = 0.5*m(i) + 0.5*mag_prev
!                !    m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev_prime


!                ! end do
!                   do i = 1, 5, 2

!                   mag_prev = m(i)
!                   mag_prev_prime = m_prime(i)

!                   call magnetization(H_total,Z,s,i,T,m(i))

!                   error(i) = abs(mag_prev - m(i))

!                   m(i) = 0.5*m(i) + 0.5*mag_prev

!                end do

!                m(2) = m(1)
!                n = 6
!                do i = 1,2
!                   m(i+n) = m(i)
!                   m(i+2*n) = m(i)
!                enddo


!                 do i = 1, num_sites

!                   ! mag_prev = m(i)
!                   mag_prev_prime = m_prime(i)

!                   ! call magnetization(H_total,Z,s,i,T,m(i))
!                   call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

!                   ! error(i) = abs(mag_prev - m(i))
!                   error(i+num_sites) = abs(mag_prev_prime - m_prime(i))

!                   ! m(i) = 0.5*m(i) + 0.5*mag_prev
!                   m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev_prime


!                end do

!             erro = maxval(error)

!             ! print*, m
!             ! print*, m_prime
!             ! print*, erro
!             ! print*, k
!             ! read*,

!             end do

!             ! call order_parameter(state,m,m_order)
!             ! call order_parameter(state,m_prime,m_second)


!             call F_helm(T,Z,F)
!             call F_helm(T,Z_prime,F_prime)

!             ! call F_helm(T_max,Z,F)
!             ! call F_helm(T_max,Z_prime,F_prime)

!             F = F + Alfa
!             F_prime = F_prime + Beta

!             F = (F + 3*F_prime)/4

!             mag = (sum(m) + 3*sum(m_prime))/(4*num_sites)

!             write(28,*) T, F, mag
!             write(21,*) T, m
!             write(22,*) T, m_prime

!             print*, T, mag
            

!             if (j==0) then
!                if (m_order<=10.d0**(-4)) then
!                   print*, '============='
!                   print*, 'Salto parâmetro de ordem (H-T)'
!                   print*, '============='
!                   write(*,10) H, (T-step/2)
!    10                format ((F8.5))
!                   j = 1
!                   !write(20,*) T_inicial, H
!                end if
!             end if



!             T = T + step


!             if ((max(T,T_max))-(min(T,T_max))<=abs(step)) then
!                T = T_max
!             endif

!          end do


!          !----------------- NOVAS ENTRADAS  -----------------


!          print*, '============='
!          Print*, 'State = ','', state
!          print*, 'H =','', H
!          print*, 'T =','',T
!          print*, 'J2 =','',J2
!          print*, 'J3 =','',J3
!          print*, '==== END ===='


!       !enddo


!       close(28)
!       close(21)
!       close(22)
! !-------------------------------------------------------------------------      

   if (H>=H_final) stop
      
      H = H + passo


   enddo

end program


