program hexa_H_varre
   use CMF
   implicit none

   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total, s_z, H_inter_prime, H_total_prime
   real*8:: m(num_sites), error(2*num_sites), mag_prev, m_prime(num_sites), mag_prev_prime
   real*8:: J2, J3, erro, Alfa, Beta, mag, passo
   real*8:: T,H,step,Z,Z_prime,tol,F,H_final,F_prime,T_final
   character(len=3):: state
   character(len=5):: temp
   integer:: j, cd, i, p, n, iter, iter_max



   tol = 10.d0**(-8); J2 = 0.806122449; J3 = -0.111564626
   s_z = 0; cd = -3; iter_max = 500

   T = 0.25d0; T_final = 0.95d0
!----------------------------BASE-------------------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
   call HAM_INTRA(J2,J3,s,H_intra)

!--------------------------------------------------------------
   ! write(*,*) 'Entre com T:'
   ! read*, T

   ! ! open(unit=20, file= 'SO_T_' // trim(state) // "_T-H.dat")

   ! write(*,*) 'Entre com o step(-3,-5):'
   ! read*, cd

   passo = 0.1


   do while (T/=T_final)

   WRITE (temp, '(F5.3)') T

   ! write(*,*) 'Entre com H e H_final:'
   ! read*, H,H_final

      ! H = 7.9; H_final = 6;

      ! if (H < H_final) then
      !    step = 10.d0**(cd)
      ! else
      !    step = -10.d0**(cd)
      ! endif

!       state = 'pm'; iter = 0
! !-------------------------- FASE PM ------------------------------------

!       open(unit=30, file= 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
!       open(unit=21, file = 'Mag_'// trim(state) // '.dat')


!          Alfa = 0.d0

!          call mag_vetor(state,m)

!          ! - - - - - - - - - - - - - - - - - - - - - - -

!          do while (H/=H_final)

!                erro = 1.d0; error = 0.d0

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

!                do i = 1, num_sites

!                   mag_prev = m(i)

!                   call magnetization(H_total,Z,s,i,T,m(i))

!                   error(i) = abs(mag_prev - m(i))
!                   error(i+num_sites) = 0
                  
!                   m(i) = 0.5*m(i) + mag_prev*0.5

!                end do

!               iter = iter + 1
!             erro = maxval(error)

!             if (iter>=iter_max) then
!                exit
!             endif
!             end do

!             if (iter>=iter_max) then
!                exit
!             endif

!             print*, iter
!             iter = 0

!             call F_helm(T,Z,F)

!             F = F + Alfa

!             mag = sum(m)/num_sites

!             print*, 'Temp:',T, 'Campo:',H, 'Mag:',mag ,'Fase: ',state
!             write(30,*) H, F, mag
!             write(21,*) H, m


!              H = H + step

!             if ((max(H,H_final))-(min(H,H_final))<=abs(step)) then
!                H = H_final
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

!       close(30)
!       close(21)
! ! !----------------------------------------------------------------------------
   
      if (T <= 0.6) then

      H = 3.0; H_final = 2.4;

      if (H < H_final) then
         step = 10.d0**(cd)
      else
         step = -10.d0**(cd)
      endif

      state = '6'; iter = 0
!-------------------------- FASE 6 ------------------------------------

      open(unit=25, file= 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
      ! open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')

         Alfa = 0.d0

         call mag_vetor(state,m)

         ! - - - - - - - - - - - - - - - - - - - - - - -

         do while (H/=H_final)
            error = 0.d0; erro = 1.d0; 

            do while(erro >= tol)

               call Ham_inter(state,J2,J3,s,m,m,H_inter)

               H_total = H_intra + H_inter - H*s_z

            !---------------------- SHIFT DA HAMILTONIANA ----------------
               if (T<=10.d0**(-1)) then

                  Alfa = minval(H_total)

                  H_total = H_total - Alfa

               endif
            !---------------------- SHIFT DA HAMILTONIANA ----------------

               call partition(H_total,T,Z)

                  mag_prev = m(1)

                  call magnetization(H_total,Z,s,1,T,m(1))

                  error(1) = abs(mag_prev - m(1))

                  ! m(1) = 0.5*m(1) + mag_prev*0.5

                  m(2) = m(1)
               
               n = 15

               do i = 3, 10

                  mag_prev = m(i)

                  call magnetization(H_total,Z,s,i,T,m(i))

                  error(i) = abs(mag_prev - m(i))

                  ! m(i) = 0.5*m(i) + mag_prev*0.5

                  m(i+n) = m(i)
                  n = n - 2

               end do

            iter = iter + 1
            erro = maxval(error)

            if (iter>=iter_max) then
               exit
            endif
            end do

            if (iter>=iter_max) then
               exit
            endif

            print*, iter
            iter = 0

            call F_helm(T,Z,F)

            F = F + Alfa

            mag = sum(m)/num_sites

            print*, 'Temp:',T, 'Campo:',H, 'Mag:',mag ,'Fase: ',state
            write(25,*) H, F, mag
            ! write(21,*) H, m

            H = H + step

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

      close(25)
      ! close(21)

!----------------------------------------------------------------------------
      endif

      H = 2.4; H_final = 2.7;

      if (H < H_final) then
         step = 10.d0**(cd)
      else
         step = -10.d0**(cd)
      endif

      state = '5'; iter = 0
!-------------------------- FASE 5 -----------------------------------
      open(unit=20, file = 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
      ! open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')
      ! open(unit=22, file = 'Mag_2_'// trim(state) // '.dat')

         Alfa = 0.d0; Beta = 0.d0

         call mag_vetor(state,m)

         j = 0

         m_prime = - m

         ! - - - - - - - - - - - - - - - - - - - - - - -

         do while (H/=H_final)

            error = 1.d0; erro = 1.d0; 

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

                  n = 7

                  error = 0.d0

               do i = 1, 4

                  mag_prev = m(i)
                  mag_prev_prime = m_prime(i)

                  call magnetization(H_total,Z,s,i,T,m(i))
                  call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

                  error(i) = abs(mag_prev - m(i))
                  error(i+num_sites) = abs(mag_prev_prime - m_prime(i))

                  ! m(i) = 0.5*m(i) + 0.5*mag_prev
                  ! m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev_prime

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

                  ! m(i) = 0.5*m(i) + 0.5*mag_prev
                  ! m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev_prime

                  m(i+n) = m(i)
                  m_prime(i+n) = m_prime(i)

                  n = n - 2
               end do

            iter = iter + 1
            erro = maxval(error)

            if (iter>=iter_max) then
               exit
            endif
            end do

            if (iter>=iter_max) then
               exit
            endif

            print*, iter
            iter = 0


            call F_helm(T,Z,F)
            call F_helm(T,Z_prime,F_prime)

            F = F + Alfa
            F_prime = F_prime + Beta

            ! write(21,*) H, m, F
            ! write(22,*) H, m_prime, F_prime

            F = (F + F_prime)/2

            mag = (sum(m) + sum(m_prime))/(2*num_sites)


            print*, 'Temp:',T, 'Campo:',H, 'Mag:',mag ,'Fase: ',state
            write(20,*) H, F, mag


   !          if (j==0) then
   !             if (m_order<=10.d0**(-4)) then
   !                print*, '============='
   !                print*, 'Salto parâmetro de ordem (H-T)'
   !                print*, '============='
   !                write(*,11) (H-step/2), T
   ! 11               format ((F8.5))
   !                j = 1
   !                !write(20,*) T_inicial, H
   !             end if
   !          end if


            H = H + step


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
      ! close(21)
      ! close(22)
!-------------------------------------------------------------------------      
      
      H = 3; H_final = 2.4;

      if (H < H_final) then
         step = 10.d0**(cd)
      else
         step = -10.d0**(cd)
      endif

      state = '4'; iter = 0
!-------------------------- FASE 4 -----------------------------------
      open(unit=28, file = 'Fase_' // trim(state) // '_Temperatura_' // trim(temp) // "_T-H.dat")
      ! open(unit=21, file = 'Mag_1_'// trim(state) // '.dat')
      ! open(unit=22, file = 'Mag_2_'// trim(state) // '.dat')

         Alfa = 0.d0; Beta = 0.d0

         call mag_vetor(state,m)

         m_prime = [1.d0, -1.d0, 1.d0, 1.d0, 1.d0, 1.d0, -1.d0, 1.d0, 1.d0, &
                  & -1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, -1.d0, 1.d0]

         j = 0
         ! - - - - - - - - - - - - - - - - - - - - - - -

         do while (H/=H_final)

            error = 0.d0; erro = 1.d0; 

            do while(erro >= tol)

               call Ham_inter(state,J2,J3,s,m,m_prime,H_inter)

               H_inter_prime = 0.d0

               do i = 1, maxConfig

                  H_inter_prime(i) = J1*((s(i,1)-m_prime(1)/2)*(m_prime(18)) + &
                  & + (s(i,2)-m_prime(2)/2)*(m_prime(15)) + (s(i,3)-m_prime(3)/2)*(m(14)) &
                  & + (s(i,6)-m_prime(6)/2)*(m(13)) + (s(i,7)-m_prime(7)/2)*(m_prime(12)) &
                  & + (s(i,8)-m_prime(8)/2)*(m_prime(9)) &
                  & + (s(i,9)-m_prime(9)/2)*(m_prime(8)) + (s(i,12)-m_prime(12)/2)*(m_prime(7)) &
                  & + (s(i,13)-m_prime(13)/2)*(m(6)) &
                  & + (s(i,14)-m_prime(14)/2)*(m(3)) + (s(i,15)-m_prime(15)/2)*(m_prime(2)) &
                  & + (s(i,18)-m_prime(18)/2)*(m_prime(1))) &
                  & + J2*((s(i,1)-m_prime(1)/2)*(m_prime(1)+m_prime(1)+m_prime(17)+m_prime(15)) &
                  & + (s(i,2)-m_prime(2)/2)*(m_prime(18)+m_prime(16)+m_prime(14)+m(14)) &
                  & + (s(i,3)-m_prime(3)/2)*(m(15)+m_prime(15)+m(13)) + (s(i,4)-m_prime(4)/2)*(m(14)) &
                  & + (s(i,5)-m_prime(5)/2)*(m(13)) + (s(i,6)-m_prime(6)/2)*(m(14)+m_prime(12)+m(12)) &
                  & + (s(i,7)-m_prime(7)/2)*(m(13)+m_prime(13)+m_prime(11)+m_prime(9)) &
                  & + (s(i,8)-m_prime(8)/2)*(m_prime(12)+m_prime(10)+m_prime(8)+m_prime(8)) &
                  & + (s(i,9)-m_prime(9)/2)*(m_prime(9)+m_prime(9)+m_prime(7)) &
                  & + (s(i,10)-m_prime(10)/2)*(m_prime(8)) + (s(i,11)-m_prime(11)/2)*(m_prime(7)) &
                  & + (s(i,12)-m_prime(12)/2)*(m_prime(8)+m(6)+m_prime(6)) &
                  & + (s(i,13)-m_prime(13)/2)*(m(7)+m_prime(7)+m(5)+m(3)) &
                  & + (s(i,14)-m_prime(14)/2)*(m(6)+m(4)+m_prime(2)+m(2)) &
                  & + (s(i,15)-m_prime(15)/2)*(m_prime(3)+m(3)+m_prime(1)) &
                  & + (s(i,16)-m_prime(16)/2)*(m_prime(2)) + (s(i,17)-m_prime(17)/2)*(m_prime(1)) &
                  & + (s(i,18)-m_prime(18)/2)*(m_prime(2)+m_prime(18)+m_prime(18))) &
                  & + J3*((s(i,1)-m_prime(1)/2)*(m_prime(18)+m_prime(16)) &
                  & + (s(i,2)-m_prime(2)/2)*(m_prime(17)+m(15)) &
                  & + (s(i,3)-m_prime(3)/2)*(m_prime(14)) + (s(i,4)-m_prime(4)/2)*(m(13)) &
                  & + (s(i,5)-m_prime(5)/2)*(m(14)) + (s(i,6)-m_prime(6)/2)*(m_prime(13)) &
                  & + (s(i,7)-m_prime(7)/2)*(m(12)+m_prime(10)) &
                  & + (s(i,8)-m_prime(8)/2)*(m_prime(11)+m_prime(9)) &
                  & + (s(i,9)-m_prime(9)/2)*(m_prime(8)) &
                  & + (s(i,10)-m_prime(10)/2)*(m_prime(7)) + (s(i,11)-m_prime(11)/2)*(m_prime(8)) &
                  & + (s(i,12)-m_prime(12)/2)*(m(7)) &
                  & + (s(i,13)-m_prime(13)/2)*(m_prime(6)+m(4)) &
                  & + (s(i,14)-m_prime(14)/2)*(m(5)+m_prime(3)) &
                  & + (s(i,15)-m_prime(15)/2)*(m(2)) &
                  & + (s(i,16)-m_prime(16)/2)*(m_prime(1)) + (s(i,17)-m_prime(17)/2)*(m_prime(2)) &
                  & + (s(i,18)-m_prime(18)/2)*(m_prime(1)))

               enddo

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
               error = 0.d0

               call partition(H_total,T,Z)
               call partition(H_total_prime,T,Z_prime)


                do i = 1, num_sites

                  mag_prev = m(i)
                  mag_prev_prime = m_prime(i)

                  call magnetization(H_total,Z,s,i,T,m(i))
                  call magnetization(H_total_prime,Z_prime,s,i,T,m_prime(i))

                  error(i) = abs(mag_prev - m(i))
                  error(i+num_sites) = abs(mag_prev_prime - m_prime(i))

                  ! m(i) = 0.5*m(i) + 0.5*mag_prev
                  ! m_prime(i) = 0.5*m_prime(i) + 0.5*mag_prev_prime


               end do

               iter = iter + 1
            erro = maxval(error)

            if (iter>=iter_max) then
               exit
            endif
            end do

            if (iter>=iter_max) then
               exit
            endif

            print*, iter
            iter = 0

            call F_helm(T,Z,F)
            call F_helm(T,Z_prime,F_prime)

            F = F + Alfa
            F_prime = F_prime + Beta

            ! write(21,*) H, m, F
            ! write(22,*) H, m_prime, F_prime

            F = (F + 3*F_prime)/4

            mag = (sum(m) + 3*sum(m_prime))/(4*num_sites)

            print*, 'Temp:',T, 'Campo:',H, 'Mag:',mag ,'Fase: ',state
            write(28,*) H, F, mag


   !          if (j==0) then
   !             if (m_order<=10.d0**(-4)) then
   !                print*, '============='
   !                print*, 'Salto parâmetro de ordem (H-T)'
   !                print*, '============='
   !                write(*,10) (H-step/2), T
   ! 10                format ((F8.5))
   !                j = 1
   !                !write(20,*) T_inicial, H
   !             end if
   !          end if

            H = H + step

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
      ! close(21)
      ! close(22)
!-------------------------------------------------------------------------      

      if (T>=T_final) stop
      
      T = T + passo

   enddo

end program


