program square_T
   use CMF
   implicit none

   integer, dimension(:,:),allocatable:: s, s_sub
   real*8, dimension(maxConfig,8):: N
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_z
   real*8:: m(12), error(8), mag_prev(8)
   real*8:: J2, erro, Alfa, passo, T_min, T_max
   real*8:: T_inicial,T_final,H,step,Z,m_order,tol,F,H_max
   character(len=3):: state
   character(len=5):: temp
   integer:: j, cd, i, p, ef
   character(8)  :: date
   character(10) :: time
   character(5)  :: zone
   integer,dimension(8) :: values

   allocate(s(maxConfig,num_sites) , s_sub(maxConfig,12))

   tol = 10.d0**(-8); J2 = -0.33d0; s_z = 0;
!----------------------------BASE-------------------------------
   call base(s,s_sub)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
   call HAM_INTRA(J2,s,H_intra,N)

   deallocate(s)
!--------------------------------------------------------------
   write(*,*) 'Entre com H e H_max:'
   read*, H,H_max

   write(*,*) 'Entre com a fase:'
   read*, state

   ! open(unit=20, file= 'SO_T_' // trim(state) // "_T-H.dat")

   write(*,*) 'Entre com o step(-3,-5) e o passo de H(-1,-2,-3):'
   read*, cd, ef

         passo = 10.d0**(ef)
   
          if ( state=="AF" ) then
            step = 10.d0**(cd)
         else
            step = -10.d0**(cd)
         end if

         write(*,*) 'Entre com a temperatura MÍN e MÁX:'
         read*, T_min, T_max

   do

      do while(H<=H_max)

         j = 0; Alfa = 0.d0 

        
         if ( state=="AF" ) then
            T_inicial = T_min
            T_final = T_max
         else
            T_inicial = T_max
            T_final = T_min
         end if  


         ! CALL CPU_TIME ( tempo_inicial )

         ! -
         if ( T_inicial>T_final ) then
            step = -10.d0**(cd)
         else
            step = 10.d0**(cd)
         end if


         call mag_vetor(state,m)

         ! - - - - - - - - - - - - - - - - - - - - - - -

         !open(unit=20, file=trim(state) // "_H_T-F-m.dat")

         WRITE (temp, '(F5.3)') H

         open(unit=20, file= trim(temp) // '_SO_T_' // trim(state) // "_T-H.dat")

         do while (T_inicial/=T_final)

            error = 1.d0; erro = 1.d0

            do while(erro >= tol)

               call Ham_inter_state(J2,N,m,H_inter)

               H_total = H_intra + H_inter - H*s_z

               !---------------------- SHIFT DA HAMILTONIANA ----------------
               if (T_inicial<=10.d0**(-1)) then

                  Alfa = minval(H_total)

                  H_total = H_total - Alfa

               endif
               !---------------------- SHIFT DA HAMILTONIANA ----------------

               call partition(H_total,T_inicial,Z)

               do i = 1, 8

                  mag_prev(i) = m(i)

                  call magnetization(H_total,Z,s,i,T_inicial,m(i))

                  error(i) = abs(mag_prev(i) - m(i))

                  ! m1 = m1*0.5 + 0.5*m(1)


                  erro = maxval(error)

               end do


            end do

               do i = 1, 12

                  call magnetization(H_total,Z,s,i,T_inicial,m(i))

               enddo

            call order_parameter(state,m,m_order)

            call F_helm(T_inicial,Z,F)

            F = F + Alfa


            write(20,*) T_inicial, F, m_order

            call date_and_time(date,time,zone,values)
            call date_and_time(DATE=date,ZONE=zone)
            call date_and_time(TIME=time)
            call date_and_time(VALUES=values)

            write(*, '(F8.5,A,F8.5,A,I2,A,I2,A,I2,A,I2)') T_inicial,' ',&
            m_order,' - ', values(5),':',values(6),' do dia ',values(3),'/',values(2)

            if (j==0) then
               if (m_order<=10.d0**(-4)) then
                  print*, '\/---------\/'
                  write(*,18) H, T_inicial
                  print*, '/\---------/\'
18                format ((F8.5))
                  j = 1
                  !write(20,*) T_inicial, H
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




         !----------------- NOVAS ENTRADAS  -----------------


         print*, '------------'
         Print*, 'State =','', state
         print*, 'H =','', H
         print*, 'T','',T_inicial
         print*, 'J2','',J2
         print*, '----END-----'

         !CALL CPU_TIME ( tempo_final )

         ! minutos = int(tempo_final - tempo_inicial)/60
         ! segundos = int(mod((tempo_final - tempo_inicial),60.0))

         ! WRITE (*, '(A,I3,A,I2,A)') 'Demorou ',minutos,' minutos e ',segundos,' segundos'

         ! call system('paplay /usr/share/sounds/gnome/default/alerts/drip.ogg')

         H = H + passo

      enddo

      !call system('paplay /usr/share/sounds/sound-icons/trumpet-12.wav')

      print*, '=== FIM ==='

      close(20)
      stop

   enddo

end program


