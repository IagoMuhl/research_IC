program square_H
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_z
   real*8:: m(8), error(6), mag_prev(6)
   real*8:: J2, erro, Alfa
   real*8:: H_inicial,H_final,T,step,Z,m_order,tol,F
   real*4:: tempo_inicial, tempo_final
   character(len=3):: state
   integer:: j, cd, i, p, minutos, segundos



   tol = 10.d0**(-8); J2 = 0.0d0; s_z = 0
!--------------------BASE---------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
   call HAM_INTRA(J2,s,H_intra)
!--------------------------------------------
   do
      j = 0; Alfa = 0.d0
      write(*,*) 'Entre com T, step:' 
      read*, T, cd

      write(*,*) 'Entre com H_inicial:'
      read*, H_inicial

      write(*,*) 'Entre com H_final:'
      read*, H_final

      ! H_inicial = 2.7d0
      ! H_final = 4.4d0

      write(*,*) 'Entre com a fase:'
      read*, state

      CALL CPU_TIME ( tempo_inicial )

      ! -
      if ( H_inicial>H_final ) then
         step = -10.d0**(cd)
      else
         step = 10.d0**(cd)
      end if

      call mag_vetor(state,m)

      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_T_H-F-m.dat")


      do while (H_inicial/=H_final)

         error = 1.d0; erro = 1.d0

         do while(erro >= tol)

            ! print*, m1, m2, erro
            ! read(*,*)

            call Ham_inter_state(J2,s,m,H_inter)


            H_total = H_intra + H_inter - H_inicial*s_z

         !---------------------- SHIFT DA HAMILTONIANA ----------------
            if (T<=10.d0**(-3)) then

               Alfa = minval(H_total)

               H_total = H_total - Alfa

            endif
         !---------------------- SHIFT DA HAMILTONIANA ----------------

            call partition(H_total,T,Z)

            ! print*, m(1)


            do i = 1, 6

               mag_prev(i) = m(i)

               call magnetization(H_total,Z,s,i,T,m(i))

               error(i) = abs(mag_prev(i) - m(i))

               ! m1 = m1*0.5 + 0.5*m(1)


               erro = maxval(error)

            end do

         end do

            do i = 7, 8

            call magnetization(H_total,Z,s,(i+10),T,m(i))

            enddo

            call order_parameter(state,m,m_order)

         call F_helm(T,Z,F)

         F = F + Alfa

         !print*, H_inicial, m_order

         write(20,*) H_inicial, F, m_order!,m(1),m(2)!,m(3),m(4),m(5),m(6),m(7),m(8)
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '\/---------\/'
               write(*,18) H_inicial, T
               print*, '/\---------/\'
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

      print*, '------------'
      Print*, 'State =','', state
      print*, 'H =','', H_inicial
      print*, 'T','',T
      print*, 'J2','',J2
      print*, '----END-----'

      CALL CPU_TIME ( tempo_final )

      minutos = int(tempo_final - tempo_inicial)/60
      segundos = int(mod((tempo_final - tempo_inicial),60.0))

      WRITE (*, '(A,I2,A,I2,A)') 'Demorou ',minutos,' minutos e ',segundos,' segundos'

      call system('paplay /usr/share/sounds/gnome/default/alerts/drip.ogg')

   enddo

end program


