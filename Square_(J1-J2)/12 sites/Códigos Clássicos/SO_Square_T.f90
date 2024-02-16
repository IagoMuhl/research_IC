program square_T
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_z
   real*8, dimension(6):: m, error, mag_prev
   real*8:: J2, erro, Alfa
   real*8:: T_inicial,T_final,H,step,Z,m_order,tol,F
   real*4:: tempo_inicial, tempo_final
   character(len=3):: state
   integer:: j, cd,i,p, minutos, segundos


   tol = 10.d0**(-8); J2 = -0.35d0; s_z = 0;
!----------------------------BASE-------------------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
   call HAM_INTRA(J2,s,H_intra)
!--------------------------------------------------------------
   write(*,*) 'Entre com H:'
   read*, H

   write(*,*) 'Entre com a fase:'
   read*, state

   open(unit=20, file= 'SO_T_' // trim(state) // "_T-H.dat")

   do

      do while(H>=0.5d0**(-2))

      j = 0; Alfa = 0.d0 ; cd = -3

      T_inicial = 3.5d0
      T_final = 4.8d0

      
      CALL CPU_TIME ( tempo_inicial )

      ! -
      if ( T_inicial>T_final ) then
         step = -10.d0**(cd)
      else
         step = 10.d0**(cd)
      end if


      call mag_vetor(state,m)

      ! - - - - - - - - - - - - - - - - - - - - - - -

      !open(unit=20, file=trim(state) // "_H_T-F-m.dat")


      do while (T_inicial/=T_final)

         error = 1.d0; erro = 1.d0

         do while(erro >= tol)

            call Ham_inter_state(J2,s,m,H_inter)

            H_total = H_intra + H_inter - H*s_z

         !---------------------- SHIFT DA HAMILTONIANA ----------------
            if (T_inicial<=10.d0**(-3)) then

               Alfa = minval(H_total)

               H_total = H_total - Alfa

            endif
         !---------------------- SHIFT DA HAMILTONIANA ----------------

            call partition(H_total,T_inicial,Z)

            do i = 1, 6

               mag_prev(i) = m(i)

               call magnetization(H_total,Z,s,i,T_inicial,m(i))

               error(i) = abs(mag_prev(i) - m(i))

               ! m1 = m1*0.5 + 0.5*m(1)


               erro = maxval(error)

            end do

            call order_parameter(state,m,m_order)


         end do

         call F_helm(T_inicial,Z,F)

         F = F + Alfa

         !print*, T_inicial, m_order


         !write(20,*) T_inicial, F, m_order,m(1),m(2)!,m(3),m(4),m(5),m(6),m(7),m(8)
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '\/---------\/'
               write(*,18) H, T_inicial
               print*, '/\---------/\'
18             format ((F8.5))
               j = 1
               write(20,*) T_inicial, H
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

      CALL CPU_TIME ( tempo_final )

      minutos = int(tempo_final - tempo_inicial)/60
      segundos = int(mod((tempo_final - tempo_inicial),60.0))

      WRITE (*, '(A,I2,A,I2,A)') 'Demorou ',minutos,' minutos e ',segundos,' segundos'

      call system('paplay /usr/share/sounds/gnome/default/alerts/drip.ogg')

      H = H - 0.15

   enddo

   call system('paplay /usr/share/sounds/sound-icons/trumpet-12.wav')

   print*, '=== FIM ==='

   close(20)
   stop

   enddo

end program


