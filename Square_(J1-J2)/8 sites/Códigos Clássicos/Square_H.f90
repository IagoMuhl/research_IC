program square_H
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_z
   real*8, dimension(8):: m, error, mag_prev
   real*8:: J2, erro
   real*8:: H_inicial,H_final,T,step,Z,m_order,tol,F
   character(len=3):: state
   integer:: j, cd, i,p


   tol = 10.d0**(-8); J2 = -0.345d0; s_z = 0
!---------------------------------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo

!--------------------------------------------
   do

      write(*,*) 'Entre com T, step:' ; j = 0
      read*, T, cd

      write(*,*) 'Entre com H_inicial:'
      read*, H_inicial

      write(*,*) 'Entre com H_final:'
      read*, H_final

      ! H_inicial = 1.d0
      ! H_final = 3.55d0

      write(*,*) 'Entre com a fase:'
      read*, state

      ! -
      if ( H_inicial>H_final ) then
         step = -10.d0**(cd)
      else
         step = 10.d0**(cd)
      end if
!-

      ! if(state/='2AF') then
      !    m1 = 1.d0
      !    m2 = -1.d0
      !    m3 = 1.d0
      !    m4 = -1.d0
      ! else
      !    m1 = 0.99982421570227475
      !    m2 = 1.4944181483173785E-002
      !    m3 = 0.99500816175228535
      !    m4 = 3.2426469508950429E-003
      !    ! m_fe = 0.84719110987579493
      !    ! m_af = 0.65197042353076307
      ! endif





      ! - - - - - - - - - - - - - - - - - - - - - - -

      call mag_vetor(state,m)

      call HAM_INTRA(J2,s,H_intra)

      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_T_H-F-m.dat")


      do while (H_inicial/=H_final)

         error = 1.d0; erro = 1.d0

         do while(erro >= tol)

            ! print*, m1, m2, erro
            ! read(*,*)

            call Ham_inter_state(J2,s,m,H_inter)


            H_total = H_intra + H_inter - H_inicial*s_z


            call partition(H_total,T,Z)

            do i = 1, num_sites

               mag_prev(i) = m(i)

               call magnetization(H_total,Z,s,i,T,m(i))

               error(i) = abs(mag_prev(i) - m(i))

               ! m1 = m1*0.5 + 0.5*m(1)

               erro = maxval(error)

            end do

            call order_parameter(state,m,m_order)


         end do

         call F_helm(T,Z,F)

         !print*, H_inicial, m_order


         write(20,*) H_inicial, F, m_order,m(1),m(2)!,m(3),m(4),m(5),m(6),m(7),m(8)
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

   enddo

end program


