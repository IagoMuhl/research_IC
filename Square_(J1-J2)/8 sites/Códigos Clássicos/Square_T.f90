program square_T
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_z
   real*8, dimension(8):: m, error, mag_prev
   real*8:: J2, erro
   real*8:: T_inicial,T_final,H,step,Z,m_order,tol,F
   character(len=3):: state
   integer:: j, cd,i,p


   tol = 10.d0**(-8); J2 = -0.345d0; s_z = 0
!--------------------------------------------------------------
   call base(s)

   do p = 1, num_sites
      do i = 1, maxConfig
         s_z(i) = s_z(i) + s(i,p)
      enddo
   enddo
!--------------------------------------------------------------
   do

      write(*,*) 'Entre com H, step:' ; j = 0
      read*, H, cd

      write(*,*) 'Entre com T_inicial:'
      read*, T_inicial

      write(*,*) 'Entre com T_final:'
      read*, T_final

      ! T_inicial = 2.4d0
      ! T_final = 3.8d0

      write(*,*) 'Entre com a fase:'
      read*, state

      ! -
      if ( T_inicial>T_final ) then
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

      open(unit=20, file=trim(state) // "_H_T-F-m.dat")


      do while (T_inicial/=T_final)

         error = 1.d0; erro = 1.d0

         do while(erro >= tol)

            call Ham_inter_state(J2,s,m,H_inter)

            H_total = H_intra + H_inter - H*s_z

            call partition(H_total,T_inicial,Z)

            do i = 1, num_sites

            mag_prev(i) = m(i)

            call magnetization(H_total,Z,s,i,T_inicial,m(i))

            ! call magnetization(H_total,Z,s,2,T_inicial,m2)

            ! call magnetization(H_total,Z,s,3,T_inicial,m3)

            ! call magnetization(H_total,Z,s,4,T_inicial,m4)


            error(i) = abs(mag_prev(i) - m(i))
            ! erro2 = abs(m(2) - m2)
            ! erro3 = abs(m(3) - m3)
            ! erro4 = abs(m(4) - m4)

            ! m1 = m1*0.5 + 0.5*m(1)
            ! m2 = m2*0.5 + 0.5*m(2)
            ! m3 = m3*0.5 + 0.5*m(3)
            ! m4 = m4*0.5 + 0.5*m(4)


            erro = maxval(error)

            end do 

            call order_parameter(state,m,m_order)


         end do

         call F_helm(T_inicial,Z,F)

         !print*, T_inicial, m_order


         write(20,*) T_inicial, F, m_order,m(1),m(2)!,m(3),m(4),m(5),m(6),m(7),m(8)
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '\/---------\/'
               write(*,18) H, T_inicial
               print*, '/\---------/\'
18             format ((F8.5))
               j = 1
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

      close(20)


      !----------------- NOVAS ENTRADAS  -----------------


      print*, '------------'
      Print*, 'State =','', state 
      print*, 'H =','', H
      print*, 'T','',T_inicial
      print*, 'J2','',J2
      print*, '----END-----'



   enddo

end program


