program square_H
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_x
   real*8, dimension(num_sites):: m
   real*8:: J2, erro, m1, m2
   real*8:: H_inicial,H_final,T,m_fe,m_af,step,Z,m_order,tol,F,erro1,erro2
   character(len=3):: state
   integer:: j, up, down, cd, i

   up = 1; down = 2
   tol = 10.d0**(-8); J2 = -0.1d0

   call base(s)
!--------------------------------------------------------------
   do i = 1, maxConfig
   s_x(i) = s(i,1) + s(i,2)
   enddo

   call HAM_INTRA(J2,s,H_intra)

   do

      j = 0


      write(*,*) 'Entre com T, step:'
      read*, T, cd

      write(*,*) 'Entre com H_inicial:'
      read*, H_inicial

      write(*,*) 'Entre com H_final:'
      read*, H_final

      ! H_inicial = 3.8d0
      ! H_final = 4.25d0

      write(*,*) 'Entre com a fase:'
      read*, state

      ! -
      if ( H_inicial>H_final ) then
         step = -10.d0**(cd)
      else
         step = 10.d0**(cd)
      end if
!-

      if(state/='2AF') then
         m_fe = 1.0d0;
         m_af = -1.0d0;
      else
         m_fe = 0.84719110987579493
         m_af = 0.65197042353076307
         !   m_fe = 0.99099828895786968
         !   m_af =  0.44969820018918100
         ! m_fe = 1.0d0;
         ! m_af = -1.0d0;
      endif

      ! - - - - - - - - - - - - - - - - - - - - - - -


 

      ! - - - - - - - - - - - - - - - - - - - - - - -

      call mag_vetor(state,m_fe,m_af,m,m_order)



      ! - - - - - - - - - - - - - - - - - - - - - - -

      open(unit=20, file=trim(state) // "_T_H-F-m.dat")


      do while (H_inicial/=H_final)

         erro1 = 1.d0; erro2 = 1.d0; erro = 1.d0

         do while(erro >= tol)

            ! print*, m1, m2, erro
            ! read(*,*)

            call Ham_inter_state(J2,s,m,H_inter)

            H_total = H_intra + H_inter -H_inicial*s_x

            call partition(H_total,T,Z)

            call magnetization(H_total,Z,s,up,T,m1)

            call magnetization(H_total,Z,s,down,T,m2)


            erro1 = abs(m(1) - m1)
            erro2 = abs(m(2) - m2)

            ! if(state=='PM') then
            !    erro2 = abs(m_fe - m(down))
            ! endif

            ! m1 = 0.5*m1 + 0.5*m(1)
            ! m2 = 0.5*m2 + 0.5*m(2)


            erro = max(abs(erro1),abs(erro2))

            !  print*, m_fe, m_af, m_order
            !  read(*,*)

            call mag_vetor(state,m1,m2,m,m_order)


         end do

         call F_helm(T,Z,F)

         !print*, H_inicial, m_order


         write(20,*) H_inicial, F, m_order, m1, m2
         ! print*, T, m_order, m_fe, m_af

         if (j==0) then
            if (m_order<=10.d0**(-4)) then
               print*, '\/---------\/'
               write(*,18) H_inicial, T, m1, m2
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

      print*, '------',State,'------'
      write(*,19) H_inicial, T
19    format ((F8.5))

   enddo

end program


