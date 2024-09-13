program Landscape
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_x
   real*8, dimension(num_sites):: m
   real*8:: J2, erro, m1, m2, Dif_F
   real*8:: H_inicial,H_final,T,m_fe,m_af,step,Z,tol,erro1,erro2
   character(len=3), dimension(2):: states
   character(len=3):: state
   integer:: j, up, down, cd, i, X_X, k, l
   real*8, dimension(:,:),allocatable:: m_order, F

   up = 1; down = 2
   tol = 10.d0**(-8); J2 = -0.1d0; T = 1.32d0


   call base(s)
   call HAM_INTRA(J2,s,H_intra)
!--------------------------------------------------------------
   do i = 1, maxConfig
   s_x(i) = s(i,1) + s(i,2)
   enddo




      ! write(*,*) 'Entre com T, step:'
      ! read*, T, cd

      ! write(*,*) 'Entre com H_inicial:'
      ! read*, H_inicial

      ! write(*,*) 'Entre com H_final:'
      ! read*, H_final



      ! do i = 1,2

      ! write(*,*) 'Entre com a fase ',i
      ! read*, states(i) 

      ! enddo

      states = ['AF','PM']

      do 
         j = 0

         open(unit=20, file="Land_" // trim(states(1)) // "-Dif_F-m.dat")
      ! - - - - - - - - - - - - - - - - - - - - - - -
         do i = 1,2

            state = states(i)

            H_inicial = 3.95d0
            H_final = 4.05d0

      ! -
      if ( state == 'AF' ) then
         step = 10.d0**(-5)
      else
         step = -10.d0**(-5)
      end if

      X_X = int(abs(H_inicial-H_final)/step)

      allocate(F(2,X_X), m_order(2,X_X))
!-

      if(states(1)/='2AF') then
         m_fe = 1.0e-1;
         m_af = -1.0e-1;
      else
         m_fe = 0.84719110987579493
         m_af = 0.65197042353076307
      endif

      ! - - - - - - - - - - - - - - - - - - - - - - -

      call mag_vetor(state,m_fe,m_af,m,m_order(i,1))


      ! - - - - - - - - - - - - - - - - - - - - - - -


      do while (H_inicial/=H_final)

      k = 1

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

            call mag_vetor(state,m1,m2,m,m_order(i,k))


         end do

         call F_helm(T,Z,F(i,k))

         !write(20,*) H_inicial, F, m_order, m1, m2
         ! print*, T, m_order, m_fe, m_af

!          if (j==0) then
!             if (m_order<=10.d0**(-4)) then
!                print*, '\/---------\/'
!                write(*,18) H_inicial, T
!                print*, '/\---------/\'
! 18             format ((F8.5))
!                j = 1
!             end if
!          end if


   end do

         H_inicial = H_inicial + step
         k = k + 1

         if ((max(H_inicial,H_final))-(min(H_inicial,H_final))<=abs(step)) then
            H_inicial = H_final
         endif

         ! if (max(H_final,H_inicial)>=12.d0) then
         !    exit
         ! endif

      end do

      do l = 1,X_X

      Dif_F = F(1,l) - F(2,X_X-(l-1))

      enddo

   write(20,*) Dif_F, m_order(1,X_X)

      close(20)

      print*, '------',State,'------'
      write(*,19) H_inicial, T
19    format ((F8.5))
 read*,
   enddo

end program


