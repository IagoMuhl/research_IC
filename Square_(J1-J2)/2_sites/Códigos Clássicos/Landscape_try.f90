program Landscape
   use CMF
   implicit none
   integer, dimension(maxConfig,num_sites):: s
   real*8, dimension(maxConfig):: H_intra, H_inter, H_total,s_x
   real*8, dimension(num_sites):: m
   real*8:: J2, erro, m1, m2, Dif_F,m_order,m_iter
   real*8:: H_inicial,H_final,T,m_fe,m_af,step,Z,tol,erro1,erro2,F(2)
   character(len=3), dimension(2):: states
   character(len=3):: state
   character(len=7):: eme
   integer:: j, up, down, cd, i, X_X, k, l
   !real*8, dimension(:,:),allocatable:: m_order, F

   up = 1; down = 2
   tol = 10.d0**(-8); J2 = -0.1d0


   call base(s)
   call HAM_INTRA(J2,s,H_intra)
!--------------------------------------------------------------
   do i = 1, maxConfig
   s_x(i) = s(i,1) + s(i,2)
   enddo




   ! write(*,*) 'Entre com T:'
   ! read*, T

   ! write(*,*) 'Entre com H_inicial:'
   ! read*, H_inicial

   T = 0.5
   H_inicial = 4.21558
   m_fe = 0.94d0

      do 

         if (m_fe<=0.9) then
            stop
         endif

      j = 1

      states = ['AF','PM']

      state = states(2)


      !m_af = .92655d0

      call mag_vetor(state,0.d0,m_fe,m,m_order)

      call Ham_inter_state(J2,s,m,H_inter)

      !H_inter = 0.d0

      H_total = H_intra + H_inter - H_inicial*s_x

      call partition(H_total,T,Z)

      call F_helm(T,Z,F(2))

      ! - - - - - - - - - - - - - - - - - - - - - - -


         state = states(1)

         step = 10.d0**(-5)


         ! m_fe = .92655d0
         ! m_af = -.92655d0

         m_af=-m_fe

         WRITE (eme, '(F7.5)') m_fe

         open(unit=j, file=  trim(eme) // "_Dif_F-m.dat")

         call mag_vetor(state,m_fe,m_af,m,m_order)

         do while (m_order >= 0.d0)

            m_af = m_fe - 2*m_order

            m = [m_order,m_af]

            call Ham_inter_state(J2,s,m,H_inter)

            H_total = H_intra + H_inter - H_inicial*s_x

            call partition(H_total,T,Z)

         call F_helm(T,Z,F(1))

         Dif_F = F(1) - F(2)

         !call mag_vetor(state,m_fe,m_af,m,m_order)

         write(j,*) m_order,Dif_F,m_af,m_fe,F(1),F(2)

         m_order = m_order - step

      enddo

      close(j)

!

      j = j + 1
      m_fe = m_fe - step

      
   
   enddo
!       print*, '------',State,'------'
!       write(*,19) H_inicial, T
! 19    format ((F8.5))
!  read*,


end program


