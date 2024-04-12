program quant_THfix
   use QUANTICO
   implicit none

   integer, parameter:: L = 2
   real*8, dimension(4,4)::  H_1,H_2, H_intra, H_inter, H_gamma, H_long
   real*8, dimension(4,4)::  V!s_x
   real*8 :: Z, H, step, tol, erro1,erro2, m_fe, m_af, J2, F_helm, m_order, m2, m1
   real*8, dimension(4):: W
   real*8, dimension(2):: m
   real*8, dimension(:,:), allocatable:: sigma_x, sigma_z, Id, Id_2
   real*8, dimension(:,:), allocatable:: s1_x, s2_x, s1, s2, s_x, s_z
   real*8:: Gamma_inicial,Gamma_final, Alfa, Alfa2, T, erro,F
   character(len=3):: state
   character(len=5) :: nameFileJ2
   integer:: dim,i,cd!,j

   H_1 = 0; H_2 = 0; W = 0; V = 0; dim = 2;

   tol = 10.d0**(-8); J2 = 0.d0 ;  
   !---------------------------------------------------------
! CALCULO DAS POSSIBILIDADES DE SIGMA-Z E IDENTIDADE

   allocate( sigma_x(l,l),sigma_z(l,l),Id(l,l))

   sigma_z = reshape([1,0,0,-1],[dim,dim])

   sigma_x = reshape([0,1,1,0],[dim,dim])

   Id = reshape([1,0,0,1],[dim,dim])

   !call print_matrix(sigma_x,2,2)

   allocate(s1(4,4),s2(4,4),s_z(4,4),s1_x(4,4),s2_x(4,4),s_x(4,4),Id_2(4,4))

   ! call tensorial(sigma_z,sigma_z,dim,sig_zz)
   call tensorial(Id,Id,dim,Id_2)

   call tensorial(sigma_z,Id,dim,s1)
   call tensorial(Id,sigma_z,dim,s2)

   call tensorial(sigma_x,Id,dim,s1_x)
   call tensorial(Id,sigma_x,dim,s2_x)

    s_z = s1 + s2
    s_x = s1_x + s2_x


   !---------------- HAMILTONIANA J1-------------------------

   !S1*S2
   call tensorial(sigma_z,sigma_z,dim,H_1)

   H_intra = J1*H_1

   ! call print_matrix(H_1,4,4)

   deallocate (sigma_x, sigma_z, Id)
!---------------------------------------------------------



   do

       T = 10.d0**(-5)
      ! Gamma = 10.d0**(-5)!3.33d0
      ! H = 0.d0
      i = 0
      Alfa = 0.d0 
      Alfa2 = 0.d0 
      !print*, 'Entre com T'
      !read(*,*) T
      !if ( T==-1 ) stop 'Fim da rotina'

      ! print*, 'Entre com J2, Step(-5,-3):'
      ! read(*,*) J2, cd

      print*, 'Entre com H, Step(-5,-3):'
      read(*,*) H, cd

      ! print*, 'Entre com T, Step(-5,-3):'
      ! read(*,*) T, cd

      !H = 10.d0**(-5)

      ! print*, 'Entre com Gamma, Step(-5,-3):'
      ! read(*,*) Gamma, cd

      ! print*, 'Entre com Gamma_inicial'
      ! read(*,*) Gamma_inicial

      ! print*, 'Entre com Gamma_final'
      ! read(*,*) Gamma_final

      Gamma_inicial = 2.d0
      Gamma_final = 5.d0

      print*, 'Entre com a fase (AF,SAF,SD,PM)'
      read(*,*)   state

         ! -
   if ( Gamma_inicial>Gamma_final ) then
      step = -10.d0**(cd)
   else
      step = 10.d0**(cd)
   end if
!-

      !---------------------------------------------------------
      !DECLARAÇÃO DE VALORES INICIAIS

   if(state/='2AF') then
      m_fe = 1.0d0;
      m_af = -1.0d0;
     else
        m_fe = 0.99099828895786968 !Para transições AF-->AF
        m_af =  0.44969820018918100 
     endif



      !------------------------ OPEN UNIT -----------------
      WRITE (nameFileJ2, '(F5.2)') j2
      !WRITE (nameFileJ3, '(F5.2)') j3

      open(unit=20, file=trim(state) // "_H_T_gamma-F-m.dat")
      !open(unit=20, file=trim(state) // "_J3(" // trim(adjustl(nameFileJ3)) // ")_gamma-F-m.dat")
      !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
      !----------------------------------------------------

      call chose(state,m_fe,m_af,m)

      H_long = (-1.d0)*H*s_z

      do while (Gamma_inicial/=Gamma_final) !FUNÇÃO DE PARTIÇÃO/ LOOP TEMPERATURA

         erro1 = 1.d0; erro2 = 1.d0; erro = 1.d0

         H_gamma = (-1)*Gamma_inicial*s_x

         do while(erro >= tol)

            call Ham_inter_state(J2,s1,s2,m,Id_2,H_inter)


            ! print*, m, erro1
            !  read(*,*)

            !!Ham = H_intra + H_inter + H_gama + H_long !+ H_inter_2

            !  call print_matrix(H_inter,dim,dim)
            !  read(*,*)

            !  call print_matrix(H_inter_2,dim,dim)
            !  read(*,*)

            H_inter = H_inter + H_long + H_gamma + H_intra



            call diagonalization(H_inter,V,W)

            

            !---------------------- SHIFT DA HAMILTONIANA ----------------

            if (T<=10.d0**(-3)) then

               Alfa = abs(W(1))

               W = W + Alfa


            endif

            !---------------------- SHIFT DA HAMILTONIANA ----------------

            call partition(W,T,dim*2,Z)

            !   print*, m, erro
            !   read(*,*)

            call magnetization_diag(W,Z,T,dim*2,s1,V,m1)
            call magnetization_diag(W,Z,T,dim*2,s2,V,m2)


            erro1 = abs((m1)) - abs(m(1))
            erro2 = abs((m2)) - abs(m(2))

 


            erro = max(abs(erro1),abs(erro2))
            ! erro = abs(erro1)
            ! print*, m
            ! read(*,*)
 
            call mag_vetor(state,m1,m2,m,m_order)


         end do



         call Free_nrg(Gamma_inicial,Z,F_helm)


         F = (F_helm - Alfa)

         !write(*,*) T_inicial, m_order

         write(20,*) Gamma_inicial, F, m_order, m1, m2


            if (i==0) then
               if (m_order<=10.d0**(-4)) then
                  print*, '\/------------\/'
                  write(*,18) Gamma_inicial, H
                  print*, '/\------------/\'
   18             format ((F8.5))
                  i = 1
               end if
            end if


         Gamma_inicial = Gamma_inicial + step

         if ((max(Gamma_inicial,Gamma_final))-(min(Gamma_inicial,Gamma_final))<=abs(step)) then
            Gamma_inicial = Gamma_final
         endif

      end do

      print*, '------------'
      Print*, 'State =','', state 
      print*, 'H =','', H
      print*, 'Gamma','',Gamma_inicial
      print*, 'T','',T
      print*, 'J2','',J2
      print*, '----END-----'

      close(20)

   end do


end program