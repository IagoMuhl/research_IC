program quant_TGammafix
   use QUANTICO
   implicit none

   integer, parameter:: L = 2
   real*8, dimension(2,2):: H_1, H_2, H_intra, H_inter_1,H_inter_2, H_gamma, H_long
   real*8, dimension(2,2)::  V,V2!s_x
   real*8 :: Z, T, step, tol, erro1,erro2, m_fe, m_af, J2, F_helm, F_prime, m_order, m2, m1
   real*8, dimension(2):: W,W2
   real*8, dimension(2):: m
   real*8, dimension(:,:), allocatable:: sigma_x, sigma_z, Id
   !real*8, dimension(:,:), allocatable:: s1_x, s2_x, s3_x, s4_x, F
   real*8:: H_inicial,H_final, Alfa, Alfa2, Gamma, erro,Z2,F
   character(len=3):: state
   character(len=5) :: nameFileJ2
   integer:: dim,i,cd!,j

   H_1 = 0; H_2 = 0; W = 0; V = 0; dim = 2;

   tol = 10.d0**(-8); J2 = -0.1d0 ;  
   !---------------------------------------------------------
! CALCULO DAS POSSIBILIDADES DE SIGMA-Z E IDENTIDADE

   allocate( sigma_x(l,l),sigma_z(l,l),Id(l,l))

   sigma_z = reshape([1,0,0,-1],[dim,dim])

   sigma_x = reshape([0,1,1,0],[dim,dim])

   Id = reshape([1,0,0,1],[dim,dim])

   !call print_matrix(sigma_x,2,2)

   !allocate(sig_zz(4:4,4:4),Id_2(4:4,4:4),Id_sig_z(4:4,4:4),sig_z_Id(4:4,4:4),Id_sigma_x(4:4,4:4), sigma_x_Id(4:4,4:4))

   ! call tensorial(sigma_z,sigma_z,dim,sig_zz)
   !  call tensorial(Id,Id,dim,Id_2)

   ! call tensorial(sigma_z,Id,dim,s1)
   ! call tensorial(Id,sigma_z,dim,s2)

   ! s_z = s1 + s2
   ! call tensorial(sigma_x,Id,dim,sigma_x_Id)
   ! call tensorial(Id,sigma_x,dim,Id_sigma_x)

   ! deallocate (sigma_x, sigma_z, Id)

   ! call tensorial(sig_z_Id,Id_2,dim*dim,s1)
   ! call tensorial(Id_sig_z,Id_2,dim*dim,s2)


   !  allocate(s1_x(l**4,l**4),s2_x(l**4,l**4))

   ! call tensorial(sigma_x_Id,Id_2,dim*dim,s1_x)
   ! call tensorial(Id_sigma_x,Id_2,dim*dim,s2_x)
   ! call tensorial(Id_2,sigma_x_Id,dim*dim,s3_x)
   ! call tensorial(Id_2,Id_sigma_x,dim*dim,s4_x)

   ! call tensorial(Id_2,Id_2,dim*dim,Id_4)

   ! s_x = s1_x + s2_x + s3_x + s4_x

   ! deallocate(Id_sigma_x, sigma_x_Id, s1_x, s2_x, s3_x,s4_x)


   !deallocate (sigma_x,sigma_z,Id,sig_zz,Id_2,Id_sig_z,sig_z_Id,Id_sigma_x,sigma_x_Id,s1_x, s2_x, s3_x,s4_x)

   !call print_matrix(s_x,dim**4)

   !---------------- HAMILTONIANA J1-------------------------
   ! allocate(F(4,4))
   ! !S1*S2
   ! call tensorial(sig_zz,Id_2,dim,F)

   ! H_1 = H_1 + F


   ! !S1*S3
   ! call tensorial(sig_z_Id,sig_z_Id,dim*dim,F)

   ! H_1 = H_1 + F


   ! !S2*S4
   ! call tensorial(Id_sig_z,Id_sig_z,dim*dim,F)

   ! H_1 = H_1 + F

   ! !S3*S4
   ! call tensorial(Id_2,sig_zz,dim*dim,F)

   ! H_1 = H_1 + F


   !---------------- HAMILTONIANA J2-------------------------


   ! call tensorial(sig_z_Id,Id_sig_z,dim*dim,F)

   ! H_2 = H_2 + F

   ! call tensorial(Id_sig_z,sig_z_Id,dim*dim,F)

   ! H_2 = H_2 + F

   ! deallocate(sig_zz, Id_2,Id_sig_z, sig_z_Id, F)


!---------------------------------------------------------



   do

       !T = 10.d0**(-5)
      Gamma = 0.d0
      i = 0
      Alfa = 0.d0 
      Alfa2 = 0.d0 
      !print*, 'Entre com T'
      !read(*,*) T
      !if ( T==-1 ) stop 'Fim da rotina'

      ! print*, 'Entre com J2, Step(-5,-3):'
      ! read(*,*) J2, cd

      print*, 'Entre com T, Step(-5,-3):'
      read(*,*) T, cd

      !H = 10.d0**(-5)

      ! print*, 'Entre com Gamma, Step(-5,-3):'
      ! read(*,*) Gamma, cd

      print*, 'Entre com H_inicial'
      read(*,*) H_inicial

      print*, 'Entre com H_final'
      read(*,*) H_final

      ! H_inicial = 3.9d0
      ! H_final = 4.8d0

      print*, 'Entre com a fase (AF,SAF,SD,PM)'
      read(*,*)   state

         ! -
   if ( H_inicial>H_final ) then
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
      m_fe = 0.99099828895786968!1.0d0;
      m_af =  0.44969820018918100 !-1.0d0;
   endif


      ! print_H = H_inicial

      ! call print_matrix(s1,4,4)
      ! print*, '------------------'
      ! read(*,*)
      ! call print_matrix(s2,4,4)
      ! read(*,*)
      !---------------------------------------------------------

      H_intra = 0.d0

      H_gamma = (-1)*Gamma*sigma_x



      !------------------------ OPEN UNIT -----------------
      WRITE (nameFileJ2, '(F5.2)') j2
      !WRITE (nameFileJ3, '(F5.2)') j3

      open(unit=20, file=trim(state) // "_T_gamma_H-F-m.dat")
      !open(unit=20, file=trim(state) // "_J3(" // trim(adjustl(nameFileJ3)) // ")_gamma-F-m.dat")
      !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
      !----------------------------------------------------

      call chose(state,m_fe,m_af,m)

      

      do while (H_inicial/=H_final) !FUNÇÃO DE PARTIÇÃO/ LOOP TEMPERATURA

         erro1 = 1.d0; erro2 = 1.d0; erro = 1.d0

         H_long = (-1)*H_inicial*sigma_z

         do while(erro >= tol)

            call Ham_inter_state(J2,sigma_z,[m(1),m(2)],Id,H_inter_1)

            call Ham_inter_state(J2,sigma_z,[m(2),m(1)],Id,H_inter_2)
            ! print*, m, erro1
            !  read(*,*)

            !!Ham = H_intra + H_inter_1 + H_gamma + H_long !+ H_inter_2

            ! print*, m
            ! read(*,*)

            !  call print_matrix(H_inter_1,dim,dim)
            !  read(*,*)

            !  call print_matrix(H_inter_2,dim,dim)
            !  read(*,*)

            H_inter_1 = H_inter_1 + H_long + H_gamma
            H_inter_2 = H_inter_2 + H_long + H_gamma


            call diagonalization(H_inter_1,V,W)

            call diagonalization(H_inter_2,V2,W2)
            

            !---------------------- SHIFT DA HAMILTONIANA ----------------

            if (T<=10.d0**(-3)) then

               Alfa = abs(W(1))

               W = W + Alfa

               Alfa2 = abs(W2(1))

               W2 = W2 + Alfa2

            endif

            !---------------------- SHIFT DA HAMILTONIANA ----------------

            ! print*, W, W2
            ! read(*,*)

            call partition(W,T,dim,Z)
            call partition(W2,T,dim,Z2)
            !   print*, m, erro
            !   read(*,*)

            call magnetization_diag(W,Z,T,dim,sigma_z,V,m1)
            call magnetization_diag(W2,Z2,T,dim,sigma_z,V2,m2)
            ! call magnetization_diag(W,Z,T,dim,s3,V,m3)
            ! call magnetization_diag(W,Z,T,dim,s4,V,m4)




            erro1 = abs((m1)) - abs(m(1))
            erro2 = abs((m2)) - abs(m(2))
            ! erro3 = abs((m3)) - abs(m(3))
            ! erro4 = abs((m4)) - abs(m(4))
 
            ! erro = max(abs(erro1),abs(erro2))

            erro = max(abs(erro1),abs(erro2))
            ! erro = abs(erro1)
            ! print*, m
            ! read(*,*)
 
            call mag_vetor(state,m1,m2,m,m_order)


         end do


         call Free_nrg(T,Z,F_helm)
         call Free_nrg(T,Z2,F_prime)

         !F_prime = (F_helm - Alfa)
         F = (F_helm + F_prime)/2.d0 - (Alfa + Alfa2)/2.d0
         !write(*,*) H_inicial, m_order

         write(20,*) H_inicial, F, m_order, m1, m2


            if (i==0) then
               if (m_order<=10.d0**(-4)) then
                  print*, '\/------------\/'
                  write(*,18) H_inicial, T
                  print*, '/\------------/\'
   18             format ((F8.5))
                  i = 1
               end if
            end if


         H_inicial = H_inicial + step

         if ((max(H_inicial,H_final))-(min(H_inicial,H_final))<=abs(step)) then
            H_inicial = H_final
         endif

      end do

      print*, '------------'
      Print*, 'State =','', state 
      print*, 'H =','', H_inicial
      print*, 'Gamma','',Gamma
      print*, 'T','',T
      print*, 'J2','',J2
      print*, '----END-----'

      close(20)

   end do


end program