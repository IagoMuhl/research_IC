program quant_TGammafix
   use QUANTICO
   implicit none

   integer, parameter:: L = 2
   real*8, dimension(2**4,2**4):: H_1, H_2, H_intra, Id_4, H_inter, Ham, H_gama, H_long
   real*8, dimension(2**4,2**4)::  s1, s2, s3, s4 ,s_x,V,s_z
   real*8 :: Z, T, step, tol, erro1,erro2, m_fe, m_af, J2, F_helm, F_prime, m_order
   real*8, dimension(16):: W
   real*8, dimension(4):: m
   real*8, dimension(:,:), allocatable:: sigma_x, sigma_z, Id, sig_zz, Id_2,Id_sig_z, sig_z_Id,Id_sigma_x, sigma_x_Id
   real*8, dimension(:,:), allocatable:: s1_x, s2_x, s3_x, s4_x, F
   real*8:: H_inicial,H_final, Alfa, Gamma, print_H
   character(len=3):: state
   character(len=5) :: nameFileJ2
   integer:: dim,i,up,down,cd!,j

   H_1 = 0; H_2 = 0; W = 0; V = 0; dim = 2;

   tol = 10.d0**(-8); J2 = 0.d0 ;
   !---------------------------------------------------------
! CALCULO DAS POSSIBILIDADES DE SIGMA-Z E IDENTIDADE

   allocate( sigma_x(l,l),sigma_z(l,l),Id(l,l),sig_zz(l**2,l**2) &
      ,Id_2(l**2,l**2),Id_sig_z(l**2,l**2),sig_z_Id(l**2,l**2),Id_sigma_x(l**2,l**2),sigma_x_Id(l**2,l**2))

   sigma_z = reshape([1,0,0,-1],[dim,dim])

   sigma_x = reshape([0,1,1,0],[dim,dim])

   Id = reshape([1,0,0,1],[dim,dim])

   !call print_matrix(sigma_x,2,2)

   !allocate(sig_zz(4:4,4:4),Id_2(4:4,4:4),Id_sig_z(4:4,4:4),sig_z_Id(4:4,4:4),Id_sigma_x(4:4,4:4), sigma_x_Id(4:4,4:4))

   call tensorial(sigma_z,sigma_z,dim,sig_zz)
   call tensorial(Id,Id,dim,Id_2)
   call tensorial(sigma_z,Id,dim,sig_z_Id)
   call tensorial(Id,sigma_z,dim,Id_sig_z)


   call tensorial(sigma_x,Id,dim,sigma_x_Id)
   call tensorial(Id,sigma_x,dim,Id_sigma_x)

   deallocate (sigma_x, sigma_z, Id)

   call tensorial(sig_z_Id,Id_2,dim*dim,s1)
   call tensorial(Id_sig_z,Id_2,dim*dim,s2)
   call tensorial(Id_2,sig_z_Id,dim*dim,s3)
   call tensorial(Id_2,Id_sig_z,dim*dim,s4)

   allocate(s1_x(l**4,l**4),s2_x(l**4,l**4),s3_x(l**4,l**4),s4_x(l**4,l**4))

   call tensorial(sigma_x_Id,Id_2,dim*dim,s1_x)
   call tensorial(Id_sigma_x,Id_2,dim*dim,s2_x)
   call tensorial(Id_2,sigma_x_Id,dim*dim,s3_x)
   call tensorial(Id_2,Id_sigma_x,dim*dim,s4_x)

   call tensorial(Id_2,Id_2,dim*dim,Id_4)

   s_x = s1_x + s2_x + s3_x + s4_x

   deallocate(Id_sigma_x, sigma_x_Id, s1_x, s2_x, s3_x,s4_x)

   s_z = s1 + s2 + s3 + s4

   !deallocate (sigma_x,sigma_z,Id,sig_zz,Id_2,Id_sig_z,sig_z_Id,Id_sigma_x,sigma_x_Id,s1_x, s2_x, s3_x,s4_x)

   !call print_matrix(s_x,dim**4)

   !---------------- HAMILTONIANA J1-------------------------
   allocate(F(l**4,l**4))
   !S1*S2
   call tensorial(sig_zz,Id_2,dim*dim,F)

   H_1 = H_1 + F


   !S1*S3
   call tensorial(sig_z_Id,sig_z_Id,dim*dim,F)

   H_1 = H_1 + F


   !S2*S4
   call tensorial(Id_sig_z,Id_sig_z,dim*dim,F)

   H_1 = H_1 + F

   !S3*S4
   call tensorial(Id_2,sig_zz,dim*dim,F)

   H_1 = H_1 + F


   !---------------- HAMILTONIANA J2-------------------------


   call tensorial(sig_z_Id,Id_sig_z,dim*dim,F)

   H_2 = H_2 + F

   call tensorial(Id_sig_z,sig_z_Id,dim*dim,F)

   H_2 = H_2 + F

   deallocate(sig_zz, Id_2,Id_sig_z, sig_z_Id, F)


!---------------------------------------------------------



   do

      !T = 10.d0**(-5)
      Gamma = 10.d0**(-5)!3.33d0
      i = 0

      !print*, 'Entre com T'
      !read(*,*) T
      !if ( T==-1 ) stop 'Fim da rotina'

      print*, 'Entre com T, Step(-5,-3):'
      read(*,*) T, cd

      ! print*, 'Entre com Gamma, Step(-5,-3):'
      ! read(*,*) Gamma, cd

      ! print*, 'Entre com H_inicial'
      ! read(*,*) H_inicial

      ! print*, 'Entre com H_final'
      ! read(*,*) H_final

      H_inicial = 2.5d0
      H_final = 5.d0

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

      m_fe = 1.0d0;
      m_af = -1.d0

      print_H = H_inicial

      !---------------------------------------------------------

      H_intra = J1*H_1 + J2*H_2

      H_gama = (-1)*Gamma*s_x

      !------------------------ OPEN UNIT -----------------
      WRITE (nameFileJ2, '(F5.2)') j2
      !WRITE (nameFileJ3, '(F5.2)') j3

      open(unit=20, file=trim(state) // "_T_gamma_H-F-m.dat")
      !open(unit=20, file=trim(state) // "_J3(" // trim(adjustl(nameFileJ3)) // ")_gamma-F-m.dat")
      !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
      !----------------------------------------------------

      do while (H_inicial/=H_final) !FUNÇÃO DE PARTIÇÃO/ LOOP TEMPERATURA

         erro1 = 1.d0; erro2 = 1.d0

         H_long = (-1)*H_inicial*s_z

         do while(max(erro1,erro2) >= tol)

            call Ham_inter_state(J2,s1,s2,s3,s4,m,Id_4,H_inter)

            Ham = H_intra + H_inter + H_gama + H_long

            ! call print_matrix(H_inter,dim**4,dim**4)
            ! read(*,*)

            call diagonalization(Ham,V,W)

            !---------------------- SHIFT DA HAMILTONIANA ----------------

            Alfa = abs(W(1))

            W = W + Alfa

            !---------------------- SHIFT DA HAMILTONIANA ----------------

            call partition(W,T,dim,Z)

            !  print*, m
            !  read(*,*)

            call magnetization_diag(W,Z,T,dim,s1,V,m_fe)
            call magnetization_diag(W,Z,T,dim,s2,V,m_af)

            erro1 = abs(m_fe - m(up))
            erro2 = abs(m_af - m(down))

            if(state=='PM') then
               erro2 = abs(m_fe - m(down))
            endif

            ! print*, Gamma_inicial, m_fe, m, erro
            ! read(*,*)
            call mag_vetor(state,m_fe,m_af,up,down,m,m_order)

         end do


         call Free_nrg(T,Z,F_helm)

         F_prime = (F_helm - Alfa)

         !write(*,*) H_inicial, m_order

         write(20,*) H_inicial, F_prime, m_order

         if (H_inicial>=1.1*print_H) then
            if (i==0) then
               if (m_order<=10.d0**(-4)) then
                  print*, '\/------------\/'
                  write(*,18) H_inicial, T
                  print*, '/\------------/\'
   18             format ((F8.5))
                  i = 1
               end if
            end if
           endif

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
      print*, '----END-----'

      close(20)

   end do


end program

