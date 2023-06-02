program quant_J
   use QUANTICO
   implicit none

   real*8, parameter:: J3 = 0.15d0
   real*8, dimension(2,2):: sigma_z, Id, sigma_x
   real*8, dimension(2*2,2*2):: sig_zz, Id_2, Id_sig_z, sig_z_Id,Id_sigma_x, sigma_x_Id
   real*8, dimension(2**4,2**4):: F , H_1, H_2, H_intra, Id_4, H_inter, Ham, H_gama
   real*8, dimension(2**4,2**4):: s1, s2, s3, s4 ,s1_x, s2_x, s3_x,s4_x,s_x,V!,H,H_prime
   real*8 :: Z, T, Gamma_final, step, m, tol, erro, m_guess, J2, F_helm, F_prime
   real*8, dimension(16):: W
   real*8:: Gamma,Alfa
   character(len=3):: state
   character(len=5) :: nameFileJ2, nameFileJ3
   integer:: dim!,i,j



   H_1 = 0; H_2 = 0; W = 0; V = 0; dim = 2;

   tol = 10.d0**(-8); T = 10.d0**(-5)
!---------------------------------------------------------
! CALCULO DAS POSSIBILIDADES DE SIGMA-Z E IDENTIDADE

   sigma_z = reshape([1,0,0,-1],[dim,dim])

   sigma_x = reshape([0,1,1,0],[dim,dim])

   Id = reshape([1,0,0,1],[dim,dim])

   call tensorial(sigma_z,sigma_z,dim,sig_zz)
   call tensorial(Id,Id,dim,Id_2)
   call tensorial(sigma_z,Id,dim,sig_z_Id)
   call tensorial(Id,sigma_z,dim,Id_sig_z)


   call tensorial(sigma_x,Id,dim,sigma_x_Id)
   call tensorial(Id,sigma_x,dim,Id_sigma_x)

   call tensorial(sig_z_Id,Id_2,dim*dim,s1)
   call tensorial(Id_sig_z,Id_2,dim*dim,s2)
   call tensorial(Id_2,sig_z_Id,dim*dim,s3)
   call tensorial(Id_2,Id_sig_z,dim*dim,s4)

    call tensorial(sigma_x_Id,Id_2,dim*dim,s1_x)
    call tensorial(Id_sigma_x,Id_2,dim*dim,s2_x)
    call tensorial(Id_2,sigma_x_Id,dim*dim,s3_x)
    call tensorial(Id_2,Id_sigma_x,dim*dim,s4_x)

   call tensorial(Id_2,Id_2,dim*dim,Id_4)

   s_x = s1_x + s2_x + s3_x + s4_x


   !call print_matrix(s_x,dim**4)
   
   
!---------------- HAMILTONIANA J1-------------------------

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


!---------------------------------------------------------
   !call print_matrix(s1,dim**4)

   do

      J2 = 0.0
   
      print*, 'Entre com J2'
         read(*,*) J2
         if ( J2==-1 ) stop 'Fim da rotina'
            
         print*, 'Entre com a fase (AF,SAF,SD,PM)'
         read(*,*)   state

!---------------------------------------------------------
!DECLARAÇÃO DE VALORES INICIAIS

         Gamma = 1.1d0; 

         Gamma_final = 1.8d0;

         step = 10.d0**(-5); 

         m_guess = 1.d0;
      
!---------------------------------------------------------

         H_intra = J1*H_1 + J2*H_2

!------------------------ OPEN UNIT -----------------
   WRITE (nameFileJ2, '(F5.2)') j2
   WRITE (nameFileJ3, '(F5.2)') j3

   open(unit=20, file=trim(state) // "_J3_gamma-F-m.dat")
   !open(unit=20, file=trim(state) // "_J3(" // trim(adjustl(nameFileJ3)) // ")_gamma-F-m.dat")
   !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
!----------------------------------------------------

   do while (Gamma <= Gamma_final) !FUNÇÃO DE PARTIÇÃO/ LOOP TEMPERATURA

      ERRO = 1.d0
      
      H_gama = (-1)*Gamma*s_x


      do while (ERRO >= tol)

         call Ham_inter_state(state,J2,J3,s1,s2,s3,s4,m_guess,Id_4,H_inter)
 

         Ham = H_intra + H_inter + H_gama

         ! call print_matrix(H_inter,dim**4,dim**4)
         ! read(*,*)

         call diagonalization(Ham,V,W)



         !---------------------- SHIFT DA HAMILTONIANA ----------------

         Alfa = abs(W(1))

         W = W + Alfa

         !---------------------- SHIFT DA HAMILTONIANA ----------------


         call partition(W,T,dim,Z)

         call magnetization_diag(W,Z,T,dim,s1,V,m)

         ERRO = abs(m_guess - m)

         ! print*, Gamma, m_guess, m, erro
         ! read(*,*)

         m_guess = m

      end do



      call Free_nrg(T,Z,F_helm)

      F_prime = (F_helm - Alfa)

         write(*,*) Gamma

         write(20,*) Gamma, F_prime, m

      Gamma = Gamma + step

   end do
   
   print*, '------------'
   write(*,*) State, J2

   close(20)

   end do


end program
