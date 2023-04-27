program quant_J
   use QUANTICO
   implicit none

   real*8, parameter:: J3 = 0.45d0
   real*8, dimension(2,2):: sigma_z, Id, sigma_x
   real*8, dimension(2*2,2*2):: sig_zz, Id_2, Id_sig_z, sig_z_Id ! sigma_x_Id, Id_sigma_x
   real*8, dimension(2**4,2**4):: F , H_1, H_2, H_intra, Id_4, H_inter, Ham
   real*8, dimension(2**4,2**4):: s1, s2, s3, s4 !s1_x, s2_x, s3_x, s4_x, s_x, D_x, V_x
   real*8 :: Z, T, T_final, step, m, tol, erro, m_guess, J2, F_helm
   character(len=3):: state
   character(len=5) :: nameFileJ2, nameFileJ3
   integer:: dim


   H_1 = 0; H_2 = 0; dim = 2;

   step = 10.d0**(-5); tol = 10.d0**(-10)
!---------------------------------------------------------
! CALCULO DAS POSSIBILIDADES DE SIGMA-Z E IDENTIDADE

   sigma_z = reshape([1,0,0,-1],[dim,dim])

   sigma_x = reshape([0,1,1,0],[dim,dim])

   Id = reshape([1,0,0,1],[dim,dim])

   call tensorial(sigma_z,sigma_z,dim,sig_zz)
   call tensorial(Id,Id,dim,Id_2)
   call tensorial(sigma_z,Id,dim,sig_z_Id)
   call tensorial(Id,sigma_z,dim,Id_sig_z)

   ! call tensorial(sigma_x,Id,dim,sigma_x_Id)
   ! call tensorial(Id,sigma_x,dim,Id_sigma_x)

   call tensorial(sig_z_Id,Id_2,dim*dim,s1)
   call tensorial(Id_sig_z,Id_2,dim*dim,s2)
   call tensorial(Id_2,sig_z_Id,dim*dim,s3)
   call tensorial(Id_2,Id_sig_z,dim*dim,s4)

   ! call tensorial(sigma_x_Id,Id_2,dim*dim,s1_x)
   ! call tensorial(Id_sigma_x,Id_2,dim*dim,s2_x)
   ! call tensorial(Id_2,sigma_x_Id,dim*dim,s3_x)
   ! call tensorial(Id_2,Id_sigma_x,dim*dim,s4_x)

   call tensorial(Id_2,Id_2,dim*dim,Id_4)
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



   do

      J2 = 0.0
   
      print*, 'Entre com J2'
         read(*,*) J2
         if ( J2==-1 ) stop 'Fim da rotina'
            
         print*, 'Entre com a fase (AF,SAF,SD,PM)'
         read(*,*)   state

!---------------------------------------------------------
!DECLARAÇÃO DE VALORES INICIAIS

         T = 0.1d0; T_final = 2.d0;
      
         m_guess = 1.d0;

!-----------------------------------------

   H_intra = J1*H_1 + J2*H_2


   WRITE (nameFileJ2, '(F5.2)') j2
   WRITE (nameFileJ3, '(F5.2)') j3

   open(unit=20, file=trim(state) // "_T-F_J3(" // trim(adjustl(nameFileJ3)) // ").dat")
   !open(unit=20, file=trim(state) // "_T-F_J2(" // trim(adjustl(nameFileJ2)) // ")_J3(" // trim(adjustl(nameFileJ3)) // ").dat")

   do while (T<=T_final) !FUNÇÃO DE PARTIÇÃO/ LOOP TEMPERATURA

      call Ham_inter_state(state,J2,J3,s1,s2,s3,s4,m_guess,Id_4,H_inter)

      Ham = H_intra + H_inter

      call partition(Ham,T,dim,Z)

      call magnetization(Ham,Z,T,dim,s1,m)


      ERRO = abs(m_guess - m)

      do while (ERRO >= tol)

         m_guess = m

         call Ham_inter_state(state,J2,J3,s1,s2,s3,s4,m_guess,Id_4,H_inter)

         Ham = H_intra + H_inter

         call partition(Ham,T,dim,Z)

         call magnetization(Ham,Z,T,dim,s1,m)

         ERRO = abs(m_guess - m)

      end do


      !  if (m<=10.d0**(-5)) then
      !     print*, T, m
      !     exit
      !  end if

      call Free_nrg(T,Z,F_helm)

      write(20,*) T, F_helm, m

      T = T + step

   end do

      close(20)

   end do

end program
