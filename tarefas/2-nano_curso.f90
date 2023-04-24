program nano2
   use QUANTICO
   implicit none

   real*8, parameter:: J2 = 0.5d0
   real*8, dimension(2,2):: sigma_z, Id
   real*8, dimension(2*2,2*2):: sig_zz, Id_2, Id_sig_z, sig_z_Id
   real*8, dimension(2**4,2**4):: F , H_1, H_2, H_intra, Id_4, H_inter, Ham
   real*8, dimension(2**4,2**4):: s1, s2, s3, s4
   real*8:: Z, T, T_final, step, m, tol, erro, m_guess
   character(len=3):: state
   integer:: dim

   ! write(*,*) 'Qual estado desejas calcular?(AF,PM,SAF)'
   ! read(*,*) state
   state = 'AF'
!---------------------------------------------------------
!DECLARAÇÃO DE VALORES

   H_1 = 0; H_2 = 0;
   
   T = 0.5d0; T_final = 4.d0; 

   step = 10.d0**(-3); tol = 10.d0**(-8)

   m_guess = 1.d0; dim = 2; 



   sigma_z = reshape([1,0,0,-1],[dim,dim])

   Id = reshape([1,0,0,1],[dim,dim])

!---------------------------------------------------------
! CALCULO DAS POSSIBILIDADES DE SIGMA-Z E IDENTIDADE

 call tensorial(sigma_z,sigma_z,dim,sig_zz)
 call tensorial(Id,Id,dim,Id_2)
 call tensorial(sigma_z,Id,dim,sig_z_Id)
 call tensorial(Id,sigma_z,dim,Id_sig_z)

 call tensorial(sig_z_Id,Id_2,dim*dim,s1)
 call tensorial(Id_sig_z,Id_2,dim*dim,s2)
 call tensorial(Id_2,sig_z_Id,dim*dim,s3)
 call tensorial(Id_2,Id_sig_z,dim*dim,s4)

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

   H_intra = J1*H_1 + J2*H_2
   

do while (T<=T_final) !FUNÇÃO DE PARTIÇÃO/ LOOP TEMPERATURA

   call Ham_inter_state(state,J2,s1,s2,s3,s4,m_guess,Id_4,H_inter)

   Ham = H_intra + H_inter
   
   call partition(Ham,T,dim,Z)

   call magnetization(Ham,Z,T,dim,s1,m)


   ERRO = abs(m_guess - m)

      do while (ERRO >= tol)
               
         m_guess = m

         call Ham_inter_state(state,J2,s1,s2,s3,s4,m_guess,Id_4,H_inter)
         
         Ham = H_intra + H_inter
         
         call partition(Ham,T,dim,Z)
        
         call magnetization(Ham,Z,T,dim,s1,m)
         !print*, m
         !read(*,*)
         ERRO = abs(m_guess - m)
 
      end do

      !print*, T,m

      if (m<=10.d0**(-5)) then
         print*, T, m
         exit
      end if
   T = T + step

   write(10,*) T, m

end do


!print *, 'T=',T,'Z=',Z, 'm=',m
end program
