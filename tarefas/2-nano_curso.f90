program nano2
   use QUANTICO
   implicit none

   real*8, parameter:: J2 = 0.d0
   real*8, dimension(2,2):: sigma_z, Id
   real*8, dimension(2*2,2*2):: sig_zz, Id_2, sig_Idz, sig_zId
   real*8, dimension(2**4,2**4):: F , H_1, H_2, H_intra, Id_4, H_inter, Ham
   real*8, dimension(2**4,2**4):: s1, s2, s3, s4
   real*8:: Z, T, T_final, step, mag, tol, erro
   real*8, dimension(4):: m_guess
   character(len=3):: state
   integer:: m

   ! write(*,*) 'Qual estado desejas calcular?(AF,PM,SAF)'
   ! read(*,*) state
   state = 'AF'
!---------------------------------------------------------
!DECLARAÇÃO DE VALORES

   H_1 = 0; H_2 = 0; T = 0.5d0; T_final = 4.d0; 
   step = 10.d0**(-3); tol = 10.d0**(-8)
   mag = 1.d0; m = 2



   sigma_z = reshape([1,0,0,-1],[m,m])

   Id = reshape([1,0,0,1],[m,m])

!---------------------------------------------------------
! CALCULO DAS POSSIBILIDADES DE SIGMA-Z E IDENTIDADE

 call tensorial(sigma_z,sigma_z,m,sig_zz)
 call tensorial(Id,Id,m,Id_2)
 call tensorial(sigma_z,Id,m,sig_zId)
 call tensorial(Id,sigma_z,m,sig_Idz)

 call tensorial(sig_zId,Id_2,m*m,s1)
 call tensorial(sig_Idz,Id_2,m*m,s2)
 call tensorial(Id_2,sig_zId,m*m,s3)
 call tensorial(Id_2,sig_Idz,m*m,s4)

 call tensorial(Id_2,Id_2,m*m,Id_4)
!---------------- HAMILTONIANA J1-------------------------

                   
 call tensorial(sig_zz,Id_2,m*m,F)

   H_1 = H_1 + F



   call tensorial(sig_zId,sig_zId,m*m,F)

   H_1 = H_1 + F

 

   call tensorial(sig_Idz,sig_Idz,m*m,F)

   H_1 = H_1 + F

   call tensorial(Id_2,sig_zz,m*m,F)

   H_1 = H_1 + F


!---------------- HAMILTONIANA J2-------------------------
   

   call tensorial(sig_zId,sig_Idz,m*m,F)

   H_2 = H_2 + F

   call tensorial(sig_Idz,sig_zId,m*m,F)

   H_2 = H_2 + F


!---------------------------------------------------------

   H_intra = J1*H_1 + J2*H_2

   

do while (T<=T_final) !FUNÇÃO DE PARTIÇÃO/ LOOP TEMPERATURA
   
   call mag_vetor(state,mag,m_guess)
   
   call Ham_inter_state(state,J2,s1,s2,s3,s4,m_guess,Id_4,H_inter)

   Ham = H_intra + H_inter
   
   call partition(Ham,T,m,Z)
   
   call magnetization(Ham,Z,T,m,s1,mag)
   
   ERRO = abs(mag - m_guess(1))

   !print*, T, ERRO, mag, m_guess(1)
   call print_matrix(s1,m**4)
   call print_matrix(Ham,m**4)
   read(*,*)
      do while (ERRO >= tol)

         !print*, T, mag, m_guess(1), erro
         

         call mag_vetor(state,mag,m_guess)

         call Ham_inter_state(state,J2, s1,s2,s3,s4,m_guess,Id_4,H_inter)

         Ham = H_intra + H_inter

         call partition(Ham,T,m,Z)


         call magnetization(H_inter,Z,T,m,s1,mag)

         ERRO = abs(mag - m_guess(1))

         print*, T, mag, m_guess(1), erro
         read(*,*)
      end do

   T = T + step

   !write(10,*) T, mag

   !print *, 'T =',T,'Z =',Z, 'MAG =',abs(m_guess(1))

end do


!print *, 'T=',T,'Z=',Z, 'MAG=',mag

end program
