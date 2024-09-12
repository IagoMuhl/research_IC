! Landscape: Landscape of the free energy as function of the order parameter
! Descrição: construímos landscapes da energia livre em função do parâmetro
!de ordem (magnetizaão). Não ha autoconsistência neste caso.
! Fixamos: Temperatura, J2 e Jp, e então registramos !num arquivo de dados a diferença
!F(SAF) - F(PM) para o espectro |m|<0.5.
!
!Matheus Roos, 26/01/2024
!Last update:05/09/2023.

PROGRAM Landscape
   use CMF
   implicit none
   character(len=18), parameter :: path = './Plots/Landscape/'
   character(len=4), parameter :: ext = '.dat'
   character(len=5), parameter :: cols = '_f-m'
   integer, parameter :: nStates = 2
   integer, parameter :: loTc = 1, Tc = 2, hiTc = 3 !menor que Tc, igual a Tc e maior que Tc

   integer, dimension(maxConfig, numSitios) :: S
   real(kind=db), dimension(maxConfig) :: H_intra, H_inter, H
   real(db), dimension(3) :: T
   real(kind=db), dimension(nStates) :: f

   character(len=7) :: StrTemp
   character(len=5) :: StrJ2
   character(len=6) :: StrJp
   character(len=3) :: state
   character(len=3), dimension(nStates) :: states
   character(len=12), dimension(3) :: Ttype
   real(kind=db) :: mag_min, mag_max, mag_iter, mag, delta_f
   real(kind=db) :: J2, Jp, Z, step
   integer :: iT


   !Seleção do estado.
   states = ['SAF', 'PM ']

   call setup()

   !!! GERAÇÃO DA BASE INDEPENDENTE DA TEMPERATURA !!!
   !Retorna uma matriz com todas as configurações possíveis.
   call baseSpins(S, maxConfig, numSitios)

   varJ2: do
      write(*,*) 'Entre com J2: [0 stop]'
      read(*,*) J2
      if (J2 == 0) exit

      write(StrJ2, '(F5.2)') J2

      call Hamilton_intra(J2, JP, S, H_intra)

      write(*,*) 'Entre com <Tc, =Tc, >Tc:'
      read(*,*) T(loTc), T(Tc), T(hiTc)

      varT: do iT = 1, 3
         write(StrTemp,'(F7.4)') T(iT)

         ! Criamos e escrevemos o cabeçalho do dataFile.
         open(unit=20, file=trim(path) // &
         & "Jp(" // trim(adjustl(StrJp)) // ")/" // "J2_" // trim(adjustl(StrJ2)) // &
         & trim(cols) // "_T_" // trim(adjustl(StrTemp)) // ext)

         write(20,*) '# Bilayer square lattice: J1-J2-Jp'
         write(20,*) '# Landscape da energia livre: T ' // trim(Ttype(iT))
         write(20,10) Jp, J2, T(iT)

         ! Começamos com a fase mais fácil, PM.
         ! Não há necessidade de um loop sob a mag., pois ela é nula.
         ! Também não precisamos calcular a hamiltoniana interclsuter, já que todos termos são
         !proporcionais à mag., que é zero.
         ! Assim, basta calcular a função de partição (Z) e a energia livre [f(2)]
         mag = 0.D0
         state = states(2)
         call Hamilton_inter(J2, S, mag, state, H_inter)

         H = H_intra + H_inter

         call Zfunction(T(iT), H, Z)

         f(2) = -T(iT)*dlog(Z)

         ! Resetamos o contador:
         mag_iter = 0
         varMag: do
            state = states(1)

            mag = mag_min + mag_iter*step

            if (mag > mag_max) exit

            call Hamilton_inter(J2, S, mag, state, H_inter)

            H = H_intra + H_inter

            call Zfunction(T(iT), H, Z)

            !freeEnergy
            f(1) = -T(iT)*dlog(Z)

            ! Calculamos a diferença entre a energia livre da fases ordenada e PM.
            delta_f = f(1) - f(2)

            write(20,*) mag, delta_f

            mag_iter = mag_iter + 1
         end do varMag

         close(20)
      end do varT
   end do varJ2

10 format (' # Jp:', F4.1, ', J2:', F6.3, ' ,T:', F8.5, / )
contains
   subroutine setup()
      implicit none
      character(len=1) :: precision
      integer :: n

      !header
      print *, '=========================='
      write(*,*) 'Fases selecionadas: ' // trim(states(1)) // "-" // trim(states(2))
      print *, '=========================='

      !Escolha dos parâmetros.
      write(*,*) 'Escolha o parâmetro Jp:'
      read(*,*) Jp
      if ( Jp <= 1 ) then
         write(StrJp, '(F4.1)') Jp
      else
         write(StrJp, '(F6.1)') Jp
      end if

      write(*,*) 'Escolha o parâmetro +/- mag?'
      read(*,*) mag_max
      mag_min = -mag_max

      !Step com que varia J2.
      write(*,*) 'Qual a precisão dos pontos low/hight? [l/h]'
      read(*,*) precision
      precision = trim(upperCase(precision))
      if ( precision == 'L' ) then
         n = 3
      else
         n = 5
      end if
      step = 10.d0**(-n)

      !temperatura qualis
      Ttype(1) = 'menor que Tc'
      Ttype(2) = 'maior que Tc'
      Ttype(3) = 'igual a Tc'

      call system('mkdir -p ' // trim(path) // 'Jp\(' // trim(adjustl(StrJp)) // '\)/')
   end subroutine
END PROGRAM Landscape
