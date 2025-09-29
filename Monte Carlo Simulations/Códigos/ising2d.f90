program ising2d
   use mt19937
   implicit none

  integer:: x,y, L, N, i
  integer:: N_steps, N_equilibrio
  integer, allocatable :: s(:,:), seed
  real*8:: r, J_1, E, M, T, T_start, T_end, T_step
  real*8:: E_acum, E2_acum, M_acum, M2_acum, M4_acum, Chi, C, M_abs_acum, U_Binder



  ! ... Parâmetros do sistema
   J_1 = 1.d0
   L = 32
   N = L*L
   T_start = 2.3d0
   T = T_start
   T_end = 2.2d0
   T_step = -0.001d0
   N_equilibrio = 2*10**(6)
   N_steps = 4*10**(6)
   seed = 65465712

   open(L)

   ! ... alocação de memória e inicialização de spins
   allocate(s(L,L))

   call sgrnd(seed)

   ! Loop principal sobre a temperatura
   do while(T>T_end)!T = T_start, T_end, T_step


      ! 1. Inicializa spins aleatoriamente para cada nova T
      do x = 1, L
         do y = 1, L
            r = grnd()
            if (r < 0.5d0) then
               s(x,y) = 1
            else
               s(x,y) = -1
            end if
         end do
      end do

      ! call print_matrix(s,L)
      ! read*,

      ! 2. Fase de Equilíbrio
      do i = 1, N_equilibrio
         call metropolis_step(s, L, J_1, T)
      end do


      ! 3. Reinicia acumuladores
      E_acum = 0.d0
      E2_acum = 0.d0
      M_acum = 0.d0
      M2_acum = 0.d0
      M4_acum = 0.d0
      M = 0.d0
      M_abs_acum = 0.d0
      ! x = 0; y = 0

      ! 4. Fase de Coleta de Dados
      do i = 1, N_steps
         call metropolis_step(s, L, J_1, T)
         call calculate_energy(s, L, E, J_1)
         E_acum = E_acum + E
         E2_acum = E2_acum + E*E

         M = sum(s)
         M_acum = M_acum + abs(M)
         M2_acum = M2_acum + M*M
         M4_acum = M4_acum + M**4

      end do

      ! 5. Calcula as médias e flutuações
      E = E_acum / N_steps
      M = M_acum/ N_steps ! abs para magnetização média

      C = ((E2_acum/N_steps) - (E_acum/N_steps)*(E_acum/N_steps)) / (T*T)
      Chi = (abs(M2_acum/N_steps) - ((M_acum/N_steps)*(M_acum/N_steps))) / T
      ! ... Calcule a capacidade calorífica e a susceptibilidade aqui

      U_Binder = 1.d0 - ((M4_acum/N_steps)/(3*((M2_acum/N_steps)*(M2_acum/N_steps))))

      E = E/N
      M = M/N
      C = C/N
      Chi = Chi/N

      write(L,*) T, E, M, C, Chi, U_Binder

      ! 6. Imprime os resultados para a temperatura T
      write(*,*) 'T = ', T, ' E = ', E, ' M = ', M


      T = T + T_step

   end do ! Fim do loop de temperatura

   close(L)

end program ising2d

subroutine metropolis_step(s, L, J_1, T)
   use mt19937
   implicit none
   integer, intent(in):: L
   real*8, intent(in):: J_1,T
   integer, intent(inout):: s(L,L)
   integer:: site_i, site_j, i1, i2, j1, j2, neigh_sum
   real*8:: r, dE

   r = grnd()
      site_i = int(r*L) + 1
   r = grnd()
      site_j = int(r*L) + 1

      i1 = site_i - 1
   if (i1 < 1) i1 = L
      i2 = site_i + 1
   if (i2 > L) i2 = 1
      j1 = site_j - 1
   if (j1 < 1) j1 = L
      j2 = site_j + 1
   if (j2 > L) j2 = 1

      neigh_sum = s(i1,site_j) + s(i2,site_j) + s(site_i,j1) + s(site_i,j2)

      dE = 2.d0 * J_1 * s(site_i,site_j) * neigh_sum ! dE in {-8,-4,0,4,8} for J=1


   if (dE <= 0) then
      s(site_i,site_j) = -s(site_i,site_j)
   else
      r = grnd()
         if (r < exp(-dE/T)) then
            s(site_i,site_j) = -s(site_i,site_j)   
         endif
   endif

end subroutine

subroutine calculate_energy(s, L, E, J_1)
  implicit none
  integer, intent(in) :: s(L,L)
  integer, intent(in) :: L
  real*8, intent(out) :: E
  real*8, intent(in) :: J_1
  integer:: i, j, i1, j1
  
  E = 0.d0
  !M = 0.d0

  ! Percorre toda a rede para calcular M e E
  do i = 1, L
     do j = 1, L
        ! Magnetização é a soma dos spins
        !M = M + s(i,j)

        ! Para calcular a energia, some as interações com os vizinhos
        ! Para evitar contagem dupla, some apenas com o vizinho à direita e para baixo
        ! Lembre-se das condições de contorno periódicas
        
        i1 = i + 1
        if (i1 > L) i1 = 1
        j1 = j + 1
        if (j1 > L) j1 = 1

        ! Interação horizontal
        E = E + s(i,j) * s(i1,j)
        
        ! Interação vertical
        E = E + s(i,j) * s(i,j1)
        
     end do
  end do
  
  ! Energia é E = -J * sum(si*sj)
  E = -J_1 * E

end subroutine

subroutine print_matrix(A,dim)
   implicit none
   integer, intent(in) :: dim
   integer, dimension(dim,dim), intent(in) :: A
   integer :: i, j

   do i = 1, dim
      write(*,20) (A(i,j), j=1,dim)
   20    format(16(I6))   ! até 16 inteiros por linha, cada um ocupando 6 colunas
   end do
end subroutine

subroutine write_spin_matrix(s, L, T, y)
   implicit none
   integer, intent(in) :: s(L,L), L, y
   real*8, intent(in):: T
   character(len=5):: filename, file_Y
   integer :: i, j, visual_unit

   visual_unit = 21
   write (file_Y, '(I5)') y
   write (filename, '(F5.3)') T
   open(unit=visual_unit, file= 'T='//trim(filename)//'_'//trim(file_Y)//'.dat')
   do i = 1, L
      write(visual_unit,*) (s(i,j), j=1,L)
   end do
   close(visual_unit)
end subroutine
