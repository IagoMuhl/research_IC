program ising2d_blocos
   use mt19937
   implicit none

  integer:: x, y, L, N, i, j, k
  integer:: N_steps, N_equilibrio, N_blocks, block_size
  integer, allocatable :: s(:,:), seed
  real*8:: r, J_1, E, M, T, T_start, T_end, T_step
  real*8:: E_total, M_total, M_abs_total, E2_total, M2_total
  real*8:: E_block_sum, M_abs_block_sum, E2_block_sum, M2_block_sum
  real*8:: Chi_total, C_total, M_abs
  integer:: L_unit
  
  real*8, allocatable :: E_block_avg(:), M_abs_block_avg(:), &
                      C_block_avg(:), Chi_block_avg(:)
  real*8:: E_err, M_abs_err, C_err, Chi_err
  
  ! ... Parâmetros do sistema
   J_1 = 1.d0
   L = 20
   N = L*L
   T_start = 2.5d0
   T = T_start
   T_end = 2.d0
   T_step = -0.01d0
   N_equilibrio = 2000000
   N_steps = 4000000
   N_blocks = 200
   block_size = N_steps/N_blocks
   seed = 43574357

   L_unit = 10
   open(unit=L_unit, file='ising.dat', status='replace')

   ! Aloca as matrizes para as médias dos blocos
   allocate(s(L,L))
   allocate(E_block_avg(N_blocks), M_abs_block_avg(N_blocks), &
            C_block_avg(N_blocks), Chi_block_avg(N_blocks))

   ! Inicializa o gerador de numeros aleatorios com a semente
   call sgrnd(seed)

   ! Loop principal sobre a temperatura
   do while(T>T_end)

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

      ! 2. Fase de Equilíbrio
      do i = 1, N_equilibrio
         call metropolis_step(s, L, J_1, T)
      end do

      ! 3. Fase de Coleta de Dados com amostragem por blocos
      do i = 1, N_blocks
         ! Zera os acumuladores para este bloco
         E_block_sum = 0.d0
         E2_block_sum = 0.d0
         M_abs_block_sum = 0.d0
         M2_block_sum = 0.d0

         do j = 1, block_size
            ! Realiza um passo de Monte Carlo
            call metropolis_step(s, L, J_1, T)
            
            ! Calcula as quantidades fisicas no estado atual
            call calculate_energy(s, L, E, J_1)
            M = sum(s)
            
            ! Acumula os valores para o bloco
            E_block_sum = E_block_sum + E
            E2_block_sum = E2_block_sum + E*E
            M_abs_block_sum = M_abs_block_sum + abs(M)
            M2_block_sum = M2_block_sum + M*M
         end do
         
         ! Calcula e armazena a média de cada bloco
         E_block_avg(i) = E_block_sum / real(block_size)
         M_abs_block_avg(i) = M_abs_block_sum / real(block_size)
         C_block_avg(i) = (E2_block_sum/real(block_size) - E_block_avg(i)**2) / T**2
         Chi_block_avg(i) = (M2_block_sum/real(block_size) - M_abs_block_avg(i)**2) / T
      end do
      
      ! 4. Calcula as médias globais e os erros a partir das médias dos blocos
      call calculate_mean_and_error(E_block_avg, N_blocks, E_total, E_err)
      call calculate_mean_and_error(M_abs_block_avg, N_blocks, M_abs_total, M_abs_err)
      call calculate_mean_and_error(C_block_avg, N_blocks, C_total, C_err)
      call calculate_mean_and_error(Chi_block_avg, N_blocks, Chi_total, Chi_err)

      ! Normaliza por spin
      E_total = E_total/N
      M_abs_total = M_abs_total/N
      C_total = C_total/N
      Chi_total = Chi_total/N

      ! 5. Escreve o resultado no arquivo, incluindo os erros
      write(L_unit,*) T, E_total, E_err, M_abs_total, M_abs_err, C_total, Chi_total
      
      ! 6. Imprime os resultados para a tela
      write(*,*) 'T = ', T, ' E = ', E_total, ' M = ', M_abs_total

      T = T + T_step

   end do ! Fim do loop de temperatura

   close(L_unit)
   deallocate(s, E_block_avg, M_abs_block_avg, C_block_avg, Chi_block_avg)

end program ising2d_blocos

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

   dE = 2.d0 * J_1 * s(site_i,site_j) * neigh_sum

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
  
  do i = 1, L
     do j = 1, L
        i1 = i + 1
        if (i1 > L) i1 = 1
        j1 = j + 1
        if (j1 > L) j1 = 1
        E = E + s(i,j) * s(i1,j)
        E = E + s(i,j) * s(i,j1)
     end do
  end do
  
  E = -J_1 * E

end subroutine

! Esta subrotina calcula a media e o erro estatistico de um array de medias de blocos
subroutine calculate_mean_and_error(block_averages, N_blocks, mean, error)
  implicit none
  integer, intent(in) :: N_blocks
  real*8, dimension(N_blocks), intent(in) :: block_averages
  real*8, intent(out) :: mean, error
  real*8:: squared_diff_sum
  integer:: i

  ! Calcula a media global
  mean = sum(block_averages) / real(N_blocks)

  ! Calcula a soma dos quadrados das diferencas para o desvio padrao
  squared_diff_sum = 0.d0
  do i = 1, N_blocks
    squared_diff_sum = squared_diff_sum + (block_averages(i) - mean)**2
  end do
  
  ! Calcula o erro (desvio padrao da media)
  if (N_blocks > 1) then
      error = sqrt(squared_diff_sum / real(N_blocks * (N_blocks - 1)))
  else
      error = 0.d0
  endif
end subroutine