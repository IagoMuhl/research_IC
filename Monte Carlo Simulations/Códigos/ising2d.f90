program ising2d
  implicit none
  integer:: x,y, L, N, i
  integer, allocatable :: s(:,:)
  real*8:: r

   L = 2
   N = L*L

   allocate(s(L,L))

   do x = 1, L
     do y = 1, L
        call random_number(r)
        if (r < 0.5d0) then
           s(x,y) = 1
        else
           s(x,y) = -1
        end if
     end do
  end do

  call print_matrix(s,L)

 do i = 1, N
       call random_number(r)
       site_i = int(r*x) + 1
       call random_number(r)
       site_j = int(r*y) + 1

       ! periodic neighbors
       i1 = site_i - 1
       if (i1 < 1) i1 = nx
       i2 = site_i + 1
       if (i2 > nx) i2 = 1
       j1 = site_j - 1
       if (j1 < 1) j1 = ny
       j2 = site_j + 1
       if (j2 > ny) j2 = 1

       neigh_sum = s(i1,site_j) + s(i2,site_j) + s(site_i,j1) + s(site_i,j2)
       dE = 2 * s(site_i,site_j) * neigh_sum    ! delta E = 2 * s_i * sum_neighbors * J (J multiplicative below)

       ! convert to actual energy difference *J
       dE = int(dE * 1)  ! dE in {-8,-4,0,4,8} for J=1

       if (dE <= 0) then
          ! aceita
          s(site_i,site_j) = -s(site_i,site_j)
          E = E + real(dE, dp) * J
          M = M + 2.0_dp * real(s(site_i,site_j), dp) * (-1.0_dp) ! update magnetization
          acceptance_count = acceptance_count + 1.0_dp
       else
          ! probabilidade exp(-beta * dE)
          idx = dE/4 + 3
          call random_number(rnd)
          if (rnd .lt. expd(idx)) then
             s(site_i,site_j) = -s(site_i,site_j)
             E = E + real(dE, dp) * J
             M = M + 2.0_dp * real(s(site_i,site_j), dp) * (-1.0_dp)
             acceptance_count = acceptance_count + 1.0_dp
          end if
       end if
    end do






end program ising2d

subroutine print_matrix(A,dim)
      implicit none
      integer, intent(in):: dim
      integer, dimension(dim,dim), intent(in):: A
      integer:: i, j, n
      n = dim
      do i=1,dim
         write(*,20)(A(i,j),j=1,dim) !para numeros inteiros (4(I3))
      20   format(16(I3))
      enddo
end subroutine