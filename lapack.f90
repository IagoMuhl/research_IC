program lapack
    implicit none
    integer , parameter :: db = 8, dim = 2
    real ( KIND = db ) , dimension ( dim , dim ) :: A
    ! !!! Para diagonalizar !!!!!
    character ( len =1) :: JOBZ , UPLO
    integer :: N , LDA , INFO , lwork
    integer , parameter :: LWMAX = 1000
    real ( KIND = db ) , dimension ( dim ) :: W
    real ( kind = db ) , dimension ( LWMAX ) :: WORK
    real ( KIND = db ) , dimension ( dim , dim ) :: sigma_z, sigma_x, Id 
    real ( KIND = db ) , dimension ( dim*2 , dim*2 ) :: sig_zz, sig_z_Id, Id_2, Id_sig_z 
    real ( KIND = db ) , dimension ( dim*2 , dim*2 ) :: sigma_x_Id, Id_sigma_x
    real(kind=db), dimension(dim**4, dim**4):: s1_x,s2_x,s3_x, s4_x
    
    sigma_z = reshape([1,0,0,-1],[dim,dim])

    sigma_x = reshape([0,1,1,0],[dim,dim])
 
    Id = reshape([1,0,0,1],[dim,dim])

    !call tensorial(sigma_z,sigma_z,dim,sig_zz)
    !call tensorial(Id,Id,dim,Id_2)
    !call tensorial(sigma_z,Id,dim,sig_z_Id)
    !call tensorial(Id,sigma_z,dim,Id_sig_z)

    call tensorial(sigma_x,Id,dim,sigma_x_Id)
    call tensorial(Id,sigma_x,dim,Id_sigma_x)

    call tensorial(sigma_x_Id,Id_2,dim*dim,s1_x)
    call tensorial(Id_sigma_x,Id_2,dim*dim,s2_x)
    call tensorial(Id_2,sigma_x_Id,dim*dim,s3_x)
    call tensorial(Id_2,Id_sigma_x,dim*dim,s4_x)



    ! Consulte o espa ç o de trabalho ideal .
    JOBZ = 'V'; UPLO = 'U'
    N = dim**4 ; LDA = dim**4 ; lwork = -1

    call print_matrix(Id_2,dim*2)

    print*, 'cheguei aqui'
    read(*,*)

    !call dsyev ( JOBZ , UPLO , N , A , LDA , W , WORK , LWORK , INFO )

    call dsyev (JOBZ, UPLO, s2_x, N , LDA , WORK, lwork, info)

    LWORK = MIN ( LWMAX , INT ( WORK (1) ) )

    



    ! Resolvendo o problema de autovalores .
    call dsyev ( JOBZ , UPLO , N , A , LDA , W , WORK , LWORK , INFO )





    ! Checando a converg ê ncia .
    if ( info > 0 ) then
    write (* ,*) 'O algoritmo falhou em encontrar os autovalores'
    stop
    end if



    ! Print autovalores .
    !call print_matrix ( 'eigenvalues' , 1 , n , w , 1)


    ! Print autovetores .
    !call print_matrix ( 'Eigenvectors' , n , n , A , LDA )
    end program lapack

    subroutine print_matrix(A,dim)
        implicit none
        integer, intent(in):: dim
        real*8, dimension(dim,dim), intent(in):: A
        integer:: i, j
        do i=1,dim
           write(*,20)(A(i,j),j=1,dim) !para numeros inteiros (4(I3))
        20   format (16(F3.0))
        enddo
  end subroutine


    
    subroutine tensorial(A,B,dim,E)
        implicit none
        integer, intent(in):: dim
        real*8, dimension(dim,dim), intent(in):: A
        real*8, dimension(dim,dim), intent(in):: B
        real*8, dimension(dim**2,dim**2), intent(out):: E
        integer:: i, j, k, l, alfa, beta
  
        do i = 1,dim
           do j = 1,dim
              do k = 1,dim
                 do l = 1,dim
  
                    alfa = dim*(i-1) + k
                    beta = dim*(j-1) + l
  
                    E(alfa,beta) = A(i,j) * B(k,l)
  
                 end do
              end do
           end do
        end do
  end subroutine
