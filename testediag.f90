program testeDiagonalization
implicit none
integer , parameter :: dim = 2 , db = 8
real ( KIND = db ) , dimension ( dim**4 , dim**4 ) :: A,Ub,mat
! !!! Para diagonalizar !!!!!
character ( len =1) :: JOBZ , UPLO
integer :: N , LDA , INFO , lwork,i,j
integer , parameter :: LWMAX = 1000
real ( KIND = db ) , dimension ( dim**4 ) :: W
real ( kind = db ) , dimension ( LWMAX ) :: WORK
real*8, dimension(dim*dim,dim*dim):: Id_2, Sigma_x_Id
real*8, dimension(2,2):: Id, sigma_x

! A = reshape ([3. , -2. ,4. , -2. ,6. &
!              ,2. , 4. , 2. , 3., 1.&
!              ,1. , 2. , 1. , -1., 1.&
!              ,2. , 4. , 2. , 3., 1.&
!              ,2. , 4. , 2. , 3., 1.] , [ dim , dim ])

             sigma_x = reshape([0,1,1,0],[2,2])

             Id = reshape([1,0,0,1],[2,2])

             call tensorial(Id,Id,dim,Id_2)
             call tensorial(sigma_x,Id,dim,sigma_x_Id)
             
             call tensorial(sigma_x_Id,Id_2,dim*dim,A)

! Consulte o espa ç o de trabalho ideal .
JOBZ = 'V'; UPLO = 'U'
N = dim**4 ; LDA = dim**4 ; lwork = -1

call dsyev ( JOBZ , UPLO , N , A , LDA , W , WORK , LWORK , INFO )

LWORK = MIN ( LWMAX , INT ( WORK (1) ) )

! Resolvendo o problema de autovalores .
call dsyev ( JOBZ , UPLO , N , A , LDA , W , WORK , LWORK , INFO )

! Checando a converg ê ncia .
if ( info > 0 ) then
write (* ,*) 'O algoritmo falhou em encontrar os autovalores'
stop
end if

! Print autovalores .
call print_matrix ( 'eigenvalues' , n , 1 , w , 1)

read(*,*)



! Print autovetores .
call print_matrix ( 'Eigenvectors' , n , n , A , LDA )

Ub = transpose(A)

mat= matmul(A,Ub)
print*, '-------------------------------------'
do i = 1, dim**4

        write (* ,10) ( mat (i , j ) , j =1 , dim**4 )
        10 format (16( F6 .2) )

end do
print*, '-------------------------------------'
do i = 1, dim**4

    write (* ,20) ( A (i , j ) , j =1 , dim**4 )
    20 format (16( F6 .2) )

end do

end program testeDiagonalization

subroutine print_matrix ( desc , m , n , a , LDA )
implicit none
integer , parameter :: db = 8
character ( len =*) , intent ( in ) :: desc
integer , intent ( in ) :: m ,n , lda
real ( kind = db ) , dimension ( LDA ,*) , intent ( in ) :: A
integer :: i , j

write (* ,*) desc

do i = 1 , m
write (* ,20) ( A (i , j ) , j =1 , n )
end do
20 format (5( F6 .2) )
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