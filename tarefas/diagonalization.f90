PROGRAM diagonalization
    IMPLICIT NONE
    integer:: dim
    real:: A(3,3),V(3,3),D(3,3)

    A = reshape([1,2,3,4,1,2,0,1,1],[dim,dim])
     
    call diagona(A,V,D,dim)

end PROGRAM

SUBROUTINE diagona(A, V, D, N)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    REAL, INTENT(INOUT) :: A(N,N)
    REAL, INTENT(INOUT) :: V(N,N)
    REAL, INTENT(INOUT) :: D(N,N)
    INTEGER :: i, j, k, p, q, ip, iq
    REAL :: c, s, t, tau, term1, term2

    ! Inicialização das matrizes D e V
    D = A
    V = 0.0
    DO i = 1, N
        V(i,i) = 1.0
    END DO

    ! Algoritmo de Jacobi
    DO k = 1, 1000*N**2 ! Número máximo de iterações permitidas
        ! Encontrar os índices p e q do elemento não-diagonal de maior módulo
        term1 = 0.0
        DO i = 1, N-1
            DO j = i+1, N
                IF (ABS(A(i,j)) > term1) THEN
                    term1 = ABS(A(i,j))
                    p = i
                    q = j
                END IF
            END DO
        END DO
        IF (term1 == 0.0) EXIT ! Matriz já é diagonal
        ! Cálculo dos ângulos de rotação
        tau = (A(q,q) - A(p,p))/(2.0*A(p,q))
        t = SIGN(1.0, tau)*1.0/(ABS(tau) + SQRT(1.0 + tau**2))
        c = 1.0/SQRT(1.0 + t**2)
        s = c*t
        ! Atualização das matrizes D e V
        term1 = A(p,p)
        term2 = A(q,q)
        A(p,p) = c**2*term1 + s**2*term2 - 2.0*c*s*A(p,q)
        A(q,q) = s**2*term1 + c**2*term2 + 2.0*c*s*A(p,q)
        A(p,q) = 0.0
        A(q,p) = 0.0
        DO i = 1, N
            IF (i /= p .AND. i /= q) THEN
                term1 = A(i,p)
                term2 = A(i,q)
                A(i,p) = c*term1 - s*term2
                A(p,i) = A(i,p)
                A(i,q) = s*term1 + c*term2
                A(q,i) = A(i,q)
            END IF
        END DO
        DO i = 1, N
            term1 = V(i,p)
            term2 = V(i,q)
            V(i,p) = c*term1 - s*term2
            V(i,q) = s*term1 + c*term2
        END DO
    END DO
    IF (k == 1000*N**2) STOP "Jacobi algorithm did not converge"
    D = A

END SUBROUTINE


! SUBROUTINE diag(A, V, D, N)

!     IMPLICIT NONE
!     INTEGER, INTENT(IN) :: N
!     REAL*8, INTENT(INOUT) :: A(N,N)
!     REAL*8, INTENT(INOUT) :: V(N,N)
!     REAL*8, INTENT(INOUT) :: D(N,N)
!     INTEGER :: i, j, k, l
!     REAL*8 :: c, s, t, tau, term1, term2
 
!     ! Inicialização das matrizes D e V
!     D = A
!     V = 0.0
!     DO i = 1, N
!         V(i,i) = 1.0
!     END DO
 
!     ! Algoritmo QR com shift de Wilkinson
!     DO l = 1, N-1
!         DO k = 1, 1000*N**2 ! Número máximo de iterações permitidas
!             ! Teste de convergência
!             term1 = ABS(D(l,l)) + ABS(D(l+1,l+1))
!             term2 = 0.0
!             DO i = l+1, N
!                 term2 = term2 + ABS(D(i,l))
!             END DO
!             IF (term2 <= 1e-12*term1) EXIT
!             ! Shift de Wilkinson
!             tau = (D(l,l) - D(l+1,l+1))/(2.0*D(l+1,l))
!             IF (tau >= 0.0) THEN
!                 t = 1.0/(tau + SQRT(1.0 + tau**2))
!             ELSE
!                 t = -1.0/(-tau + SQRT(1.0 + tau**2))
!             END IF
!             c = 1.0/SQRT(1.0 + t**2)
!             s = t*c
!             ! Atualização das matrizes D e V
!             DO i = l+1, N
!                 term1 = D(l,i)
!                 term2 = D(l+1,i)
!                 D(l,i) = c*term1 - s*term2
!                 D(l+1,i) = s*term1 + c*term2
!             END DO
!             DO i = 1, N
!                 term1 = V(i,l)
!                 term2 = V(i,l+1)
!                 V(i,l) = c*term1 - s*term2
!                 V(i,l+1) = s*term1 + c*term2
!             END DO
!             DO j = l+2, N
!                 term1 = D(j,l)
!                 term2 = D(j,l+1)
!                 D(j,l) = c*term1 - s*term2
!                 D(j,l+1) = s*term1 + c*term2
!             END DO
!         END DO
!         IF (k == 1000*N**2) STOP "QR algorithm did not converge"
!     END DO
 
 
!     call print_matrix(A,N)
!     read(*,*)
 
!     call print_matrix(V,N)
!     read(*,*)
 
!     call print_matrix(D,N)
!     read(*,*)
 
!  END SUBROUTINE


