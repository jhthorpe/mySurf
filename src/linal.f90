!linear algebra routines
MODULE linal

CONTAINS

!---------------------------------------------------------------------
! linal_sym_diag
!    - diagonalizes a symmetric matrix
!---------------------------------------------------------------------
! N             : int, number of dimensions
! A             : 2D real*8, matrix to diagonalize
! l             : 1D real*8, array of eigenvalues
! error         : int, exit code

SUBROUTINE linal_sym_diag(N,A,l,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:N-1,0:N-1), INTENT(INOUT) :: A
  REAL(KIND=8), DIMENSION(0:N-1), INTENT(INOUT) :: l
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  INTEGER :: i,j
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: WORK
  CHARACTER(LEN=1) :: JOBZ, UPLO,RNGE
  INTEGER :: LDA,LWORK,IL,UL

  error = 0

  JOBZ = 'V'
  RNGE='I'
  UPLO = 'L'
  LWORK = -1
  LDA = N
  IL = 0
  UL = N
  ALLOCATE(WORK(0:1))

  CALL DSYEV(JOBZ,UPLO,N,A(0:N-1,0:N-1),LDA,l(0:N-1),WORK(0:1),LWORK,error)

  LWORK = CEILING(WORK(0))
  DEALLOCATE(WORK)
  ALLOCATE(WORK(0:LWORK-1))

  CALL DSYEV(JOBZ,UPLO,N,A(0:N-1,0:N-1),LDA,l(0:N-1),WORK(0:LWORK-1),LWORK,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "linal_sym_diag : error from DSYEV",error
    DEALLOCATE(WORK)
    RETURN
  END IF

  !Standardize the eigenvectors
  DO j=0,N-1
    i = MAXLOC(ABS(A(0:N-1,j)),1)-1
    IF (A(i,j) .LT. 0) A(0:N-1,j) = -1.0D0*A(0:N-1,j)
  END DO

  DEALLOCATE(WORK)

END SUBROUTINE linal_sym_diag

!---------------------------------------------------------------------
! linal_invert
!    - invert a matrix
!---------------------------------------------------------------------
! N             : int, size of square matrix
! A             : 2D real*8, matrix to invert
! Ai            : 2D real*8, inverted matrix
! error         : int, error code

SUBROUTINE linal_invert(N,A,Ai,error)
  IMPLICIT NONE

  REAL(KIND=8), DIMENSION(0:N-1,0:N-1), INTENT(INOUT) :: Ai
  REAL(KIND=8), DIMENSION(0:N-1,0:N-1), INTENT(IN) :: A
  INTEGER, INTENT(INOUT) :: error
  INTEGER, INTENT(IN) :: N

  REAL(KIND=8), DIMENSION(0:N-1) :: WORK
  INTEGER, DIMENSION(0:N-1) :: IPIV

  error = 0
  Ai = A

  CALL DGETRF(N,N,Ai(0:N-1,0:N-1),N,IPIV(0:N-1),error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "linal_invert : error from DGETRF", error
    RETURN
  END IF

  CALL DGETRI(N,Ai(0:N-1,0:N-1),N,IPIV(0:N-1),WORK(0:N-1),N,error)
  IF (error .NE. 0) THEN
    WRITE(*,*) "linal_invert : error from DGETRI", error
    RETURN
  END IF

END SUBROUTINE linal_invert

!---------------------------------------------------------------------
! linal_matvec
!    - performs matrix vector multiplication
!      y = A.x
!---------------------------------------------------------------------
! n             : int, size of square matrix and vector
! A             : 2D real*8, matrix
! x             : 1D real*8, input vector 
! y             : 1D real*8, output vector

SUBROUTINE linal_matvec(n,A,x,y)
  IMPLICIT NONE
  REAL(KIND=8), DIMENSION(0:n-1,0:n-1), INTENT(IN) :: A
  REAL(KIND=8), DIMENSION(0:n-1), INTENT(INOUT) :: y
  REAL(KIND=8), DIMENSION(0:n-1), INTENT(IN) :: x
  INTEGER, INTENT(IN) :: n

  CHARACTER(LEN=1) :: TRANS
  REAL(KIND=8) :: ALPHA,BETA
  INTEGER :: INCX,INCY

  TRANS = 'N'
  ALPHA = 1.0D0
  BETA = 0.0D0
  INCX = 1
  INCY = 1

  CALL DGEMV(TRANS,n,n,ALPHA,A(0:n-1,0:n-1),n,x(0:n-1),INCX,BETA,y(0:n-1),INCY) 

END SUBROUTINE linal_matvec

!---------------------------------------------------------------------
! linal_diag_inplace
!       - uses Jacobi Transforms to block-diagonalize a sym. real matrix 
!         in the form specified by `blocks`
!       - very poorly coded, but it will do for small matrices
!       - based off algorithms described in Numerical Recipies
!       - uses lower triangular form
!       - standardizes all eigenvectors to have their largest value
!         as positive
!---------------------------------------------------------------------
! Variables
! n             : int, size of matrix
! A             : 2D real*8, real, symmetric matrix to be block-diag
! V             : 2D real*8, transform matrix 
! w             : 1D real*8, eigenvalues
! info          : int, information about how calculation went

SUBROUTINE linal_diag_inplace(n,A,w,info)
  IMPLICIT NONE
  !inout
  REAL(KIND=8), DIMENSION(0:n-1,0:n-1), INTENT(INOUT) :: A
  REAL(KIND=8), DIMENSION(0:n-1),INTENT(INOUT) :: w
  INTEGER, INTENT(INOUT) :: info
  INTEGER, INTENT(IN) :: n
  !internal
  REAL(KIND=8), DIMENSION(0:n-1,0:n-1) :: V
  REAL(KIND=8) :: conv,tol
  INTEGER :: i,j

  info = 0
  tol = 1.0D-16

  !initialize
  V = 0
  DO i=0,n-1
    V(i,i) = 1.0D0
  END DO

  DO i=0,49 !max of 50 sweeps
    !check convergence
    CALL get_conv(conv,A,n)
    IF (conv .LT. tol) EXIT
    !otherwise, do Jacobi sweep
    CALL jacobi_sweep(n,A,V,i)
  END DO

  !Standardize the eigenvectors
  DO j=0,n-1
    i = MAXLOC(ABS(V(0:n-1,j)),1)-1
    IF (V(i,j) .LT. 0) V(0:n-1,j) = -1.0D0*V(0:n-1,j)
  END DO

  DO i=0,n-1
    w(i) = A(i,i)
  END DO
  
  A = V

  IF (conv .LT. tol) THEN
    info = 0
  ELSE
    info = 1
  END IF
END SUBROUTINE linal_diag_inplace

!---------------------------------------------------------------------
!       get_conv
!               James H. Thorpe
!       - checks convergence of block-diagonal matrix
!---------------------------------------------------------------------
! conv          : real*8, convergence
! A             : 2D real*8, matrix
! n             : int, rank

SUBROUTINE get_conv(conv,A,n)
  IMPLICIT NONE   !inout
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: A
  REAL(KIND=8), INTENT(INOUT) :: conv
  INTEGER, INTENT(IN) :: n
  !internal
  INTEGER :: i,j

  conv = 0.0D0

  DO j=0,n-1
    DO i=j+1,n-1
      conv = conv + 2*(A(i,j)**2.0D0)
      IF (ABS(A(i,j)) .GT. 1.0D-15) conv = conv + 100 
    END DO
  END DO

END SUBROUTINE get_conv

!---------------------------------------------------------------------
!       jacobi_sweep
!               James H. Thorpe
!       -performs jacobi sweep of matrix to diagonalize elements in 
!        blocks, in lower triangular form
!       - uses code from NR 77,adapted for lower triangular
!---------------------------------------------------------------------
! n             : int, rank of matrix
! A             : 2D real*8, matrix to be reduced
! V             : 2D real*8, eigenvectors 
! w1,w2         : 1D real*8, working vectors
! p,q           : int, indicies
! c,s,t         : real*8, same as in Num. Rec.
! th            : real*8, theta(θ) in Num. Rec. 
! ta            : real*8, tau(τ) in Num. Rec.
! h             : real*8, differences between A(q,q) and A(p,p)
! itr           : int, iteration of sweep

SUBROUTINE jacobi_sweep(n,A,V,itr)
  IMPLICIT NONE
  !inout
  REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: A,V
  INTEGER, INTENT(IN) :: n,itr
  !internal
  REAL(KIND=8) :: c,s,t,o,th,ta,tol,h,g,f
  INTEGER :: p,q,i,j,k,r

  tol = 1.0D-16

  !sweep through indices
  !DO p=0,n-2     !original
  !  DO q=p+1,n-1 !original
  DO q=0,n-2
    DO p=q+1,n-1
      g = 100.0D0*ABS(A(p,q))       !numerical stuff

      IF ( (itr .GT. 4) .AND. (ABS(A(p,p))+g .EQ. ABS(A(p,p))) &
           .AND. (ABS(A(q,q))+g .EQ. ABS(A(q,q))) ) THEN
        A(p,q) = 0.0D0
      ELSE IF (ABS(A(p,q)) .GT. tol) THEN !element is not zero, eliminate
        h = A(q,q) - A(p,p)
        IF ((ABS(h)+g) .EQ. ABS(h)) THEN !machine overflow
          t = A(p,q)/h
        ELSE !normal case
          th = 0.5D0*h/A(p,q)
          t = 1.0D0/(ABS(th)+SQRT(1.0D0 + th**2.0D0))
          IF (th .LT. 0) t = -t
        END IF !machine overflow

        !get parameters
        c = 1.0D0/SQRT(1.0D0 + t**2.0D0)
        s = t*c
        ta = s/(1.0D0 + c)
        f = t*A(p,q)

        !change A matrix 
        DO r=0,p-1 !1 <= r < p         
          g = A(r,p)
          h = A(r,q)
          A(r,p) =  g - s*(h + ta*g)
          A(r,q) =  h + s*(g - ta*h)
        END DO

        A(p,p) = A(p,p) - f !r = p, points p,p and p,q 
        A(p,q) = 0.0D0

        DO r=p+1,q-1 !p < r < q
          g = A(p,r)
          h = A(r,q)
          A(p,r) =  g - s*(h + ta*g)
          A(r,q) =  h + s*(g - ta*h)
        END DO

        !A(q,p) = 0.0D0 !points q,q and q,p
        A(q,q) = A(q,q) + f

        DO r=q+1,n-1 !q < j <= n
          g = A(p,r)
          h = A(q,r)
          A(p,r) =  g - s*(h + ta*g)
          A(q,r) =  h + s*(g - ta*h)
        END DO

        !change V matrix
        DO r=0,n-1
          g = V(r,p)
          h = V(r,q)
          V(r,p) = g - s*(h + ta*g)
          V(r,q) = h + s*(g - ta*h)
        END DO

      END IF !nonzero elements

    END DO
  END DO

END SUBROUTINE jacobi_sweep
!---------------------------------------------------------------------



END MODULE linal
