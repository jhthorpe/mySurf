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

END MODULE linal
